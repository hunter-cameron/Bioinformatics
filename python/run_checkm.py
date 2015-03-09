#!/usr/bin/env python
from __future__ import print_function

import argparse
import subprocess
import os
import sys
import json
from mypyli import taxtree       # mypyli built using python3 so this may not work

class TaxonomyError(Exception):
    def __init__(self, message, target=None):
        self.message = message
        self.target = target


def main(fasta_f, tax_f, taxonly, out, processors, tree=None):
    # change default printing to err for working steps
    stdout = sys.stdout
    sys.stdout = sys.stderr

    # check the output directory
    if not os.path.isdir(out):
        print("Output path {} is not a directory".format(out))
        sys.exit(1)

 
    # get lists of the files to be used (by reading directory if necessary)
    fasta_files = read_file_or_dir(fasta_f, ext="fasta")

    # make single element tuples if gbk isn't specified
    if not tax_f:
        fasta_tax = [(fasta, None) for fasta in fasta_files]
    else:
        if taxonly:
            fasta_tax = link_fasta_w_tax(fasta_files, tax_f, tree=tree)
        else:
            gbk_files = read_file_or_dir(tax_f, ext="gbk")
    
            # link each fasta with its genbank derived taxonomy
            fasta_tax = link_fasta_w_gbk(fasta_files, gbk_files, tree=tree)
   

    output_paths = []
    for fasta, taxonomy in fasta_tax:
        print("Processing genome: {}".format(fasta))
        
        # make output directory
        out_path = out + "/" + get_basename(fasta) + "/"
        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        # if this file exists, the run for that isolate has already been done
        if os.path.isfile("{}/storage/bin_stats_ext.tsv".format(out_path)):
            print("    Genome already complete. Skipping.")
            output_paths.append(out_path)
            continue

        # link the fasta to the folder if it isn't already -- this will overwrite files by the same name
        subprocess.call(["ln", "-sf", fasta, out_path])
        
        if tax_f:
            # begin trying to run CheckM, starting at the lowest taxonomic level
            for level in reversed(taxonomy):
                print("  Running checkm with taxonomy: {} {}.".format(level[0], level[1]))
                rtn_code = run_checkm(level, out_path, processors)
                # exit loop if no error use next tax level if error
                if not rtn_code:
                    print("  CheckM ran without error.", end="\n\n\n")
                    output_paths.append(out_path)
                    break
                else:
                    print("  CheckM encountered an error. Running at next level.", end="\n\n")
    
            else:
                print("  Already at highest level. Genome cannot be processed")
                print("Processing genome: {} failed.".format(fasta), file=stdout)

        # if no tax info, run using lineage_wf
        else:
            print("Running checkm using lineage_wf.")
            rtn_code = run_checkm(rank_taxon=None, fasta_bin=out_path, processors=processors)
            if not rtn_code:
                print("  CheckM ran without error.", end="\n\n\n")
                output_paths.append(out_path)
            else:
                print("  CheckM encountered an error. Genome cannot be processed")
                print("Processing genome: {} failed.".format(fasta), file=stdout)
                

    # merge all results into a single file 
    out_ext = out + "/bin_stats_ext.tsv"
    print("\n\nConcatenating all the output results to {}".format(out_ext), file=stdout)
    # add the file to the output root
    output_paths = [path + "/storage/bin_stats_ext.tsv" for path in output_paths]
    merge_ext_files(out_ext, output_paths)

    print("\n\nFinished.", file=stdout)

def read_file_or_dir(path, ext=None):
    """ 
    Takes a directory or file and if it is a directory, returns
    a list of files in that directory with the specified ext.
    If it is a file, returns the filename in list context.
    """

    files = []
    if os.path.isdir(path):
        if not ext:
            raise ValueError("Path {} is a directory and no extension was provided")
        files = os.listdir(path)
        files = [path + "/" + file for file in files if file.endswith(".{}".format(ext))]
        if files:
            return files
        else:
            raise ValueError("No files with ext {} found at path {}".format(ext, path)) 

    elif os.path.isfile(path):
        return [path]
    else:
        raise IOError("Path {} is neither a file nor directory.".format(path))



# functions to link the fasta with the gbk and get the taxonomy
def link_fasta_w_gbk(fasta_files, gbk_files, tree=None):
    """ 
    Matches the fasta and gbk files based on name similarity and parses the taxonomy from the gbk files
    Returns a list of tuples, (fasta_name, [taxonomy])

    Raises: TaxonomyError
    """
    fasta_tax = []
    for fasta in fasta_files:

        f_base = get_basename(fasta)
        for gbk in gbk_files:
            g_base = get_basename(gbk)
            if f_base == g_base:
                try:
                    fasta_tax.append((fasta, get_tax_from_gbk(gbk, tree=tree)))
                except:     # I just want to catch all errors with taxonomy parsing
                    print("  Taxonomy could not be found in gbk: {}".format(gbk))
                    print("  Skipping genome...")
                break
        else:
            print("  No gbk match found for fasta: {}".format(fasta))
            print("  Skipping genome...")

    return fasta_tax

def get_basename(path):
    basename = os.path.basename(path)
    return os.path.splitext(basename)[0]

def get_tax_from_gbk(gbk_file, tree=None):
    """ Returns a taxstring from the info in a gbk file """
    tax_lines = []
    append = False
    with open(gbk_file, 'r') as IN:
        for line in IN:
            if line.startswith("  ORGANISM"):
                append = True
            if line.startswith("FEATURES"):
                break

            if append:
                tax_lines.append(line)

    # remove the first line because it's junk and remove leading and trailing spaces from lines
    clean_lines = [line.strip() for line in tax_lines[1:]]

    # concatenate the lines and remove that period they put at the end
    tax_string = " ".join(clean_lines)
    tax_string = tax_string.replace(".", "")

    # if it didn't have the ORGANISM line, look for the organism in source
    if not tax_string:
        if tree is None:
            print("  Organism tag not found and tree is not available. Running CheckM will fail.")
            return [("none", "none")]   # this should throw an error not return
        print("  ORGANISM tag not found, looking for '/organism' in Source")
        with open(gbk_file, 'r') as IN:
            for line in IN:
                line = line.strip()
                if line.startswith("/organism="):
                    # this line just gets the first word after a double quote
                    tax_string = line.split('"')[1].split(" ")[0]
                    print("  Organism found: {}.".format(tax_string))
                    tax_obj = tree.lookup_single_tax(tax_string)
                    rank_tax = []
                    # look up taxonomy by rank, break at first undefined
                    for rank in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
                        
                        tax = tax_obj.get_tax_at_rank(rank, null="")
                        if tax:
                            rank_tax.append((rank, tax))
                        else:
                            break
                    print("  Found lineage in taxtree {}".format(tax_obj.get_tax_string()))
                    return rank_tax
    else:
        rank_tax = []
        for rank, tax in zip(["domain", "phylum", "class", "order", "family", "genus", "species"], tax_string.split("; ")):
            rank_tax.append((rank, tax))

        return rank_tax

# function to link fasta with tab delimited tax
def link_fasta_w_tax(fasta_files, tax_f, tree=None):
    name_to_tax = {}
    with open(tax_f, 'r') as IN:
        for line in IN:
            fasta, taxstring = line[:-1].split("\t")
            name_to_tax[fasta] = taxstring

    fasta_tax = []
    for fasta in fasta_files:
        try:
            taxstring = name_to_tax[get_basename(fasta)]
        except KeyError:
            print("  No taxonomy match found for fasta: {}".format(fasta))
            print("  Skipping genome...")
            continue

        if tree:
            rank_tax = []
            node = tree.lookup_taxstring(taxstring)
            for rank in ["domain", "phylum", "class", "order", "family", "genus", "species"]:
                        
                tax = node.get_tax_at_rank(rank, null="")
                if tax:
                    rank_tax.append((rank, tax))
                else:
                    break
            #print("  Found lineage in taxtree {}".format(node.get_tax_string()))
            fasta_tax.append((fasta, rank_tax))
        else:
            raise NotImplementedError("Currently a taxtree is required to parse tax files")

    return fasta_tax

def merge_ext_files(output_f, ext_paths):
    """ Merges any number of checkm ext files to a single path"""
    
    full_dict = {}
    for path in ext_paths:
        with open(path, 'r') as IN:
            print("Merging {}...".format(path))
            # json uses double quotes while the CheckM folks used singles.
            full_dict.update(json.loads(IN.read().replace("'", '"')))
            #print(full_dict)

    with open(output_f, 'w') as OUT:
        json.dump(full_dict, OUT)

def run_checkm(rank_taxon, fasta_bin, processors=1):

    """ Tries to run checkm tax_wf and returns the error code generated. """

    if rank_taxon:

        return subprocess.call(["checkm", "taxonomy_wf",
                            "-x", "fasta",
                            "-t", str(processors),
                            rank_taxon[0], rank_taxon[1],
                            fasta_bin, fasta_bin],
                            stdout=open(fasta_bin + "stdout_checkm", 'w'),
                            stderr=open(fasta_bin + "stderr_checkm", 'w'))

    else:
        return subprocess.call(["checkm", "lineage_wf",
                            "-x", "fasta",
                            "-t", str(processors),
                            fasta_bin, fasta_bin],
                            stdout=open(fasta_bin + "stdout_checkm", 'w'),
                            stderr=open(fasta_bin + "stderr_checkm", 'w'))



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Runs CheckM in taxonomy_wf mode for the lowest taxonomy possible based on a gbk file. Alternatively, if the -gbk option is left blank, this program will run CheckM in lineage_wf. Both workflows check first if the isolate has already been completed (and skip it is so) and yield a concatenated file of all data at the end.")

    parser.add_argument("-fasta", type=str, help="single fasta file or folder of fasta files (ext: .fasta)", required=True)
    parser.add_argument("-tax", type=str, help="single gbk file or folder of gbk files (ext: .gbk) to use for parsing taxonomy")
    parser.add_argument("-taxonly", help="flag to specify that the tax argument isn't a gbk file but a tab delim file of fasta and taxstring", action='store_true')
    parser.add_argument("-out", type=str, help="directory prefix for results (defaults to current directory)", default=os.getcwd())
    parser.add_argument("-n", type=int, help="number of processors", required=True)
    parser.add_argument("-tree", type=str, help="taxtree pickled tree to load to lookup taxonomies of genomes without the ORGANISM line in the GenBank")
    args = parser.parse_args()

    if args.tree:
        tree = taxtree.TaxTree.load_tree(args.tree)
    else:
        tree=None


    main(args.fasta, args.tax, args.taxonly, args.out, args.n, tree)
