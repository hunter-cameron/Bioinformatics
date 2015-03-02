#!/usr/bin/env python
from __future__ import print_function

import argparse
import subprocess
import os
import sys
import json
from mypyli import taxtree       # mypyli built using python3 so this may not work

def main(fasta, gbk, out, processors, tree=None):
    # change default printing to err for working steps
    stdout = sys.stdout
    sys.stdout = sys.stderr

    # check the output directory
    if not os.path.isdir(out):
        print("Output path {} is not a directory".format(out))
        sys.exit(1)

    # get lists of the files to be used (by reading directory if necessary)
    
    fasta_files = []
    if os.path.isdir(fasta):
        fasta_files = process_dir(fasta, "fasta")
    elif os.path.isfile(fasta):
        fasta_files.append(fasta)
    else:
        print("Path {} is neither a file nor directory.".format(fasta))
        sys.exit(1)

    # make single element tuples if gbk isn't specified
    if not gbk:
        genomes = [(fasta,) for fasta in fasta_files]
    else:
        gbk_files = []
        if os.path.isdir(gbk):
            gbk_files = process_dir(gbk, "gbk")
        elif os.path.isfile(gbk):
            gbk_files.append(gbk)
        else:
            print("Path {} is neither a file nor directory.".format(gbk))
            sys.exit(1)
    
        # link each fasta with its genbank
        genomes = link_fasta_w_gbk(fasta_files, gbk_files)
    
    output_paths = []
    for genome in genomes:
        print("Processing genome: {}".format(genome[0]))
            
        # make output directory
        out_path = out + "/" + get_basename(genome[0]) + "/"
        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        # if this file exists, the run for that isolate has already been done
        if os.path.isfile("{}/storage/bin_stats_ext.tsv".format(out_path)):
            print("    Genome already complete. Skipping.")
            output_paths.append(out_path)
            continue

        # link the fasta to the folder if it isn't already
        test_fasta = process_dir(out_path, "fasta")
        subprocess.call(["ln", "-sf", genome[0], out_path])
        
        if gbk:
            # begin trying to run CheckM, starting at the lowest taxonomic level
            for level in reversed(get_tax_from_gbk(genome[1], tree)):
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
                print("Processing genome: {} failed.".format(genome[0]), file=stdout)

        else:
            print("Running checkm using lineage_wf.")
            rtn_code = run_checkm(rank_taxon=None, fasta_bin=out_path, processors=processors)
            if not rtn_code:
                print("  CheckM ran without error.", end="\n\n\n")
                output_paths.append(out_path)
            else:
                print("  CheckM encountered an error. Genome cannot be processed")
                print("Processing genome: {} failed.".format(genome[0]), file=stdout)
                

    # merge all results into a single file 
    out_ext = out + "/bin_stats_ext.tsv"
    print("\n\nConcatenating all the output results to {}".format(out_ext), file=stdout)
    # add the file to the output root
    output_paths = [path + "/storage/bin_stats_ext.tsv" for path in output_paths]
    merge_ext_files(out_ext, output_paths)

    print("\n\nFinished.", file=stdout)

def process_dir(directory, ext):
    """ Returns a list of files with the given ext from a directory """
    local_files = os.listdir(directory)
    return [directory + "/" + file for file in local_files if file.endswith(".{}".format(ext))]

def link_fasta_w_gbk(fasta_files, gbk_files):
    """ Matches the fasta and gbk files based on name similarity and returns a list of tuples of matched files"""
    files = []
    for fasta in fasta_files:

        f_base = get_basename(fasta)
        for gbk in gbk_files:
            g_base = get_basename(gbk)
            if f_base == g_base:
                files.append((fasta, gbk))
                break
        else:
            print("No gbk match found for fasta: {}".format(fasta))
            print("Skipping...")

    return files

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
                    try:
                        tax_obj = tree.lookup_single_tax(tax_string)
                    except LookupError:
                        print(" Lookup failed. Running CheckM will fail.")
                        return [("none", "none")]
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
    parser.add_argument("-gbk", type=str, help="single gbk file or folder of gbk files (ext: .gbk) to use for parsing taxonomy")
    parser.add_argument("-out", type=str, help="directory prefix for results (defaults to current directory)", default=os.getcwd())
    parser.add_argument("-n", type=int, help="number of processors", required=True)
    parser.add_argument("-tree", type=str, help="taxtree pickled tree to load to lookup taxonomies of genomes without the ORGANISM line in the GenBank")
    args = parser.parse_args()

    if args.tree:
        tree = taxtree.TaxTree.load_tree(args.tree)
    else:
        tree=None

    main(args.fasta, args.gbk, args.out, args.n, tree)
