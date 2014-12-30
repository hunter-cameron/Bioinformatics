#!/usr/bin/env python
from __future__ import print_function

import argparse
import subprocess
import os
import sys

def main(fasta, gbk, out, processors):
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
        fasta.files.append(fasta)
    else:
        print("Path {} is neither a file nor directory.".format(fasta))
        sys.exit(1)

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

    for genome in genomes:
        print("Processing genome: {}".format(genome[0]))
            
        # make output directory
        out_path = out + "/" + get_basename(genome[0]) + "/"
        if not os.path.isdir(out_path):
            os.mkdir(out_path)

        # link the fasta to the folder if it isn't already
        test_fasta = process_dir(out_path, "fasta")
        if not get_basename(genome[0]) in [get_basename(test) for test in test_fasta]:
            subprocess.call(["ln", "-s", genome[0], out_path])

        # begin trying to run CheckM, starting at the lowest taxonomic level
        for level in reversed(get_tax_from_gbk(genome[1])):
            print("  Running checkm with taxonomy: {} {}.".format(level[0], level[1]))
            rtn_code = run_checkm(level, out_path, processors)
            # exit loop if no error use next tax level if error
            if not rtn_code:
                print("  CheckM ran without error.", end="\n\n\n")
                break
            else:
                print("  CheckM encountered an error. Running at next level.", end="\n\n")

        else:
            print("  Already at highest level. Genome cannot be processed")
            print("Processing genome: {} failed.", file=stdout)

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

def get_tax_from_gbk(gbk_file):
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

    rank_tax = []
    for rank, tax in zip(["domain", "phylum", "class", "order", "family", "genus", "species"], tax_string.split("; ")):
        rank_tax.append((rank, tax))

    return rank_tax

def run_checkm(rank_taxon, fasta_bin, processors=1):

    """ Tries to run checkm tax_wf and returns the error code generated. """

    return subprocess.call(["checkm", "taxonomy_wf",
                            "-x", "fasta",
                            "-t", str(processors),
                            rank_taxon[0], rank_taxon[1],
                            fasta_bin, fasta_bin],
                            stdout=open(fasta_bin + "stdout_checkm", 'w'),
                            stderr=open(fasta_bin + "stderr_checkm", 'w'))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Runs CheckM in taxonomy_wf mode for the lowest taxonomy possible based on a gbk file")

    parser.add_argument("-fasta", type=str, help="single fasta file or folder of fasta files (ext: .fasta)", required=True)
    parser.add_argument("-gbk", type=str, help="single gbk file or folder of gbk files (ext: .gbk) to use for parsing taxonomy", required=True)
    parser.add_argument("-out", type=str, help="directory prefix for results (defaults to current directory)", default=os.getcwd())
    parser.add_argument("-n", type=int, help="number of processors", required=True)
    args = parser.parse_args()

    main(args.fasta, args.gbk, args.out, args.n)

