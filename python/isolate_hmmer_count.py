#!/usr/bin/env python

import sys
import argparse
import subprocess
import os
from pprint import pprint as pprint

def main():
    """
    The main body of the hmmer_count program. Basically reads in arguments and calls functions.
    """
    
    parser = argparse.ArgumentParser(description="Generates a table of counts for the number of times each gene is found in the genome.")

    parser.add_argument("genome", help="FAA of genes in the input genome.")
    parser.add_argument("-db", help="Path to the HMM database", default="/proj/dangl_lab/hunter/wordle/jgi_universal_db/universal")
    parser.add_argument("-hmmscan", help="Path to hmmscan executable if it is not in system path", default="hmmscan")
    parser.add_argument("-e", help="E-value cut off for profile matches (default = 1e-10)", type=int, default=1e-10)

    args = parser.parse_args()

    dictionary = generate_empty_gene_dictionary(args.db)
    hmmscan_out = run_hmmscan(args.genome, args.db, args.hmmscan, args.e)
    counts = parse_hmmscan_tbl(hmmscan_out, dictionary)

    out_file = os.path.basename(args.genome)
    out_file = os.path.splitext(out_file)[0]
    print_dictionary_table(counts, "{}_hmmcount.txt".format(out_file))

def run_hmmscan(query, db, hmmscan, min_e):
    """
    Runs hmmscan and returns the path to the output table file.
    """

    #need to check for completion somehow perhaps?
    tblout = "{query}-hmmscan_out.txt".format(query=query)
    p = subprocess.Popen([hmmscan, "--tblout", tblout, "-E", str(min_e), db, query], stdout=subprocess.PIPE)
    out, err = p.communicate()

    return tblout
    
def parse_hmmscan_tbl(table, dictionary):
    """
    Parses the hmmscan table and keeps a count using the dictionary.
    """
    
    with open(table, 'r') as IN:
        for line in IN:
            if line.startswith("#"):
                continue
            elements = line.split()
            gene = elements[0]
            dictionary[gene] += 1
            
    return dictionary
            
def generate_empty_gene_dictionary(file):
    """ Generates a dictionary with all the genes in a given HMM file set to 0
        
"""
    dictionary = {}
    with open(file, 'r') as IN:
        for line in IN:
            if line.startswith("NAME"):
                _, _, name = line.partition("  ")      #split string on two spaces, return before, sep, after (all we want is the after)                
                dictionary[name.strip()] = 0

    return dictionary
                            
def print_dictionary_table(dictionary, file):
    """ Prints a dictionary as a tab-delimited file """

    with open(file, 'w') as OUT:
        OUT.write("Gene\tCounts\n")
        keys = sorted(dictionary.keys())
        for key in keys:
            OUT.write("{}\t{}\n".format(key, dictionary[key]))

            



if __name__ == "__main__":
    main()

