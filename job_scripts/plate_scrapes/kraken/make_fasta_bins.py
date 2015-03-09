
import sys
import argparse
import os
import re
from Bio import SeqIO
import numpy as npy
from mypyli import taxtree, kraken


def get_lengths(fasta_f):
    """ Returns as dict with the sequences and their lengths """
    contig_data = {}

    with open(fasta_f, 'rU') as IN:
        for record in SeqIO.parse(IN, "fasta"):
            length = len(record.seq)

            contig_data[record.id] = {'length': length}

    return contig_data
    
def get_taxonomic_classifications(contig_data, kraken_f, tree):
    """ Parses (currently) the kraken file from contigs and updates the hash with contig name: TaxNode, length
    """

    for krec in kraken.KrakenRecord.parse_kraken_file(kraken_f, iterate=True):
        if krec.classified:
            tax_node = tree.lookup_taxid(krec.taxid)
        else:
            tax_node = "unassigned"

        contig_data[krec.name].update({'taxonomy': tax_node})
   
def bin_sequences(contig_data, rank="genus", imerge=False):

    bins = {}
    for contig, data in contig_data.items():
        taxstring = data['taxonomy']

        if taxstring != "unassigned":
            taxstring = taxstring.get_tax_string(trim_to=rank)
        
        bins[taxstring] = bins.get(taxstring, []) + [contig]



    while imerge:
        print("   Running intelligerge...", file=sys.stderr)
        imerge = intellimerge(bins)


    print("\n\nBins\t#seqs\ttotal_bp")
    for bin in sorted(bins):
        bp = 0
        for header in bins[bin]:
            bp += contig_data[header]['length']
        print("\t".join([bin, str(len(bins[bin])), str(bp)]))

    return bins


def intellimerge(bins, minimum=10000):
    """
    Need to control for bins that get moved up from kingdom to nothing if thats possible"""
    
    modified = 0
    num_bins = len(bins)
    print(("num_bins", num_bins))
    for bin in sorted(bins.copy(), reverse=True):
        bp = 0
        for header in bins[bin]:
            bp += contig_data[header]['length']

        if bp < minimum:
            tax = bin.split("; ")
            
            if len(tax) == 1:
                continue
            else:
                new_bin = "; ".join(tax[:-1])
            old_bin = bins.pop(bin)
            #print(old_bin)
            bins[new_bin] = bins.get(new_bin, []) + old_bin
        

    
    return not num_bins == len(bins)


def write_fastas(bins, fasta_f, out):

    if not os.path.isdir(out):
        os.makedirs(out)

    # make lookup hash and initialize files
    lookup = {}
    for bin in bins:
        for contig in bins[bin]:
            lookup[contig] = bin

        #open("{}/{}.fasta".format(out, bin), 'w')
        open("{}/{}.fasta".format(out, bin.replace("; ", "-")), 'w')

    with open(fasta_f, 'rU') as IN:
        for record in SeqIO.parse(IN, "fasta"):
            bin = lookup[record.id]
            #SeqIO.write(record, open("{}/{}.fasta".format(out, bin), 'a'), "fasta")
            SeqIO.write(record, open("{}/{}.fasta".format(out, bin.replace("; ", "-")), 'a'), "fasta")



def try_to_change_name(sample, name):
    """ The SAM file has different header names than the FASTA file (ugh..) so this is an attempt to change the SAM headers to the FASTA headers without downloading the "map" file and requiring it as an argument.
    """

    #print((sample, name))

    prefix, suffix = sample.split("_")

    #print((prefix, suffix))

    name_digits = name.replace("scaffold", "")

    #print(("digits", name_digits))
    suffix = suffix[:-len(name_digits)]

    #print(("suf", suffix))
    suffix += name_digits

    #print(("new_suf", suffix))
    return "_".join([prefix, suffix])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Writes a set of fasta files, one for each bin, using kraken's classification of contigs.""")

    parser.add_argument("-fasta", help="the fna file of contigs", required=True)
    parser.add_argument("-kraken", help="kraken contigs file", required=True)
    parser.add_argument("-out", help="directory to write binned fastas in", default="fasta_bins/")
    parser.add_argument("-tree", help="a pickled TaxTree object", required=True)
    args = parser.parse_args()



    print("\nGetting sequence lengths...", file=sys.stderr)
    contig_data = get_lengths(args.fasta)

    print("\nLoading TaxTree...", file=sys.stderr)
    # set up the tree to look up taxonomies
    tree = taxtree.TaxTree.load_tree(args.tree)

    print("\nAssigning taxonomy...", file=sys.stderr)
    get_taxonomic_classifications(contig_data, args.kraken, tree)

    print("\nBinning Sequences...\n", file=sys.stderr)
    bins = bin_sequences(contig_data, rank="genus", imerge=True)
   
    print("\nWritting fasta files...", file=sys.stderr)
    write_fastas(bins, args.fasta, args.out)
