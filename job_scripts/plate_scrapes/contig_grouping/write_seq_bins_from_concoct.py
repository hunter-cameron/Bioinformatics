

import argparse
import sys
import os

from mypyli import utilities

def parse_concoct(concoct_f):
    """ Returns a hash of bins that it parses from a file """
    bins = {}
    with open(concoct_f, 'r') as IN:
        for line in IN:
            header, bin = line[:-1].split(",")
            bins[bin] = bins.get(bin, []) + [header]

    return bins


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", help="clustering file", required=True)
    parser.add_argument("-s", help="clustering software", default="CONCOCT")
    parser.add_argument("-f", help="original fasta", required=True)
    parser.add_argument("-o", help="out directory for fasta files", default=os.getcwd())

    args = parser.parse_args()

    if args.s == "CONCOCT":
        bins = parse_concoct(args.c)

    else:
        raise ValueError("Unknown argument for -s: {}".format(args.s))

   

    # this is the best way to make a directory b/c it handles the possibility that the dir does not
    # exist when it is checked but does by the time the script tries to create it
    try:
        os.makedirs(args.o)
    except OSError as e:
        if e.errno == errno.EEXIST and os.path.isdir(args.o):
            pass
        else: raise


    # note: this is a slow implementation because the file will be read in its entirety for each bin
    # I just don't think speed is necessary.
    for bin in sorted(bins):
        utilities.write_fasta_by_header(fasta=args.f, headers=bins[bin], out="{}/bin{}.fasta".format(args.o, bin))
