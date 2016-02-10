
from Bio import AlignIO
import argparse
import sys

parser = argparse.ArgumentParser(description="prints a region from a MSA")
parser.add_argument("-msa", help="msa file", required=True)
parser.add_argument("-start", help="inclusive 0 based start", required=True, type=int)
parser.add_argument("-end", help="inclusive 0 based end", required=True, type=int)

args = parser.parse_args()

with open(args.msa, 'r') as IN:
    alignment = AlignIO.read(IN, "fasta")

    AlignIO.write(alignment[:, args.start:args.end+1], sys.stdout, "fasta")
