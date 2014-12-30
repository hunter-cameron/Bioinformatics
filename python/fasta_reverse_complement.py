#!/usr/bin/env python3

from Bio import SeqIO
import sys

args = sys.argv[1:]

handle = open(args[0], 'rU')
out = open(args[1], 'w')
for record in SeqIO.parse(handle, 'fasta'):
    record.id = record.id + "_rc"
    record.seq = record.seq.reverse_complement()
    record.description = ''     # if the description is different than the header, it prints it?
    SeqIO.write(record, out, 'fasta')
