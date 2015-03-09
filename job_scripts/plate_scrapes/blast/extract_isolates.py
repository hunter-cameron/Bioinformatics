

from mypyli import blastparser, utilities
import sys
import os

MIN_ID = 90
MIN_LEN = 1000

blast_f = sys.argv[1]
fasta_f = sys.argv[2]

def write_isolate_fastas(isolates):

    for isolate in isolates:
        utilities.write_fasta_by_header(fasta=fasta_f, headers=isolates[isolate], out="meta_{}.fasta".format(isolate))

isolates = {}
with open(blast_f, 'r') as IN:
    for record in blastparser.parse(IN, "6"):
        if record.get_perc_id() >= MIN_ID and record.get_length() >= MIN_LEN:
            subject = record.get_subj()
            isolate = subject.split("_")[0]
            
            contig = record.get_query()
            
            isolates[isolate] = isolates.get(isolate, []) + [contig]


write_isolate_fastas(isolates)
            



