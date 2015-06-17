#!/usr/bin/env python3

"""
This is a module used to store a variety of utilities to avoid rewriting the same functions over and over again.
"""
import os
import sys
from Bio import SeqIO

def write_fasta_by_header(fasta, headers=[], new_headers=None, out="fasta_from_headers.fasta"):
    if new_headers:
        if len(headers) != len(new_headers):
            raise AssertionError("Number of new headers doesn't match the number of existing headers.")
        
        header_hash = {k: v for k, v in zip(headers, new_headers)}

    else:

        header_hash = {k: k for k in headers}

    total_headers = len(header_hash)
    with open(fasta, 'rU') as IN, open(out, 'w') as OUT:
        for record in SeqIO.parse(IN, "fasta"):
            if record.id in header_hash:
                record.id = header_hash.pop(record.id)
                #print(record.id)
                record.description = ""
                SeqIO.write(record, OUT, "fasta")
    
    print("New fasta {} written!\n{} of {} headers were not found.".format(out, len(header_hash.keys()), total_headers), file=sys.stderr)

def gbk2faa(gbk, out="genes_from_gbk.faa"):
 
    with open(gbk, 'r') as GBK, open(out, 'w') as OUT:
        for record in SeqIO.parse(GBK, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    # why is this done?
                    assert len(feature.qualifiers['translation']) == 1
                    
                    OUT.write(">{} | {} | start={};end={};product={}\n".format(
                        feature.qualifiers['locus_tag'][0], 
                        record.name,
                        feature.location.start,
                        feature.location.end,
                        feature.qualifiers['product'][0]
                        ))
                    OUT.write(feature.qualifiers['translation'][0] + "\n")
    
