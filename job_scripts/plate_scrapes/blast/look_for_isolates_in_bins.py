

from mypyli import blastparser, isolatedb
import argparse
import sys
import os
import logging

from Bio import SeqIO
MIN_ID = 95
MIN_LEN = 1000




def process_blast_file(blast_f):

    isolates = {'total': {'contigs': 0, 'length':0}}
    total = {'contigs': 0, 'length':0}
    with open(blast_f, 'r') as IN:
        for record in blastparser.parse(IN, outfmt="6"):
            if record.perc_id >= MIN_ID and record.length >= MIN_LEN:
                subject = record.subject
                isolate = subject.split("_")[0]
            
                isolates[isolate] = isolates.get(isolate, {'contigs': 0, 'length': 0})
                isolates[isolate]['contigs'] += 1
                isolates[isolate]['length'] += record.length


            isolates['total']['contigs'] += 1
            isolates['total']['length'] += record.length

    return isolates

def main():
    bins = {}
    for blastf in sys.argv[1:]:
        isolates = process_blast_file(blastf)

        bin = os.path.basename(blastf).split(".")[0]

        print(bin)
        print("\tTotal - Contigs: {}, Length: {}".format(str(isolates['total']['contigs']), str(isolates['total']['length'])))
        
        for isolate in isolates:
            if isolate == "total":
                continue
            if isolates[isolate]['length'] > .5 * isolates['total']['length']:

                print("\t{} - Contigs: {}, Length: {}. ***".format(isolatedb.convert_values([isolate])[0], str(isolates[isolate]['contigs']), str(isolates[isolate]['length'])))
            else:
                print("\t{} - Contigs: {}, Length: {}".format(isolatedb.convert_values([isolate])[0], str(isolates[isolate]['contigs']), str(isolates[isolate]['length'])))

        print()

def make_contig_dict(bins):
    contig_dict = {}
    for bin in bins:
        bin_name = os.path.basename(bin).rsplit(".", 1)[0]
        with open(bin, 'r') as IN:
            for line in IN:
                if line.startswith(">"):
                    if line[1:-1] in contig_dict:
                        logging.warning("Contig {} is repeated!".format(line[1:-1]))
                    else:
                        contig_dict[line[1:-1]] = bin_name
                        #print(line)

    return contig_dict

def locate_isolates(contigs, isolates):
    for isolate_f in isolates:
        bin_counts = {}
        with open(isolate_f, 'r') as IN:

            isolate_id = os.path.basename(isolate_f).rsplit(".", 1)[0]
            for seq in SeqIO.parse(isolate_f, "fasta"):
                
                try:
                    bin = contigs[str(seq.id)]
                except KeyError:
                    #print(str(seq.id))
                    bin = "unbinned"
                    
                #print(len(seq.seq))
                bin_counts[bin] = bin_counts.get(bin, 0) + len(seq.seq)

        print("Isolate {}".format(isolatedb.convert_values([isolate_id])[0]))
        for bin in sorted(bin_counts, key=lambda k: bin_counts[k], reverse=True):
            print("\t{}: {}".format(bin, str(bin_counts[bin])))

        print()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-bins", nargs="+", help="all the bin fasta files", required=True)
    parser.add_argument("-isolates", nargs="+", help="all the isolate", required=True)
    args = parser.parse_args()

    contigs = make_contig_dict(args.bins)
    locate_isolates(contigs, args.isolates)

