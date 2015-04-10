
from mypyli import blastparser
from Bio import SeqIO
import argparse
import logging
import sys

def get_seq_ids(fasta_f):
    """ Returns a list of sequence headers and makes sure they are all unique."""
    ids = []
    with open(fasta_f, 'r') as IN:
        for seq in SeqIO.parse(IN, 'fasta'):
            ids.append(seq.id)

    if len(ids) != len(set(ids)):
        raise AssertionError("One or more sequence headers are duplicates!")
    else:
        return ids


def get_top_hits_from_blast(blast_f, max_e, min_id):
    top_hits = {}
    logging.info(blast_f)
    with open(blast_f, 'r') as IN:
        for record in blastparser.parse(IN, outfmt='6'):
            if record.perc_id >= min_id and record.evalue <= max_e:
                prev_hit = top_hits.get(record.query, (0, 0))

                if prev_hit[1] < record.evalue:
                    top_hits[record.query] = (record.subject, record.evalue)

    return top_hits

def print_taxonomy_assignments(out_f, seq_ids, hits, taxonomy_assignments):
    with open(out_f, 'w') as OUT:
        for seq in seq_ids:
            if seq in hits:
                taxonomy = taxonomy_assignments[hits[seq][0]]
                evalue = hits[seq][1]
            else:
                taxonomy = "Unassigned"
                evalue = 0

            OUT.write("\t".join([seq, taxonomy, str(evalue)]) + "\n")

def get_taxonomy(seq2tax_f, seq_ids):
    taxonomy = {}
    with open(seq2tax_f, 'r', errors='ignore') as IN:
        for line in IN:
            seq, tax = line[:-1].split("\t")
            #print((seq, tax))
            if seq in seq_ids:
                taxonomy[seq] = tax
    return taxonomy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is an implementation of QIIME's assign_taxonomy.py -m blast pipeline. This is required because the QIIME BLAST wrapper required legacy BLAST and all we have is BLAST+.")
    parser.add_argument("-e", type=float, help="maximum E-value to consider (default: 1E-30)", default=1E-30)
    parser.add_argument("-id", type=float, help="minimum identity to consider (default: 90)", default=90.0)
    parser.add_argument("-t", help="tab-delimited file that maps sequences to taxonomy", required=True)
    parser.add_argument("-r", help="FASTA file of reference sequences", required=True)
    parser.add_argument("-input", "-i", help="input FASTA file", required=True)
    parser.add_argument("-out", help="output classifications file", default="blast_taxonomy_classifications.txt")
    parser.add_argument("-b", help="BLAST output file in outfmt 6")
    args = parser.parse_args()

    # turn id into a percent if necessary
    if args.id < 1:
        args.id *= 100

    # run the BLAST if needed
    if not args.b:
        raise AssertionError("Currently, you must feed the program a BLAST file!")
 
    #get_taxonomy(args.t, [])

    logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)

    logging.info("Getting seq ids from input...")
    seq_ids = get_seq_ids(args.input)


    logging.info("Getting top hit from each id...")
    hits = get_top_hits_from_blast(args.b, args.e, args.id)

    logging.info("Looking up taxonomy for top hits...")
    
    taxonomy_assignments = get_taxonomy(args.t, set([hits[hit][0] for hit in hits]))

    logging.info("Printing results to {}".format(args.out))
    print_taxonomy_assignments(args.out, seq_ids, hits, taxonomy_assignments)
