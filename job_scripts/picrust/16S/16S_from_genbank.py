
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import os


def find_16S(gbk_f):
    """ Yields all the 16S genes in a GenBank file. """

    with open(gbk_f, 'r') as IN:
        for record in SeqIO.parse(IN, "genbank"):
            for feature in record.features:
                if feature.type == "rRNA":
                    if feature.qualifiers["gene"][0] == "16S":
                        start = feature.location.start.position
                        end = feature.location.end.position
                        
                        if feature.location.strand == -1:
                            yield record.seq[start:end].reverse_complement()
                        else:
                            yield record.seq[start:end]


def main(args):
    with open(args.out, 'w') as OUT:

        for gbk in args.gbk:
            name = os.path.splitext(os.path.basename(gbk))[0]

            best_s = ''
            for s_gene in find_16S(gbk):
                if len(s_gene) > len(best_s):
                    best_s = s_gene

            # make sure a 16S sequence could be found
            if len(best_s) < 800:
                print("No full-length 16S record found in gbk: {}".format(gbk))
            else:
                seq = SeqRecord(seq=best_s, id=name, description="16S gene")
                SeqIO.write(seq, OUT, 'fasta')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates a Fasta file of 16S sequences from GenBank files")
    parser.add_argument("-gbk", help="one or more GenBank files", nargs="+")
    parser.add_argument("-out", help="path to write the Fasta", default="16S_from_GenBank.fasta")

    args = parser.parse_args()

    record = main(args)
