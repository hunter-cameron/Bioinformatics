
import argparse
import sys
import os
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class Gene(object):

    def __init__(self, name, fasta, scaffold, start, end, strand):
        self.name = name
        self.fasta = fasta
        self.scaffold = scaffold
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

        self.seq = None

    @property
    def length(self):
        return self.end - self.start + 1

    def extract_seq(self):
        if self.seq:
            return

        with open(self.fasta, 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                if record.id == self.scaffold:
                    # -1 in start for base 0 to base 1; nothing on end becuase -1 for base 0 and +1 for closed interval
                    self.seq = record.seq[self.start-1:self.end]
                    return
            else:
                raise ValueError("Header {} was not found in fasta".format(self.scaffold))

    def write(self, fh):
        if not self.seq:
            self.extract_seq()

        record = SeqRecord(self.seq, self.name, description="")
        SeqIO.write(record, fh, "fasta")


def get_16S_from_barrnap(fasta_f, barrnap_path):
    output = subprocess.check_output(["perl", barrnap_path, "--quiet", fasta_f], universal_newlines=True)

    for line in output.strip().split("\n"):
        print(line)
        if line.startswith("#"):
            continue

        scaffold, bn_version, type, start, end, evalue, strand, _, annotation = line.split("\t")

        # build a dict from the annotation
        annots = [annot.split("=") for annot in annotation.split(";")]
        annot_dict = {k: v for k, v in annots}

        if annot_dict["Name"] == "16S_rRNA":
            base_name = os.path.splitext(os.path.basename(fasta_f))[0]
            try:
                # I'm not sure the whole scope of notes that can be generated but what I have
                # seen has the structure "aligned only 64 percent of the 16S ribosomal RNA"
                # in my ignornce, I'll split by word and take the 3rd
                note = annot_dict["note"]
                aln_perc = note.split(" ")[2]
                name = base_name + " 16S partial({}%)".format(aln_perc)
            except KeyError:
                name = base_name + " 16S"

            gene = Gene(name, fasta_f, scaffold, start, end, strand)
            yield gene

def main(args):

    to_write = []
    for fasta_f in args.genome:
        genes = [gene for gene in get_16S_from_barrnap(fasta_f, args.barrnap)]

        # report message if no seqs found
        if len(genes) < 1:
            print("No 16S annotated for {}".format(fasta_f))
            continue
        
        # if multiple, right now take longest; eventually may want to do some merging
        gene = genes[0]
        for g in genes:
            if g.length > gene.length:
                gene = g

        to_write.append(gene)

    with open(args.out, 'w') as OUT:
        for gene in to_write:
            gene.write(OUT)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="uses barrnap to annotate 16S genes and then extracts them")
    parser.add_argument("-genome", help="one or more genomic fasta files", nargs="+")
    parser.add_argument("-barrnap", help="barrnap exe path", required=True)
    parser.add_argument("-out", help="output FASTA path", default="extracted_16S.fasta")

    args = parser.parse_args()
    main(args)
