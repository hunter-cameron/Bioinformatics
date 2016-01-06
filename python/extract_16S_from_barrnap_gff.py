
description = """
Author: Hunter Cameron
Date: December 2015
Description:
Small utility for extracting 16S genes from one or more 
barrnap *.gff files.
"""

import argparse
from Bio import SeqIO
import os

class Gene(object):
    """ Data container for a 16S gene """

    def __init__(self, name, contig, start, end, strand, note):
        self.name = name
        self.contig = contig
        self.start = start
        self.end = end
        self.strand = strand
        self.note = note

        # placeholder to be set later
        self.seq = None


    def generate_header(self):
        """ Puts together the data to generate a header """

        header = ";".join([
                ">" + self.name + " 16S",  
                self.contig + "[{self.start}..{self.end}]{self.strand}".format(self=self), 
                ])

        if self.note:
            header += ";" + self.note

        return header


def parse_gff(gff_f):
    """ Parses the GFF and returns a list of Gene objects for the things annotated as 16S """

    genes = []
    with open(gff_f, 'r') as IN:
        for line in IN:

            # skip comments
            if line.startswith("#"):
                continue
            
            # split out the data and make a dict of whatever keys are present in the annotation
            contig, version, type, start, end, evalue, strand, _, annotation = line.rstrip().split("\t")
            annotation_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in annotation.split(";")}

            # check if this is 16S
            if annotation_dict["Name"] == "16S_rRNA":
                # now that we know this is 16S, make it a name
                name = os.path.basename(os.path.splitext(gff_f)[0])

                # check for a note; a note indicates a partial sequence
                try:
                    note = annotation_dict["note"]
                except KeyError:
                    note = ""

                # build the Gene object and add it to the list
                gene = Gene(name, contig, int(start), int(end), strand, note)
                genes.append(gene)
            
    return genes           

def get_genes_from_fasta(fasta_f, genes):
    with open(fasta_f, 'r') as IN:
        contigs = []

        # cycle through contigs in the file and write all the genes from each
        for record in SeqIO.parse(IN, "fasta"):
            contigs.append(record.id)
            genes_in_this_contig = [gene for gene in genes if gene.contig == record.id]

            if genes_in_this_contig:
                # set the seq attribute of each gene
                for gene in genes_in_this_contig:
                    # -1 because positions are 1 based (arrays 0 based)
                    seq = record.seq[gene.start-1:gene.end]

                    if gene.strand == "-":
                        seq = seq.reverse_complement()

                    gene.seq = seq


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-gffs", help="gff files, name should be name.gff", nargs="+", required=True)
    parser.add_argument("-fastas", help="fasta files, name should be name.fasta", nargs="+", required=True)
    parser.add_argument("-out", help="path for a fasta of 16S sequences")

    args = parser.parse_args()

    # get the positions of the genes per file
    genes_per_genome = {}
    for gff in args.gffs:
        name = os.path.basename(os.path.splitext(gff)[0])
        genes = parse_gff(gff)
        
        if not genes:
            print("Found no 16S genes for {}".format(gff))

        genes_per_genome[name] = genes


    # get the seq for each gene
    for fasta in args.fastas:
        name = os.path.basename(os.path.splitext(fasta)[0])

        if name not in genes_per_genome:
            print("{} does not have a corresponding GFF. Skipping.".format(fasta))
            continue

        if genes_per_genome[name]:
            get_genes_from_fasta(fasta, genes_per_genome[name])


    # write the output fasta
    missing_genes = []
    with open(args.out, 'w') as OUT:
        for genome in genes_per_genome:
            for gene in genes_per_genome[genome]:
                if gene.seq is None:
                    missing_genes.append(gene)
                else:
                    OUT.write(gene.generate_header() + "\n" + str(gene.seq) + "\n")

    # check for genes that were not found in the fasta
    if missing_genes:
        print("Didn't find seq for the following genes.")
        for gene in missing_genes:
            print("  " + gene.generate_header())

    else:
        print("Found seqeunces for all genes.")
