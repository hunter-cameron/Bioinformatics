

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

class AnnotatedGene(object):
    """ Simple class to be used as a custom data type """
    
    def __init__(self, name, annot):
        self.name = name
        self.annot = annot
        
        self.scaffold = None
        self.start = None
        self.end = None
        self.strand = None

        self.found_in_gff = False
        self.found_in_fna = False


def get_gene_to_annot(anot_f):
    """ Returns a dict linking gene id to an AnnotatedGene object """
    
    gene2annot = {}
    with open(anot_f, 'r') as IN:
        for line in IN:
            fields = line.rstrip().split("\t")

            gene_id = fields[0]
            annot = fields[2]

            gene2annot[gene_id] = AnnotatedGene(gene_id, annot)

    return gene2annot

def get_scaffold_to_genes(gff_f, genes):
    """ Returns a dict of AnnotatedGenes indexed by scaffold """

    scaf2genes = {}

    with open(gff_f, 'r') as IN:
        for line in IN:
            scaffold, source, type, start, end, score, strand, phase, attributes = line.rstrip().split("\t")

            # first filter by CDS because that is all we are interested in
            if type == "CDS":

                # get the locus tag by converting attributes to a dict
                # each attribute is sep by ";" and then it is a key value pair sep by "="
                # there is also a trailing ";" to worry about
                attributes = attributes.split(";")[:-1]
                attribs = {name: value for name, value in [attrib.split("=") for attrib in attributes]}

                # we will just let this fail if locus_tag isn't there
                locus_tag = attribs["locus_tag"]

                # check if the locus is in the target
                if locus_tag in genes:
                    annot_gene = genes[locus_tag]
                    annot_gene.found_in_gff = True
                    annot_gene.scaffold = scaffold
                    annot_gene.start = int(start)
                    annot_gene.end = int(end)
                    annot_gene.strand = strand

                    # add the locus to the scaffold, might want to make it a simple class rather than a tuple
                    try:
                        scaf2genes[scaffold].append(annot_gene)
                    except KeyError:
                        scaf2genes[scaffold] = [annot_gene]

    return scaf2genes

def write_annotation_fasta(fna_f, scaf2genes, out_f):
    """ Writes a fasta of AnnotatedGenes """

    with open(fna_f, 'r') as IN, open(out_f, 'w') as OUT:
        for record in SeqIO.parse(IN, "fasta"):
            try:
                gene_list = scaf2genes[record.id]
            except KeyError:
                continue

            for gene in gene_list:
                # I'm assuming coordinates are 1 based, inclusive
                seq = record.seq[gene.start-1:gene.end]
                seq_rec = SeqRecord(seq, id="{};annot={}".format(gene.name, gene.annot), description="")
                SeqIO.write(seq_rec, OUT, "fasta")

                gene.found_in_fna = True

def check_completeness(genes):
    gff_missing = []
    fna_missing = []
    for gene in genes:
        if not gene.found_in_gff:
            gff_missing.append(gene.name)
            continue

        if not gene.found_in_fna:
            fna_missing.append(gene.name)

    if gff_missing:
        print("The following annotated genes were missing from the gff file:")
        for name in gff_missing:
            print("   " + name)
    else:
        print("All annotated genes found in gff.")


    if fna_missing:
        print("The following annotated genes were missing from the fna file:")
        for name in fna_missing:
            print("   " + name)
    else:
        print("All annotated genes found in fna.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Writes sequences from a JGI annotation file using the .gff file to link gene names with scaffold position.")
    parser.add_argument("-annot", help="an annotation file from JGI", required=True)
    parser.add_argument("-gff", help=".gff file from JGI", required=True)
    parser.add_argument("-fna", help="FASTA file of scaffolds", required=True)
    parser.add_argument("-out", help="output path for KO fasta", required=True)

    args = parser.parse_args()

    gene2annot = get_gene_to_annot(args.annot)
    scaf2genes = get_scaffold_to_genes(args.gff, gene2annot)
    write_annotation_fasta(args.fna, scaf2genes, args.out)

    # completeness check
    check_completeness(gene2annot.values())
