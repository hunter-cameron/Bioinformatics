

import argparse
from Bio import SeqIO
from mypyli import kraken, taxtree
import os

def parse_fasta(fasta_f, contig_data):
    """ Parses fasta for contigs and length and also the name for the metadata"""

    basen = os.path.basename(fasta_f)
    [soil, ecotype, media] = basen.split("_")[:3]

    with open(fasta_f, 'rU') as IN:
        for record in SeqIO.parse(IN, "fasta"):
            contig_data[record.description] = {'length': len(record.seq), 'soil': soil, 'ecotype': ecotype, 'media': media}

def parse_kraken(kraken_f, tree, contig_data):

    for record in kraken.KrakenRecord.parse_kraken_file(kraken_f):
        #print((record.name, record.classified))
        if record.classified:
            contig_data[record.name]['taxonomy'] = tree.lookup_taxid(record.taxid)
        else:
            contig_data[record.name]['taxonomy'] = "unassigned"

def parse_coverage(coverage_f, contig_data):
    with open(coverage_f, 'r') as IN:
        for line in IN:
            if line.startswith("contig"):
                continue

            [contig, mean_cov, med_cov, num_reads] = line[:-1].split("\t")
            
            contig_data[contig].update({'mean_cov': mean_cov, 'med_cov': med_cov, 'num_reads': num_reads})
        
    # make sure all values are initialized
    for contig in contig_data:
        for field in ["mean_cov", "med_cov", "num_reads"]:
            try:
                if not contig_data[contig][field]:
                    contig_data[contig][field] = "0"
            except KeyError:
                contig_data[contig][field] = "0"
        


def write_table(contig_data, out_f):
    with open(out_f, 'w') as OUT:
        OUT.write("\t".join(["contig", "soil", "ecotype", "media", "kingdom", "phylum", "class", "order", "family", "genus", "mean_cov", "median_cov", "length", "num_mapped_reads"]) + "\n")
        for contig in sorted(contig_data):
            tax = []
            for rank in ["kingdom", "phylum", "class", "order", "family", "genus"]:
                taxonomy = contig_data[contig].get('taxonomy', "unassigned")
                if taxonomy == "unassigned":
                    tax.append("unassigned")
                else:
                    tax.append(taxonomy.get_tax_at_rank(rank, "unassigned"))
            OUT.write("\t".join([contig,
                                 contig_data[contig]['soil'],
                                 contig_data[contig]['ecotype'],
                                 contig_data[contig]['media']] +
                                 tax +
                                [str(contig_data[contig]['mean_cov']),
                                 str(contig_data[contig]['med_cov']),
                                 str(contig_data[contig]['length']),
                                 str(contig_data[contig]['num_reads'])]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta", help="fasta file", required=True)
    parser.add_argument("-kraken", help="kraken assignment of the contigs", required=True)
    parser.add_argument("-coverage", help="final coverage file", required=True)
    parser.add_argument("-tree", help="taxtree object to use with kraken", required=True)
    parser.add_argument("-out", help="outfile path for the metadatable", default="metadata_table.txt")

    args = parser.parse_args()


    contig_data = {}

    parse_fasta(args.fasta, contig_data)

    tree = taxtree.TaxTree.load_tree(args.tree)
    parse_kraken(args.kraken, tree, contig_data)

    parse_coverage(args.coverage, contig_data)

    write_table(contig_data, args.out)
    
