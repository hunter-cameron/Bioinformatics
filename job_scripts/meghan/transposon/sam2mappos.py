

import sys
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sam = sys.argv[1]

elem_dict = {}
with open(sam, 'r') as IN:
    for line in IN:
        if line.startswith("@"):
            continue
        elements = line[:-1].split("\t")

        # skip unmapped
        if elements[1] == '4':
            continue

        name = elements[0]
        pos = elements[3]
        cigar = elements[5]
        qual = elements[4]
        elem_dict[name] = {'pos': pos, 'qual': qual, 'cigar': cigar}


def get_genes_from_gbk(gbk=sys.argv[2]):
    genes = []
    with open(gbk, 'r') as IN:
        for record in SeqIO.parse(IN, "genbank"):
            for feature in record.features:
                if feature.type in ["CDS", "Gene"]:
                    start = feature.location.start
                    end = feature.location.end
                    locus_tag = feature.qualifiers["locus_tag"][0]

                    genes.append({'start': start, 'end': end, 'locus_tag': locus_tag})

    return genes


def print_pos(elem_dict, genes, qual=20):
    gene_counts = {}
    print("\t".join(["name", "position", "cigar", "quality", "gene"]))
    for seq in sorted(elem_dict, key=lambda k: int(elem_dict[k]['pos'])):
        if int(elem_dict[seq]['qual']) < qual:
            continue
        pos = int(elem_dict[seq]['pos'])
        for gene in genes:
            #print((int(gene['start']), int(gene['end'])))
            if int(gene['start']) < pos <= int(gene['end']):
                locus = gene['locus_tag']
                gene_counts[locus] = gene_counts.get(locus, 0) + 1
                break
        else:
            locus = "None"

        print("\t".join([seq, elem_dict[seq]['pos'], elem_dict[seq]['cigar'], elem_dict[seq]['qual'], locus]))


    uniq_genes = len(gene_counts)

    # make sure all genes are included
    for gene in genes:
        gene_counts[gene['locus_tag']] = gene_counts.get(gene['locus_tag'], 0)

    # plot gene counts
    to_plt = [gene_counts[gene['locus_tag']] for gene in sorted(genes, key=lambda k: int(k['start']))]
    ax = plt.figure().add_subplot(111)
    plt.xlabel("gene")
    plt.ylabel("count")
    plt.bar(range(len(to_plt)), height=to_plt)
    plt.xlim([0, len(to_plt)])
    plt.text(.8, .8, "N={}\nUniq={}\nCov={:.2f}X".format(sum(to_plt), uniq_genes, sum(to_plt) / len(to_plt)), horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    #plt.legend(["N={}".format(sum(to_plt)),"Avg={}".format(sum(to_plt) / len(to_plt))], loc=1)
    plt.savefig("counts_per_gene.png", bbox_inches='tight', dpi=200)

def plot_by_position(elem_dict, qual=20):
    """ Scatterplot of count of alignments starting at that position and position. """
    
    pos_counts = {}
    for seq in sorted(elem_dict, key=lambda k: int(elem_dict[k]['pos'])):
        if int(elem_dict[seq]['qual']) < qual:
            continue
        pos = elem_dict[seq]['pos']
        pos_counts[pos] = pos_counts.get(pos, 0) + 1
        
    
    x_vals = [int(pos) for pos in sorted(pos_counts, key=lambda k: int(k))]
    y_vals = [pos_counts[pos] for pos in sorted(pos_counts, key=lambda k: int(k))]

    ax = plt.figure().add_subplot(111)
    plt.bar(x_vals, y_vals)
    #ax.set_yscale('log', basey=2)
    plt.xlim([0, 6169072])
    plt.ylim([0, max(y_vals) + 5])
    plt.ylabel("count")
    plt.xlabel("position on genome")
    plt.text(.8, .8, "N={}\nUniq={}".format(sum(y_vals), len(x_vals)), horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
    #plt.show()
    plt.savefig("plot_by_position.png", bbox_inches='tight', dpi=200)
    
def write_gene_hist(elem_dict, genes, qual=20):
    gene_counts = {}
    for seq in sorted(elem_dict, key=lambda k: int(elem_dict[k]['pos'])):
        if int(elem_dict[seq]['qual']) < qual:
            continue
        pos = int(elem_dict[seq]['pos'])
        for gene in genes:
            #print((int(gene['start']), int(gene['end'])))
            if int(gene['start']) < pos <= int(gene['end']):
                locus = gene['locus_tag']
                gene_counts[locus] = gene_counts.get(locus, 0) + 1
                break
        else:
            locus = "None"
  

    make_hist([gene_counts[gene] for gene in gene_counts])

def make_hist(data):
    """ Function I had to make because there isn't a single useful histogram program that Google knows about......"""

    bins = [1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 100, 10000]

    hist, h_bins = np.histogram(data, bins=bins)
    
    # begin plotting 
    ax = plt.figure().add_subplot(111)
    print(hist)
    print(h_bins)
    plt.bar(range(len(hist)), height=hist, align='center')
    plt.xlim([-1, len(hist)])
    plt.xticks(range(len(hist)), [str(bin) for bin in bins])
    plt.xlabel("bin")
    plt.ylabel("count")
    #plt.show()
    plt.savefig("hist_of_gene_counts.png", bbox_inches='tight', dpi=200)

def make_qual_hist(elem_dict):
    bins = [1, 5, 11, 21, 31, 41, 51]

    data = [int(record['qual']) for record in elem_dict.values()]
    
    hist, h_bins = np.histogram(data, bins=bins)

    # begin plotting 
    ax = plt.figure().add_subplot(111)
    print(hist)
    print(h_bins)
    plt.bar(range(len(hist)), height=hist, align='center')
    plt.xlim([-1, len(hist)])
    plt.xticks(range(len(hist)), [str(bin) for bin in bins])
    plt.xlabel("quality")
    plt.ylabel("count")
    #plt.show()
    plt.savefig("hist_of_mapping_quality.png", bbox_inches='tight', dpi=200)

def count_hist(data, bins):
    counts = [0] * len(bins)
    for value in data:
        for indx, bin in enumerate(bins):
            if value < bin:
                counts[indx] += 1
                break

    return counts


qual = 20

make_qual_hist(elem_dict)
genes = get_genes_from_gbk(sys.argv[2])
print_pos(elem_dict, genes, qual)
plot_by_position(elem_dict, qual)
#write_gene_hist(elem_dict, genes)
