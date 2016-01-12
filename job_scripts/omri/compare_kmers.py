
import argparse
import os
import sys
from mypyli import jellyfish
import numpy as np
import logging

import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel("DEBUG")

class KmerTree(object):
    
    def __init__(self, k, errors=0):
        self.k = k
        self.errors = errors
    
        self.tree = None

    @classmethod
    def _recurse_all_kmers(cls, k, curseq="", collected_seqs=[]):
        """ Generates all kmers to avoid having to store a large array """
        if len(curseq) == k:
            yield Kmer(curseq)
        else:
            for nt in ["A", "C", "G", "T"]:
                for kmer in cls._recurse_all_kmers(k, curseq + nt, collected_seqs):
                    yield kmer
               

    def _recurse_add_nodes_w_errors(self, kmer, remaining_seq=None, errors=0, tree={}):
        if remaining_seq is None:
            remaining_seq = kmer.seq
        
        # if error is allowed, add under all possibilities
        if errors < self.errors:
            for nt in ["A", "C", "G", "T"]:
                if len(remaining_seq) == 1:
                    node = tree.add_child(nt)
                    node.add_data(kmer)
                else:   # more to process
                    node = tree.add_child(nt)

                    # see if error is added
                    if nt == remaining_seq[0]:
                        self._recurse_add_nodes_w_errors(kmer, remaining_seq[1:], errors, node) 
                    else:
                        self._recurse_add_nodes_w_errors(kmer, remaining_seq[1:], errors+1, node) 

        else:
            if len(remaining_seq) == 1:
                node = tree.add_child(remaining_seq[0])
                node.add_data(kmer)
            else:
                node = tree.add_child(remaining_seq[0])
                self._recurse_add_nodes_w_errors(kmer, remaining_seq[1:], errors, node)

        return tree

    def build_tree(self):
        
        tree = TreeNode()
        for kmer in self._recurse_all_kmers(self.k):
            self._recurse_add_nodes_w_errors(kmer, tree=tree)

        self.tree = tree
        tree.print_tree()

class TreeNode(object):

    def __init__(self, name="", parent=None):
        self.name = name
        self.parent = parent
        self.children = {}
   
        self.data = []

    def __str__(self):
        return self.name
    
    def print_tree(self, level=0):
        print("  " * level + str(self))
        if self.data:
            for elem in self.data:
                print("  " * level + str(elem))

        for child in self.children.values():
            child.print_tree(level=level+1)


    def add_child(self, key, child=None):
        """ Adds a child to a given key and returns. If key already exists, returns the child """
        try:
            return self.children[key]
        except KeyError:
            if child is None:
                self.children[key] = TreeNode(self.name + key, self)
            else:
                self.children[key] = child

            return self.children[key]

    def add_data(self, data):
        self.data.append(data)

class Kmer(object):
    def __init__(self, seq):

        self.seq = seq

        self.genomes = []

    def __str__(self):
        return "Kmer {}".format(self.seq)


def get_taxonomy(metadata_f):

    ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    genome2taxonomy = {}
    with open(metadata_f, 'r') as IN:
        for line in IN:
            gbk_f, genome_name, taxstr = line.rstrip().split("\t")

            name = os.path.splitext(os.path.basename(gbk_f))[0]

            tax_elems = taxstr.split("; ")

            taxonomy = {rank: elem for rank, elem in zip(ranks[:len(tax_elems)-1], tax_elems)}

            genome2taxonomy[name] = taxonomy


    return genome2taxonomy


def count_conserved_per_tax(kmers, genomes2tax, tax_level):
    # make lists of genomes for each group in the tax level
    genome_groups = {}
    for genome in genomes2tax:
        try:
            group = genomes2tax[genome][tax_level]
        except KeyError:
            LOG.warning("No tax level {} for genome {}".format(tax_level, genome))
            continue

        try:
            genome_groups[group].append(genome)
        except KeyError:
            genome_groups[group] = [genome]

    # reindex the groups by number of genomes (remove groups with n<=1)
    for group in list(genome_groups.keys()):
        n = len(genome_groups[group])

        if n <= 1:
            del genome_groups[group]
        else:
            new_key = group + "(n={})".format(n)
            genome_groups[new_key] = genome_groups[group]
            del genome_groups[group]



    # go through all kmers and get conserved counts
    group_counts = {group: 0 for group in genome_groups}
    # look at each kmer
    for kmer in kmers:
        # look at each taxonomic group
        for group in genome_groups:
            # check if each genome from group was found in the kmer
            for genome in genome_groups[group]:
                if genome not in kmers[kmer]:
                    break
            else:
                group_counts[group] += 1

    return group_counts


def count_conserved_per_tax_disk(genome2kmer, genomes2tax, tax_level, core_frac):
    # make lists of genomes for each group in the tax level
    genome_groups = {}
    for genome in genomes2tax:
        try:
            group = genomes2tax[genome][tax_level]
        except KeyError:
            LOG.warning("No tax level {} for genome {}".format(tax_level, genome))
            continue

        try:
            genome_groups[group].append(genome)
        except KeyError:
            genome_groups[group] = [genome]

    # reindex the groups by number of genomes (remove groups with n<=1)
    for group in list(genome_groups.keys()):
        n = len(genome_groups[group])

        if n <= 1:
            del genome_groups[group]
        else:
            new_key = group + "(n={})".format(n)
            genome_groups[new_key] = genome_groups[group]
            del genome_groups[group]



    # go through all kmers and get conserved counts
    group_counts = {group: 0 for group in genome_groups}
    for group in genome_groups:
        LOG.info("Getting conserved kmers for {}".format(group))
        kmer_dict = {}
        num_genomes = len(genome_groups[group])
        num_required = num_genomes * core_frac
        for index, genome in enumerate(genome_groups[group], 0):
            LOG.debug("Processing genome {}".format(genome))
            # read in each kmer from that genome
            with open(genome2kmer[genome], 'r') as IN:
     
                for line in IN:
                    kmer, count = line.rstrip().split("\t")
                    try:
                        kmer_dict[kmer] += 1
                    except KeyError:
                        # skip adding if there aren't enough genomes left to make it count
                        if num_genomes - index < num_required:
                            continue
                        else:
                            kmer_dict[kmer] = 1

        # find the shared kmers
        shared_kmers = []
        for kmer in kmer_dict:
            if kmer_dict[kmer] >= num_required:
                shared_kmers.append(kmer)

        group_counts[group] = len(shared_kmers)

    return group_counts


def make_bar_graph(counts_dict, tax_rank):
    taxa = sorted(list(counts_dict.keys()))
    counts = [counts_dict[tax] for tax in taxa]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_title("Conserved 20-mers at {} level".format(tax_rank))
    try:
        ax.bar(range(len(taxa)), counts, 1, align="center", log=True)
    except:
        LOG.warning("Making Bar graph failed for tax_rank '{}'. Printing dict of counts to stdout".format(tax_rank))
        print(counts_dict)

    plt.xlim(-1, len(taxa))
    plt.xticks(range(len(taxa)), taxa, rotation="vertical", ha='center', fontsize='small')
    plt.xlabel("Taxa")

    plt.ylabel("# Conserved 20mers")
    #ax.yaxis.set_ticks_position("right")

    fig.savefig("conserved_kmer.{}.png".format(tax_rank))


def main(fastas, k, metadata_f, core_frac):
    # count the kmers
    fasta2kmer_counts = {}
    for fasta in fastas:
        fasta_name = os.path.splitext(os.path.basename(fasta))[0]
        if os.path.isfile(fasta_name + ".counts.txt"):
            LOG.info("Found existing counts file for {}".format(fasta_name))
            fasta2kmer_counts[fasta_name] = fasta_name + ".counts.txt"
            continue

        LOG.info("Counting kmers for {}".format(fasta))

        jf_indx = jellyfish.count_kmers(fasta, k, prefix=fasta_name)
        kmer_counts = jellyfish.dump_kmer_counts(jf_indx, output_f=fasta_name + ".counts.txt")

        fasta2kmer_counts[fasta_name] = kmer_counts

    
    """ REMOVED BECAUSE THIS TAKES TOO MUCH MEMORY
    # read all kmer profiles into memory
    all_kmers = {}
    for fasta_name in fasta2kmer_counts:
        LOG.info("Reading kmers from {} into memory.".format(fasta_name))
        with open(fasta2kmer_counts[fasta_name], 'r') as IN:
            for line in IN:
                kmer, count = line.rstrip().split("\t")
                
                try:
                    all_kmers[kmer][fasta_name] = int(count)
                except KeyError:
                    all_kmers[kmer] = {fasta_name: int(count)}

    """


    # cycle through kmers by taxonomy
    genome2taxonomy = get_taxonomy(metadata_f)
    for taxonomic_level in ["kingdom", "phylum", "class", "order", "family"]:
        LOG.info("Generating counts by taxonomic rank for {}".format(taxonomic_level))

        #counts_per_tax = count_conserved_per_tax(all_kmers, genome2taxonomy, taxonomic_level)
        counts_per_tax = count_conserved_per_tax_disk(fasta2kmer_counts, genome2taxonomy, taxonomic_level, core_frac)
        make_bar_graph(counts_per_tax, taxonomic_level)

def test_bargraph():
    data = {"samp1": 0}
    #make_bar_graph(data, "test1")

    #data["samp1"] = 1
    #make_bar_graph(data, "test2")

    data.update({"samp2": 50, "samp3": 4, "samp4": 4000})
    make_bar_graph(data, "test3")
            
    
""" 
Make Kmer object that is low weight and stores basically just which genomes it is in. 

Then, make trie of all paths each kmer could be allowing x mismatches. Bottom node of trie has a list of indices to an array of kmer objects. 

Perhaps init all possible Kmer objects at once. 
"""

""" For a simple comparison, I could just make a np boolean array of found or not and then compare masks using an AND logical comparison. 
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fastas", help="a bunch of fastas to compare", nargs="+", required=True)
    parser.add_argument("-k", help="value of k to use", required=True, type=int)
    parser.add_argument("-metadata", help="file with genome in the first column and a tax string in the third", required=True)
    parser.add_argument("-core", help="fraction of genomes required for presence", type=float, default=1)
    args = parser.parse_args()

    main(args.fastas, args.k, args.metadata, args.core)
