#!/usr/bin/env python3

from Bio import SeqIO
import math
import numpy
from scipy.cluster.vq import kmeans, whiten, vq
import getopt
import sys

matrix = ''

def main(argv):
    fasta = ''
    try:
        opts, args = getopt.getopt(argv, "f:", ["fasta=", "in="])
    except getopt.GetoptError:
        print("USAGE: kmer_cluster.py -f <my.fasta>")
    for opt, arg in opts:
        if opt in ('-f', '--file', '-in'):
            fasta = arg
    print("Making sequences...")
    sequences = make_sequences(fasta)
    #print(sequences[1])

    for seq in sequences:
        print("Beginning profiling for: {}".format(seq.header))
        seq.make_profile_diff()
        #print(seq)
    global matrix

    matrix = make_difference_matrix(sequences, 5)
    #print_table(sequences)
    #clusters = cluster(sequences)
    #for i in range(len(sequences)):
    #    print("{}\t{}".format(sequences[i].header, clusters[i]))

def make_sequences(fasta):
    sequences = []
    for bioSeq in SeqIO.parse(fasta, "fasta"):
        obj = Sequence(bioSeq)
        sequences.append(obj)

    return sequences




class Sequence(object):

    data_lookup = {}        # lookup table that will relate the object sent to the user (key) with the actual object which will be pickled

    def __init__(self, bioSeq):

        self.header = bioSeq.id
        self.seq = bioSeq.seq
        self.nt_freq = {"A": self.seq.upper().count("A"),
                        "C": self.seq.upper().count("C"),
                        "G": self.seq.upper().count("G"),
                        "T": self.seq.upper().count("T")
                        }

        total = sum(self.nt_freq.values())

        for nt in self.nt_freq:
            self.nt_freq[nt] /= total

    def __str__(self):
        atr_list = []
        for key in self.__dict__:
            atr_list.append("{key} = '{value}'".format(key = key, value  = self.__dict__[key]))
        return "\n".join(atr_list)


    def make_profile_diff(self):
        """
        Stores a k-mer profile array for a given sequence.
        """

        self.kmers = {}
        for k in range(5, 6):
            self.kmers[k] = self.count_kmers(k)

    def make_profile(self):
        self.profile = {}
        for k in range(2, 6):
            observed = self.count_kmers(k)
            expected = self.generate_expected(k)
            rel_entropy = 0
            for i in range(len(observed)):
                #pass
                if observed[i] == 0:        #filter out 0's to avoid log(0)
                    pass
                else:
                    rel_entropy += observed[i] * math.log(observed[i] / expected[i])

            #print(observed)
            #print(sum(observed))
            #print(expected)
            #print(sum(expected))
            self.profile[k] = rel_entropy


    def count_kmers(self, k):
        freq = numpy.zeros(4**k, numpy.float)     #make an array of 0s
        total_mer = 0
        for x in range(len(self.seq) + 1 - k):       # +1 for base 0, -k to keep last mer in bounds
            kmer = self.seq[x:x + k]
            if kmer.upper().count("N") > 0:
                continue
            freq[kmer_to_index(kmer)] += 1
            total_mer += 1

        #get frequency per 1000 bases
        (freq / total_mer) * 1000
            # also need to do some sort of correction for different nucleotide composition
            # possibly subtracting the expected values
        return freq

    def generate_expected(self, k):
        """
        Nucleotide composition does not have a direct bearing on information content (you might imagine  I take a sequence and perfectly randomize it such that it has the same nucleotide frequencies but contains no information.

        Therefore, to get information above and beyond random organization of nucleotides, an expected k-mer count is required
        """

        expected = numpy.zeros(4**k, numpy.float)
        #print(len(expected))
        for index in range(4 ** k):
            p = 1
            kmer = index_to_kmer(index, k)
            for base in range(k):
                p *= self.nt_freq[kmer[base]]     #multiple the chance of each base together to get probability of kmer
            expected[index] = p

        return expected




#converts kmers to an alphabettically sorted list index (base 0)
def kmer_to_index(kmer):
    lookup = {"A": 0, "C": 1, "G": 2, "T": 3}
    index = 0
    for i in range(len(kmer)):
        base = kmer[-i - 1]
        index += lookup[base] * 4**i
    return index

def index_to_kmer(index, k):
    lookup = ["A", "C", "G", "T"]
    kmer = ""
    for i in range(k - 1, -1, -1):
        key = index // 4 ** i
        kmer = kmer + lookup[key]

        #print("kmer =", kmer, ", index =", index, ", key =", key)
        if key > 0:
            index %= 4 ** i
    return kmer

#this cluster comparison won't work because clusters aren't constant, it finds different clusters first
def cluster(sequences):
    results = []
    for x in range(100):
         results.append(make_clusters(sequences))

    results = numpy.array(results)
    numpy.mean(results, axis=0)
    return results

def make_clusters(sequences):
    """
    Clusters based on relative entropy (plasmids should have less)
    """

    data = []
    for seq in sequences:
        temp = []
        for k in sequences[0].profile.keys():
            temp.append(seq.profile[k])
        data.append(temp)

    #for k in sequences[0].profile.keys():
    #    temp = []
    #    for seq in sequences:
    #        temp.append(seq.profile[k])
    #    data.append(temp)
    #print(data)
    features = numpy.array(data)
    #print(features)
    whitened = whiten(features)     # scaling of the data (division by stddev) 
    #print(whitened)
    #book = numpy.array((whitened[0],whitened[2]))     # an array that holds the centroids of the cluster (number of clusters = num rows
    #print(book)
    centroids,_ = kmeans(whitened,2)       #call w/ data and clusters
    #print(centroids)

    clusters,_ = vq(whitened, centroids)
    #print(numpy.array(clusters))
    return clusters

def dbscan(sequences):
    pass

def print_table(sequences):
    """
    Print the data to look for trends"
    """

    OUT = open("kmer_cluster_out.txt", 'w')

    #write the headers
    for seq in sequences:
        OUT.write("{}\t".format(seq.header))


    for k in range(2, 6):
        OUT.write("\n")
        OUT.write("{}\t".format(k))
        for seq in sequences:
            OUT.write("{}\t".format(seq.profile[k]))

def make_consensus_matrix(clusters):
    """
    Makes a 3-D matrix that compares the clustering of relative entropy at different values of k
    """
    matrix = []
    for k in clusters:
        columns = []
        for x in clusters[k]:
            row = []
            for y in clusters[k]:
                if x == y:
                    row.append(1)
                else:
                    row.append(0)
            columns.append(row)
        matrix.append(columns)

def calc_genomic_sig_diff(seq1, seq2, k):
    pass



def make_difference_matrix(sequences, k):
    """
    Generates a 2d matrix of the genomic signature difference (at some k) for each sequence against every other sequence.
    """
    matrix = []
    header = [""]
    for seqx in sequences:
        row = []
        header.append(seqx.header)
        for seqy in sequences:
            if row == []:
                row.append(seqx.header)
            row.append(pairwise_difference(seqx, seqy, k))
        matrix.append(row)

    matrix.insert(0, header)

    OUT = open("diff_matrix.txt", "w")
    for row in matrix:
        for element in row:
            OUT.write(str(element) + "\t")
        OUT.write("\n")
    return matrix

def pairwise_difference(sequence1, sequence2, k):
    """
    Takes two sequences and returns the difference in genomic signature at some length k.
    """
    diff = sequence1.kmers[k] - sequence2.kmers[k]

    return abs((1 / 4 ** k) * numpy.absolute(diff).sum())


if __name__ == "__main__":
    main(sys.argv[1:])
