
import argparse
import sys
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

class MarkerGene(object):

    def __init__(self, aln_marker_f):
        self.fasta = aln_marker_f

    def _read_seqs(self):
        """ Returns a dict of seqs from the file header: seq """
        # read all the seqs in
        seqs = {}
        with open(self.fasta, 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                genome = record.id.split("&&", 1)[0]
                seqs[genome] = str(record.seq)
        return seqs


    def make_distance_matrix(self):
        """ Returns a list of names and a distance matrix in lower triangular format """
        
        # get all the seqs as a dict 'header: seq'
        seqs = self._read_seqs()

        genomes = sorted(seqs)

        distance_matrix = []
        for g1 in genomes:
            row = []
            for g2 in genomes:
                # look for the cue to end the row
                if g2 == g1:
                    row.append(0)
                    break
                else:
                    try:
                        row.append(self._calc_distance(seqs[g1], seqs[g2]))
                    except ZeroDivisionError:
                        print(seqs[g1])
                        print()
                        print(seqs[g2])
                        sys.exit()

            distance_matrix.append(row)


        return genomes, distance_matrix


    def _calc_distance(self, seq1, seq2, max_gap=3):
        """ 
        Counts the number of mismatches between two aligned sequences 
        
        Ignores gaps over the max_gap param.

        Gaps over max_gap will not count towards the length of the alignment or as mismatches,
        gaps under the max_gap will count towards length and mismatches if both aren't gapped.

        What should I do if there are no conserved regions??
        """
        mismatches = 0
        length = 0

        # binary string over a gap -- 1 for both gapped (match) 0 for only one gapped
        gap = ""
        for char1, char2 in zip(seq1, seq2):

            # check if either is a gap
            if char1 == "-" or char2 == "-":
                if char1 == "-" and char2 == "-":
                    gap += "1"
                else:
                    gap += "0"

            else:
                # add the gap if it was short enough, otherwise omit it
                if gap:
                    if len(gap) <= max_gap:
                        length += len(gap)
                        mismatches += gap.count("0")
                    gap = ""

                # check if the current chars match
                if char1 != char2:
                    mismatches += 1

                length += 1

        # add a terminal gap if it was short enough, otherwise omit it
        if gap:
            if len(gap) <= max_gap:
                length += len(gap)
                mismatches += gap.count("0")

        # need to conv to float to force python2 to report the result as float (not int)
        return float(mismatches) / length


def write_newick_tree(names, dist_matrix, outfile="newick.tre"):
    """ Makes a tree from the distance matrix """

    dm = _DistanceMatrix(names, dist_matrix)
    tree_constructor = DistanceTreeConstructor()

    tree = tree_constructor.nj(dm)
    
    Phylo.write(tree, outfile, "newick") 

def test(args):
    mg = MarkerGene(args.markers[0])
    names, matrix = mg.make_distance_matrix()

    #print(names)
    #print(matrix)

    write_newick_tree(names, matrix) 

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Builds a distance matrix using the aligned genes from CheckM ---- This program is untested and remains as a concept that could be further developed.")
    parser.add_argument("-markers", help="fasta files that correspond to aligned markers. CheckM's *.masked.faa files", nargs="+")

    args = parser.parse_args()


    test(args)
