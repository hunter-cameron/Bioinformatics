
import argparse

from Bio import AlignIO
from Bio import SeqUtils
from Bio.Align import MultipleSeqAlignment 
import os
import sys
import math

class Primer(object):
    def __init__(self, seq):
        self.seq = seq
        
        self.melting_temp = self.calc_melting_temp()


    def calc_melting_temp(self):
        """ Calculates melting temp according to X's paper and unified values """
        SeqUtils.MeltingTemp.Tm_NN(self.seq, dnac2=0)


class Consensus(object):
    def __init__(self, msa):
        self.msa = msa

        # placeholders set in find_consensus()
        self.seq = None
        self.entropy = None

        self.find_consensus()

    def __str__(self):
        return self.seq

    def find_consensus(self, max_entropy=0):
        consensus = ""
        entropy = []
        
        # process each column of the MSA
        for i in range(len(self.msa[0])):
            col = self.msa[:, i]
          
            # get the most frequent character and the entropy (I should probably decouple this...)
            consensus_nt, ent = self.shannon_entropy(col)
            
            # set the consensus character if entropy is low enough
            if ent <= max_entropy:
                consensus_symbol = consensus_nt
            else:
                consensus_symbol = "-"

            consensus += consensus_symbol
            entropy.append(ent)

        self.seq = consensus
        self.entropy = entropy

    @staticmethod
    def shannon_entropy(seq):
        """ Returns a tuple with the most common character and the shannon entropy for a given seq """

        # fields for selecting the character to use for the consensus
        highest_prob = 0
        nt_name = ""
            
        # entropy running sum
        entropy = 0

        # store seq_len to avoid repeated calculations
        seq_len = len(seq)
        for nt in set(seq):
            prob = seq.count(nt) / seq_len

            # see if this nt has the new highest probability
            if prob > highest_prob:
                highest_prob = prob
                nt_name = nt
            
            # add the entropy weight to the running sum
            entropy += prob * math.log(prob, 2)

        return nt_name, entropy * -1


class PotentialAmplicon(object):
    def __init__(self, start, end, score, parent):
        self.start = start
        self.end = end
        self.score = score
   
        self.bounding_left = None
        self.bounding_right = None

        self.parent = parent
        
    def __str__(self):
        return "PotentialAmplicon {self.start}..{self.end}  score: {self.score}".format(self=self)

    def get_uniqueness(self):
        """ Makes a distance matrix and returns the lowest value """
        
        # arbitrarily high starting number
        min_distance = 99999

        # get only the portion of the seq corresponding to this amplicon
        # add 1 to the end because half openness
        seqs = self.parent.get_region(self.start, self.end)

        # calculate the distance using lower triangular algorithm
        #print()
        #print()
        #print("Here1")
        for seq1 in seqs:
            #print()
            
            for seq2 in seqs:

                # lower triangular algorithm
                if seq1 == seq2:
                    break


                # calculate the distance between the two seqs and set min_dist if necessary
                distance = 0
                for nt1, nt2 in zip(seq1, seq2):
                    if nt1 != nt2:
                        distance += 1

                #print(str(distance),end="\t")

                if distance < min_distance:
                    min_distance = distance
                
                # if min = 0 we know that the uniqueness is 0 so no need to process any more
                if min_distance == 0:
                    return 0

        return min_distance


    def print_mask(self, fh=sys.stdout):
        """ Prints a "mask" of the amplicon where the bounding left and right regions are the only NT printed and the variable region is represented by Ns """

        # + 1 to get len from distance
        num_Ns = (self.end - self.start) + 1

        fwd_bound = self.parent.get_region(self.bounding_left.start, self.bounding_left.end)[0].seq
        rev_bound = self.parent.get_region(self.bounding_right.start, self.bounding_right.end)[0].seq

        mask = str(fwd_bound) + "N"*num_Ns + str(rev_bound)

        header = ">" + "PotentialAmplicon {} {}..{}  score: {}".format(self.parent.name, self.bounding_left.start, self.bounding_right.end, self.score) + "  uniqueness: {}".format(self.get_uniqueness())

        fh.write(header + "\n" + mask + "\n")
        



class ConservedRegion(object):
    """ 
    This handles a ConservedRegion as represented by a MSA.
    
    The basic workflow to find primers is:

    1. Find a consensus sequence for the MSA
    2. Identify ultra_conserved (currently no mismatches) regions
    3. Check inter_ultra regions to see if there are any that are of the right length
    4. 
    """
    
    def __init__(self, msa_f, name=None):
        # read and filter duplicates from the MSA_F
        self.msa = self.remove_duplicates(self._read_msa(msa_f))

        if name:
            self.name = name
        else:
            self.name = os.path.splitext(os.path.basename(msa_f))[0]

        self.consensus = Consensus(self.msa)


    @staticmethod
    def _read_msa(msa_f, aln_format="fasta"):
        return AlignIO.read(msa_f, aln_format)

    @ staticmethod
    def remove_duplicates(msa):
        """ There are some duplicate genomes in Omri's list. This is a custom function to remove them. 
        I want to remove: 378, CL136, and CL126
        
        """
        records = []
        for record in msa:
            header = record.description
            org_name = header.split("__")[1]

            if org_name in ["Methylobacterium_CL126", "Methylobacterium_CL136", "Methylobacterium_378"]:
                continue
            else:
                records.append(record)

        return MultipleSeqAlignment(records)


    def get_region(self, start, end):
        """ Returns a region of the msa """

        # +1 to convert full open to half open range
        return self.msa[:, start:end+1]

    def get_uniqueness(self):
        """ Gets the uniqueness over the whole region """
        amp = PotentialAmplicon(0, len(self.msa[0]), 0, self)
        uniq = amp.get_uniqueness()
        return uniq

    def find_ultra_conserved(self, min_len):
        """ Finds ultra-conserved regions of at least min_len targeted at being ideal primer sites"""
        
        class UltraConserved(object):
            """ Simple class to hold info about UC regions """
            def __init__(self, start, end):
                self.start = start
                self.end = end

            def __str__(self):
                return "UltraConserved Region {self.start}..{self.end}".format(self=self)
        
            def distance(self, other):
                """ Returns the distance between the end of one and the beginning of another as a positive value or 0 if the two overlap """

                # three possible options self is before, after, or overlaps other
                
                # first determine if they overlap
                if self.start <= other.end and other.start <=self.end:
                    return 0

                # now that we know there is no messy over lap, we need to see which is first
                elif self.start < other.start:
                    # self is before other
                    return other.start - self.end
                else:
                    # other is before self
                    return self.start - other.end
                

        ultra_conserved = []        

        start = None

        for i in range(len(self.consensus.seq)):
            # if no current region, check for start
            if start is None:
                
                # begin if consensus is ungapped
                if self.consensus.seq[i] != "-":
                    start = i
            
            # if there is a current region check for end
            else:
                if self.consensus.seq[i] == "-":
                    
                    # i - 1 because length ended at last position
                    # everything + 1 because the length is always 1 larger than the index difference
                    if ((i - 1) - start) + 1 >= min_len:

                        ultra_conserved.append(UltraConserved(start, i-1))

                    # reset the start position
                    start = None
            
        return ultra_conserved

    def find_potential_amplicons(self, ultra_conserved, amp_length, min_distance=100):
        """ 
        Locates regions between ultra_conserved regions and calculates the information content of the amplicon 
        min_distance is a computational param to limit the number of good matches to those likely to be long enough to be unique
        """


        amplicons = []

        # process each potential start (uc1) by each potential end (uc2)
        for uc1 in ultra_conserved:
            for uc2 in ultra_conserved:
                # only process each pairwise combination once using a lower-triangular algorithm
                if uc1 == uc2:
                    break

                # check if the UC regions are an appropriate distance apart
                distance = uc1.distance(uc2)

                if distance >= min_distance and distance <= amp_length:
                    if uc1.start < uc2.start:
                        start = uc1.end + 1
                        end = uc2.start - 1
                    else:
                        start = uc2.end + 1
                        end = uc1.start - 1

                    # add 1 to end because end is not inclusive 
                    score = sum(self.consensus.entropy[start:end + 1])

                    amp = PotentialAmplicon(start, end, score, self)
                    amp.bounding_left = uc2
                    amp.bounding_right = uc1
                    amplicons.append(amp)
                    
        return amplicons

def process_msa(msa_f, min_primer):
    cr = ConservedRegion(msa_f)

    ultra_conserved = cr.find_ultra_conserved(min_primer)
    [print(str(uc)) for uc in ultra_conserved]

    print()
    potential_amplicons = cr.find_potential_amplicons(ultra_conserved, 250)
    [print(str(amp) + "\t" + str(amp.get_uniqueness())) for amp in potential_amplicons]


def process_all_msa(msa_files, min_primer, amp_length):

    all_potential_amplicons = []
    num_processed = 0
    for msa_file in msa_files:
        cr = ConservedRegion(msa_file)

        # find potential primer regions
        ultra_conserved = cr.find_ultra_conserved(min_primer)

        # find amplicons between the primer regions
        potential_amplicons = cr.find_potential_amplicons(ultra_conserved, amp_length)

        all_potential_amplicons += potential_amplicons

        num_processed += 1

        print("Processed {} files.".format(num_processed), end="\r")

    print()


    sorted_amps = sorted(all_potential_amplicons, key=lambda amp: amp.score, reverse=True)
    
    print("{} Potential Amplicons Found".format(len(sorted_amps)))
    unique = 0
    with open("potential_amplicons.fasta", 'w') as OUT:
        for amp in sorted_amps[:500]:
            amp.print_mask(OUT)
            if amp.get_uniqueness():
                unique += 1

    print("{} Unique Amplicons Found".format(unique))

def uniqueness_histogram(msa_files):
    uniq_hash = {}
    num_processed = 0
    for msa_file in msa_files:
        cr = ConservedRegion(msa_file)

        uniq = cr.get_uniqueness()
        uniq_hash[uniq] = uniq_hash.get(uniq, 0) + 1

        num_processed += 1

        print("Processed {} files.".format(num_processed), end="\r")
    print()


    for key in sorted(uniq_hash):
        print(str(key) + "\t" + str(uniq_hash[key]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This script is meant to be used to find potential amplicons from a set of MSA files based on shannon entropy and a uniqueness value fround from a distance matrix that can be used to tell the group of aligned sequences apart but will amplify in them all.")
    parser.add_argument("-msa", help="multiple sequence alignment file", nargs="+")
    parser.add_argument("-amp_len", help="the maximum length of the amplicon", type=int)

    parser.add_argument("-min_primer", help="the minimum primer length (>= 10)", type=int, default=20)
    parser.add_argument("-max_primer", help="the maximum primer length (<= 40)", type=int, default=25)

    args = parser.parse_args()

    #uniqueness_histogram(args.msa)
    process_all_msa(args.msa, args.min_primer, args.amp_len)
    
