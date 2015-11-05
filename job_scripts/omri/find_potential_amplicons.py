
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

    def find_consensus(self, max_entropy=100):
        consensus = ""
        entropy = []
        
        # process each column of the MSA
        for i in range(len(self.msa[0])):
            col = self.msa[:, i]
          
            # get the most frequent character and the entropy (I should probably decouple this...)
            consensus_nt = self.get_ambiguous_code(col)
            ent = self.shannon_entropy(col)
            
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
        """ Returns the shannon entropy for a given seq """

        # entropy running sum
        entropy = 0

        # store seq_len to avoid repeated calculations
        seq_len = len(seq)
        for nt in set(seq):
            prob = seq.count(nt) / seq_len

            # add the entropy weight to the running sum
            entropy += prob * math.log(prob, 2)

        return entropy * -1

    @staticmethod
    def get_ambiguous_code(seq):
        """ Returns the ambiguous code for a given seq """

        code_dict = {
                # basic nt
                'A': 'A',
                'C': 'C',
                'G': 'G',
                'T': 'T',

                'AC': 'M',
                'AG': 'R',
                'AT': 'W',

                'CG': 'S',
                'CT': 'Y',

                'GT': 'K',
                
                'ACG': 'V',
                'ACT': 'H',
                'AGT': 'D',
                'CGT': 'B',

                'ACGT': 'N'
                }

        bases = set(seq)

        key = "".join(sorted(bases))

        if "-" in key:
            return "-"
        else:
            return code_dict[key]


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

    def get_indistinguishable(self):
        """ 
        Returns a list of lists of seqs that are indistinguishable at this amplicon. 
        
        Each list within a list are genomes that have the same seq.
        """

        # get only the portion of the seq corresponding to this amplicon
        seqs = self.parent.get_region(self.start, self.end)

        no_dist = []
        # calculate the distance using lower triangular algorithm
        for seq1 in seqs:
            for seq2 in seqs:

                # lower triangular algorithm
                if seq1 == seq2:
                    break

                # calculate the distance between the two seqs
                distance = 0
                for nt1, nt2 in zip(seq1, seq2):
                    if nt1 != nt2:
                        distance += 1

                # add the seqs to a group if dist = 0 
                if distance == 0:
                    for group in no_dist:
                        if seq1.id in group:
                            if seq2.id not in group:
                                group.append(seq2.id)
                            break
                        elif seq2.id in group:
                            if seq1.id not in group:
                                group.append(seq1.id)
                            break
                    
                    # if didn't break out, must need to make a new group
                    else:
                        no_dist.append([seq1.id, seq2.id])

        return no_dist

    def print_mask(self, fh=sys.stdout):
        """ Prints a "mask" of the amplicon where the bounding left and right regions are the only NT printed and the variable region is represented by Ns """

        # + 1 to get len from distance
        num_Ns = (self.end - self.start) + 1

        fwd_bound = self.parent.get_region(self.bounding_left.start, self.bounding_left.end)[0].seq
        rev_bound = self.parent.get_region(self.bounding_right.start, self.bounding_right.end)[0].seq

        mask = str(fwd_bound) + "N"*num_Ns + str(rev_bound)

        header = ">" + "PotentialAmplicon {} {}..{}  score: {}".format(self.parent.name, self.bounding_left.start, self.bounding_right.end, self.score) + "  uniqueness: {}".format(self.get_uniqueness())

        fh.write(header + "\n" + mask + "\n")
        
    def print_bad_amp_summary(self, fh=sys.stdout):
        """ Prints a line detailing which genomes were not unique for the amplicon """
        
        name = "PotentialAmplicon {} {}..{} score: {}".format(self.parent.name, self.bounding_left.start, self.bounding_right.end, self.score)

        no_dist = self.get_indistinguishable()

        genome_count = 0
        genomes = ""
        for group in no_dist:
            genome_count += len(group)

            if genomes:
                genomes += " " + ",".join(group)
            else:
                genomes = ",".join(group)
        

        fh.write("\t".join([name, str(genome_count), str(len(no_dist)), genomes]) + "\n")



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
        self.msa = self._read_msa(msa_f)

        if name:
            self.name = name
        else:
            self.name = os.path.splitext(os.path.basename(msa_f))[0]

        self.consensus = Consensus(self.msa)


    @staticmethod
    def _read_msa(msa_f, aln_format="fasta"):
        return AlignIO.read(msa_f, aln_format)

    def filter_msa(self, genomes):
        """ Filters the MSA to only include the genomes supplied."""
        records = []
        for record in self.msa:
            header = record.description
            
            # name in second set of |, surrounded with brackets
            name = header.split("|")[1][1:-1]

            if name in genomes:
                records.append(record)

        self.msa = MultipleSeqAlignment(records)
        self.consensus = Consensus(self.msa)

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

    def find_ambig_ultra_conserved(self, min_len, max_ambig):
        """ 
        Finds ultra conserved regions 
      
        Strategy:

        1. Start a potential UC region at the beginning and after each bad character.
        2. When a bad character is found, check if we need to terminate potential UC.
        3. Add terminated UC regions that are long enough to finished 
        4. Return all the finished UC regions

        """


        class UltraConserved(object):
            """ Simple class to hold info about UC regions """
            def __init__(self, start, end=0, ambiguities=0):
                self.start = start
                self.end = end
                self.ambiguities = ambiguities

            def __str__(self):
                return "UltraConserved Region {self.start}..{self.end}".format(self=self)

            @property
            def length(self):
                # +1 to get length from distance
                return (self.end - self.start) + 1
        
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

        print(self.consensus.seq)

        uc_finished = []
        uc_in_progress = []
        last_was_terminal = True        # stores whether the last character was a possible terminator
        for i in range(len(self.consensus.seq)):
            # check for gaps, if found, kill all uc
            if self.consensus.seq[i] == "-":
                last_was_terminal = True
                
                for uc in uc_in_progress:
                    uc.end = i - 1

                    if uc.length >= min_len:
                        uc_finished.append(uc)

                uc_in_progress = []
                    
            # check for ambiguous code
            elif self.consensus.seq[i] in ["R", "Y", "K", "M", "S", "W", "B", "D", "H", "V", "N"]:
                last_was_terminal = True

                # update regions in progress
                next_in_progress = []
                for uc in uc_in_progress:
                    if uc.ambiguities < max_ambig:
                        uc.ambiguities += 1
                        next_in_progress.append(uc)

                    else:   # maximum ambig already reached
                        uc.end = i - 1
                        
                        # add to final if long enough
                        if uc.length >= min_len:
                            uc_finished.append(uc)

                uc_in_progress = next_in_progress

            # not gap or ambig; must be a good position
            else:
                if last_was_terminal:
                    last_was_terminal = False
                    uc_in_progress.append(UltraConserved(start=i))

        # process any open regions
        for uc in uc_in_progress:
            uc.end = i
            if uc.length >= min_len:
                uc_finished.append(uc)

        return uc_finished

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


def process_all_msa(msa_files, min_primer, primer_ambig, amp_length, subset=None):

    all_potential_amplicons = []
    num_processed = 0
    for msa_file in msa_files:
        cr = ConservedRegion(msa_file)
        
        if subset:
            genomes = []
            with open(subset, 'r') as IN:
                for line in IN:
                    genomes.append(line.strip())
            cr.filter_msa(genomes)

        # find potential primer regions
        ultra_conserved = cr.find_ambig_ultra_conserved(min_primer, primer_ambig)
        print("{} ultra conserved regions found for {}".format(len(ultra_conserved), msa_file))
        for uc in ultra_conserved[:10]:
            print(str(uc))

        if ultra_conserved:

            # find amplicons between the primer regions
            potential_amplicons = cr.find_potential_amplicons(ultra_conserved, amp_length)

            all_potential_amplicons += potential_amplicons

        num_processed += 1
        print("Processed {} files.".format(num_processed), end="\n")

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


    # if no unique amps, we want to print a "bad amplicon summary" to let the user know how close they were
    with open("bad_amplicon_summary.txt", 'w') as OUT:
        OUT.write("\t".join(["amplicon", "# genomes indistinguishable", "# groups", "genomes"]) + "\n")
        for amp in sorted_amps[:500]:
            amp.print_bad_amp_summary(OUT)

    print("Wrote bad amplicon summary to 'bad_amplicon_summary.txt'")

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
    parser.add_argument("-ambig", help="maximum number of ambiguous bases to allow in primers", type=int, default=0)
    parser.add_argument("-min_primer", help="the minimum primer length (>= 10)", type=int, default=20)
    #parser.add_argument("-max_primer", help="the maximum primer length (<= 40)", type=int, default=25)
    parser.add_argument("-subset", help="filter MSA to only include these genome names as listed in the GenBank file")

    args = parser.parse_args()

    #uniqueness_histogram(args.msa)
    process_all_msa(args.msa, args.min_primer, args.ambig, args.amp_len, args.subset)
    
