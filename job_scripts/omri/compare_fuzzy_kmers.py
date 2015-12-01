
import argparse
import os
import sys
import math
import time
from mypyli import jellyfish
import numpy as np
import logging
from Bio import SeqIO

import multiprocessing

import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})


logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel("DEBUG")

"""
Python allocates enough memory (28 bytes) in a standard int to hold 9 digits before it must be re allocated.

So, it would be good to be able to index genome, contig, and kmer in just 9 digits.

Want to allow 4 digits for genome index (each will be indexed as a number) to compare 1000+ genomes. 
Want 4 more for contig (Hopefully there are no more than 1000 contigs...)
Need to allow 8 digits for starting position to allow a 10Mb contig. 

A string can allow all 3 of these with ; between for 67 bytes.
"""


class KmerCounter(object):
    def __init__(self, mismatches=1, stable_start=2, stable_end=2):
        self.mismatches = mismatches

        # stables are # bases at beginning and end that cannot be variable
        self.stable_start = stable_start
        self.stable_end = stable_end

        # make partitions based on stable_start
        self.queues = {}
        self.partitions = {}
        self.results_queue = multiprocessing.Queue()
        for kmer in self._recurse_all_kmers(self.stable_start):
            self.queues[kmer] = multiprocessing.Queue()
            self.partitions[kmer] = KmerPartition(kmer, self.queues[kmer], self.results_queue)
            self.partitions[kmer].start()

        # set genome id to starting value
        self.genome_id = 0

    @classmethod
    def _recurse_all_kmers(cls, k, curseq="", collected_seqs=[]):
        """ Generates all kmers to avoid having to store a large array """
        if len(curseq) == k:
            yield curseq
        else:
            for nt in ["A", "C", "G", "T"]:
                for kmer in cls._recurse_all_kmers(k, curseq + nt, collected_seqs):
                    yield kmer
 
    def count_kmers(self, fasta_f, k):
        LOG.info("Counting {}-mers for {}".format(k, fasta_f))

        fasta_name = os.path.splitext(os.path.basename(fasta_f))[0]
        if os.path.isfile(fasta_name + ".counts.txt"):
            LOG.info("Found existing counts file for {}".format(fasta_name))
            kmer_counts = fasta_name + ".counts.txt"

        else:
            # count kmers using jellyfish
            jf_indx = jellyfish.count_kmers(fasta_f, k, prefix=fasta_name)
            kmer_counts = jellyfish.dump_kmer_counts(jf_indx, output_f=fasta_name + ".counts.txt")

        # read in the counts file
        num_processed = 0
        with open(kmer_counts, 'r') as IN:
            for line in IN:
                kmer, count = line.rstrip().split("\t")

                # send the kmer to the correct partition for processing
                self.partitions[kmer[:self.stable_start]].add_kmer(kmer, count)
                
                
                num_processed += 1
                #if num_processed == 100000:
                    #sys.exit()

        [print((par.name, len(par.kmers))) for par in self.partitions.values()]

        conserved_kmers = []
        for par in self.partitions.values():
            print((par.name, len(par.kmers)))
            par.update_kmer_file()
            conserved_kmers += par.get_conserved_kmers()

    def count_kmers_custom_single(self, fasta_f, k):
    
        self.genome_id += 1
        genome_str = str(self.genome_id)

        kmers_counted = 0

        contig_id = 0
        for seq in self._parse_fasta(fasta_f):
            contig_id += 1
            gen_ctg_str = ";".join((genome_str, str(contig_id)))

            # count kmers in a sliding window
            # must add 1 to the length because range has an exclusive end
            for indx in range(len(seq)+1 - k):
                kmers_counted += 1
                kmer = seq[indx:indx+k]

                # check if there is an N in the kmer and discard it if so
                if "N" in kmer:
                    continue
                
                self.partitions[kmer[:self.stable_start]].add_kmer(kmer, ";".join((gen_ctg_str, str(indx))))

                if kmers_counted % 10000 == 0:
                    print(kmers_counted)

        for partition in self.partitions.values():
            partition.update_kmer_file()
    
    def count_kmers_custom(self, fasta_f, k):
    
        self.genome_id += 1
        genome_str = str(self.genome_id)

        kmers_counted = 0

        contig_id = 0
        for seq in self._parse_fasta(fasta_f):
            contig_id += 1
            gen_ctg_str = ";".join((genome_str, str(contig_id)))

            # count kmers in a sliding window
            # must add 1 to the length because range has an exclusive end
            for indx in range(len(seq)+1 - k):
                kmers_counted += 1
                kmer = seq[indx:indx+k]

                # check if there is an N in the kmer and discard it if so
                if "N" in kmer:
                    continue
               
                self.partitions[kmer[:self.stable_start]].kmer_queue.put((kmer, ";".join((gen_ctg_str, str(indx))))) 

                if kmers_counted % 10000 == 0:
                    # check if the queues are getting too big, wait if so
                    for p in self.partitions.values():
                        if p.kmer_queue.qsize() > 10000:
                            time.sleep(2)

                    print(kmers_counted)

        for p in self.partitions.values():
            p.kmer_queue.put("DUMP")


    def get_conserved_kmers_single(self, num_genomes):
        conserved_kmers = []
        for p in self.partitions.values():
            conserved_kmers += p.get_conserved_kmers()

        return conserved_kmers

    def get_conserved_kmers(self, num_genomes):
        for p in self.partitions.values():
            p.kmer_queue.put("CONSERVED{}".format(num_genomes))

        for p in self.partitions.values():
            p.kmer_queue.close()
            p.kmer_queue.join_thread()
            p.join()

            yield self.results_queue.get()

        self.results_queue.close()


    def _parse_fasta(self, fasta_f):
        """ Yields sequences from a fasta file """

        with open(fasta_f, 'r') as IN:
            for record in SeqIO.parse(IN, 'fasta'):
                yield str(record.seq)


class KmerPartition(multiprocessing.Process):
    """ A free-standing class (for multithreading purposes) that processes a subset of kmers 
    
    Responsible for holding kmer information and dumping it to a file when it gets the signal.

    I think the best way to hold kmer information is in a trie. 
    """

    # binary here is merely idealistic because the binary gets converted to ints
    #encode_map = {"A": 0b00 "a": 0b00, "C": 0b01, "c": 0b01, "G": 0b10, "g": 0b10, "T": 0b11, "t":0b11}
    encode_map = {"A": "00", "C": "01", "G": "10", "T": "11"}


    def __init__(self, name, kmer_queue, results_queue):
        super().__init__()
        self.name = name
        self.kmer_queue = kmer_queue
        self.results_queue = results_queue

        # set up the storage file
        self.out_dir = "."
        self.file = "/".join([self.out_dir, self.name + "_data_dump.txt"])
        
        if not os.path.isfile(self.file):
            with open(self.file, 'w'):
                pass

        # set up some vars to use later
        self.kmers = {}
    
        self.mismatches = 1
        self.stable_start = 2
        self.stable_end = 2
        self.k = 20


    def run(self):
        print("Starting thread {}".format(self.name)) 
        # wait for kmers to count
        while True:
            itm = self.kmer_queue.get()

            if type(itm) is tuple:
                self.add_kmer(itm[0], itm[1])

            else:

                if itm == "DUMP":
                    print("{} got dump command".format(self.name))
                    self.expand_and_update_kmer_file()

                elif itm.startswith("CONSERVED"):
                    # split the number of genomes from the itm
                    self.get_conserved_kmers(int(itm.split("D")[1]))
                    break

        return

    @classmethod
    def _recurse_all_groups(cls, k, curseq="", mismatches=0):
        """ Generates all kmer groups for a given length k """
        if len(curseq) == k:
            yield curseq

        # don't allow N's in stable regions
        elif len(curseq) <= self.stable_start:
            pass

            for nt in ["N", "A", "C", "G", "T"]:
                for kmer in cls._recurse_all_kmers(k, curseq + nt):
                    yield kmer
 
    def _digest_kmer(self, kmer, current_seq="", mismatches=0):
        """ This is a slow point of my program right now. I need to speed up. """
        # check if we are to the end
        if len(current_seq) >= self.k - self.stable_end:
            # not sure why I can't just return rather than yield, then return
            # might have something to do with the fact that some will be yielded and others will be returned
            yield "".join((current_seq, kmer))
            return
        else:

            # check if the current base can be fuzzy (not too many previous errors)
            if mismatches < self.mismatches:
                # run the ambiguous result
                for result in self._digest_kmer(kmer[1:], current_seq+"N", mismatches+1):
                    yield result
                    break
 
            # base cannot be fuzzy, already enough errors
            else:
                # once again, not sure why I can't just return this...
                yield "".join((current_seq, kmer))
                return
        
        # run the next iteration if not returned
        for result in self._digest_kmer(kmer[1:], current_seq+kmer[0], mismatches):
            #print(("result", result))
            yield result

    def add_kmer(self, kmer, payload=None):
        """ Adds a kmer to the dict with a payload """
        for group in self._digest_kmer(kmer):
            try:
                self.kmers[group].append(payload)
            except KeyError:
                self.kmers[group] = [payload]

    def update_kmer_file(self):
        """ Loads/creates a kmer file and writes all the current kmer information to it. """

        LOG.debug("Updating kmer file for {}".format(self.file))

        tmp_file = self.file + "_tmp"

        # iterate through the input file and stored kmers and dump ordered results
        with open(self.file, 'r') as IN, open(tmp_file, 'w') as OUT:

            stored_line = IN.readline()[:-1]
           

            for kmer in sorted(self.kmers):
                # make sure something is on the line; assume file over if not
                if stored_line:
                    #stored_kmer = int(stored_line.split("\t")[0])
                    stored_kmer = stored_line.split("\t")[0]
                else:
                    line = "\t".join([str(kmer), "\t".join(self.kmers[kmer])])
                    OUT.write(line + "\n")
                    continue

                if kmer < stored_kmer:
                    # write kmer
                    line = "\t".join([str(kmer),"\t".join(self.kmers[kmer])])
                elif kmer == stored_kmer:
                    line = "\t".join([stored_line, "\t".join(self.kmers[kmer])])
                    stored_line = IN.readline()[:-1]

                else:
                    # read in some more lines until the line is greater than 
                    line = stored_line
                    stored_line = IN.readline()[:-1]

                OUT.write(line + "\n")

        os.rename(tmp_file, self.file)
        self.kmers = {}

    ###
    # Checking analyzing stored files
    ###
    def get_conserved_kmers(self, num_genomes):
        """ Look for kmers that have a position from each genomes """
       
        LOG.debug("Finding conserved kmers for {}...".format(self.name))
        conserved_kmers = []
        genomes_not_found = {}
        with open(self.file, 'r') as IN:
            for line in IN:
                elems = line[:-1].split("\t")
                kmer = elems[0]
                locations = elems[1:]
    
                # let's assume that most kmers will NOT be found in all genomes
                # so, I don't want to make objects (slow) until I know they are conserved. 
                # this will be slower for kmers that ARE in all genomes
    
                # first check if there are enough locations -- this will remove most kmers
                if len(locations) >= num_genomes:
                    
                    # now check that each genome is represented and build a location list
                    genomes_found = [False]*num_genomes
                    for location in locations:
                        try:
                            genome = int(location.split(";")[0]) - 1
                        except:
                            print(locations)
                            raise
                        try:
                            genomes_found[genome].append(location)
                        except AttributeError:
                            genomes_found[genome] = [location]
                        
                    not_found = []
                    for indx, pos in enumerate(genomes_found):
                        if pos is False:
                            not_found.append(indx)
                    
      
                    if not_found:
                        # build a histogram like data structure with kmers not found
                        try:
                            genomes_not_found[len(not_found)].append(kmer)
                        except KeyError:
                            genomes_not_found[len(not_found)] = [kmer]
                    else:
                        # if all genomes found add to list of conserved regions
                        conserved_kmers.append((kmer, locations))

        print("Found {} conserved kmers".format(len(conserved_kmers)))
        self.results_queue.put(conserved_kmers)
        return conserved_kmers


    def expand_and_update_kmer_file(self):
        """ Loads/creates a kmer file and writes all the current kmer information to it. """

        LOG.debug("Updating kmer file for {}".format(self.file))

        tmp_file = self.file + "_tmp"

        # iterate through the input file and stored kmers and dump ordered results
        with open(self.file, 'r') as IN, open(tmp_file, 'w') as OUT:

            stored_line = IN.readline()[:-1]
           

            for kmer in sorted(self.kmers):
                # make sure something is on the line; assume file over if not
                if stored_line:
                    #stored_kmer = int(stored_line.split("\t")[0])
                    stored_kmer = stored_line.split("\t")[0]
                else:
                    line = "\t".join([str(kmer), "\t".join(self.kmers[kmer])])
                    OUT.write(line + "\n")
                    continue

                if kmer < stored_kmer:
                    # write kmer and all expansions
                    line = "\t".join([str(kmer),"\t".join(self.kmers[kmer])])

                elif kmer == stored_kmer:
                    line = "\t".join([stored_line, "\t".join(self.kmers[kmer])])
                    stored_line = IN.readline()[:-1]

                else:
                    # read in some more lines until the line is greater than 
                    line = stored_line
                    stored_line = IN.readline()[:-1]

                OUT.write(line + "\n")

        os.rename(tmp_file, self.file)
        self.kmers = {}


class Amplicon(object):

    existing_amplicons = {}    

    def __init__(self, name):
        self.name = name
        self.locations = []

    @classmethod
    def from_kmer_pair(cls, k1, loc1, k2, loc2):
        key = tuple(sorted((k1, k2)))
        try:
            cls.existing_amplicons[key].add_new_location(loc1, loc2)
        except KeyError:
            amp = cls("-".join((k1, k2)))
            amp.add_new_location(loc1, loc2)
            cls.existing_amplicons[key] = amp

    def add_new_location(self, loc1, loc2):
       self.locations.append((loc1, loc2))




#
## Find kmers that are the right distance apart
#
def get_properly_spaced_pairs(conserved_kmers, k=20, max_dist=400, min_dist=100):
    ppp = {}
    for kmer, locations in conserved_kmers:
        for location in locations:
            genome, contig, start = location.split(";")
            
            try:
                ppp[(genome, contig)].append((kmer, genome, contig, start))
            except KeyError:
                ppp[(genome, contig)] = [(kmer, genome, contig, start)]

    # this part could be split up over multiple processes
    results_queue = multiprocessing.Queue()
    bins_processed = 0
    for data_bin in ppp.values():
        print("Processed {} bins".format(bins_processed))
        bins_processed += 1


        _find_properly_spaced(results_queue, data_bin, k, max_dist, min_dist)

    while True:
        if results_queue.qsize() > 0:
            yield results_queue.get()

        else:
            break

def _find_properly_spaced(queue, data_bin, k=20, max_dist=400, min_dist=100):
    # sort kmers by position
    sorted_kmers = sorted(data_bin, key=lambda k: k[-1])

    print("# kmers: {}".format(len(sorted_kmers)))

    for indx1, k1 in enumerate(sorted_kmers):
        # begin at the next kmer
        indx2 = indx1 + 1
        while indx2 < len(sorted_kmers):
            k2 = sorted_kmers[indx2]
            indx2 += 1

            # locations are starting values so need to add k to k1 to get the end of it
            distance = int(k2[-1]) - (int(k1[-1]) + k) 

            # break at first one that fails max dist test
            if distance > max_dist:
                break

            # skip if the distance is too small
            elif distance < min_dist:
                continue

            # distance must be acceptable
            else:
                loc1 = ";".join(k1[1:])
                loc2 = ";".join(k2[1:])

                # puts a tuple(first, second) of tuple(kmer, location)
                queue.put(((k1[0], loc1), (k2[0], loc2)))




def main(fastas, k):

    kcounter = KmerCounter()
   
    # count kmers
    for fasta in fastas:
        #kcounter.count_kmers_custom(fasta, 20)
        pass

    # map the conserved kmers back to genomes
    conserved_kmers = []
    for kmer_bundle in kcounter.get_conserved_kmers(len(fastas)):
        conserved_kmers += kmer_bundle

    # send kcounter for garbage collection
    del kcounter

    for pair in get_properly_spaced_pairs(conserved_kmers):
        start, end = pair

        s_kmer, s_loc = start
        e_kmer, e_loc = end

        Amplicon.from_kmer_pair(s_kmer, s_loc, e_kmer, e_loc)


    print("Properly spaced amplicons = {}".format(len(Amplicon.existing_amplicons)))

    """
    All that is left to do is to make a table of which genomes have unique sequences for each amplicon.

    I'm thinking there will be multiple files ideally. 

    1. Amplicon stats - sorted file that has # of Ns in fwd and rev primers and the number of genomes it can tell appart

    2. Amp X genome matrix - matrix that tells which genomes each amplicon can tell apart

    3. Fasta file for each amplicon that can tell apart any number of genomes uniquely (n > 1).
    """



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fastas", help="a bunch of fastas to compare", nargs="+", required=True)
    parser.add_argument("-k", help="value of k to use", required=True, type=int)
    args = parser.parse_args()

    main(args.fastas, args.k)
