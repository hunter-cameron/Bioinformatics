
import sys
import argparse
import logging
import os
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
import subprocess
import time
import math

from mypyli import samparser

logging.basicConfig()
LOG = logging.getLogger()
LOG.setLevel("DEBUG")

class QueryMapper(object):

    def __init__(self, query_f, out_dir=None, split_len=50, min_id=.90, min_len=40):
        self.query_f = query_f
        self.q_basename = os.path.splitext(os.path.basename(query_f))[0]


        # set up the outdir
        if out_dir is None:
            self.out_dir = self.q_basename
        else:
            self.out_dir = out_dir

        if not self.out_dir.endswith("/"):
            self.out_dir += "/"

        if not os.path.isdir(self.out_dir):
            os.mkdir(self.out_dir)

        self.split_len = split_len
        self.min_id = min_id
        self.min_len = min_len

   
        self.split_f = self.out_dir + self.q_basename + ".split.fasta"
        self.regions = []

    def process(self, ref_f):

        if not self.regions:
            # either  er split the query or read a previous split
            self.split_query()
        

        r_basename = os.path.splitext(os.path.basename(ref_f))[0]

        bbmap_f = self.out_dir + self.q_basename + "__--__" + r_basename + ".sam"

        # skip if this is done
        if not os.path.isfile(bbmap_f):
            LOG.info("Mapping {} to {}...".format(self.split_f, ref_f))
            self.run_bbmap(ref_f, bbmap_f)
        else:
            LOG.info("Skipping mapping {} to {}...".format(self.split_f, ref_f))

        # parse each hit in the SAM file
        alns_found = 0
        quality_alns = 0
        for record in self.parse_SAM_file(bbmap_f):
            alns_found += 1
            if record.mapped and record.perc_id >= self.min_id:
                qry_region = GenomeRegion.from_header(record.qname)
    
                if record.length >= self.min_len:
                    self.regions[qry_region.index].report_hit(r_basename, record.rname, record.pos, record.length, record.perc_id)
                    quality_alns += 1

        LOG.info("{}({} quality) hits among {} splits found in {}".format(alns_found, quality_alns, len(self.regions), r_basename))


    def split_query(self):
        """ 
        Splits query into fragments of around split_len

        TODO: Add a small algorithm to calculate the actual split_len based on the
        length of the contig to get it as close as possible to the supplied split_len

        Returns the total length of the sequences.
        """

        # read in existing split file if it exists
        if os.path.isfile(self.split_f):
            with open(self.split_f, 'r') as IN:
                for record in SeqIO.parse(IN, "fasta"):
                    
                    # make an object to store the data
                    self.regions.append(QueryRegion.from_header(record.description))
                    
            LOG.info("Loaded split sequence from {}".format(self.split_f))

        else:

            LOG.info("Splitting {} into pieces ~{}bp. Storing as {}".format(self.query_f, self.split_len, self.split_f))

            split_number = 0 
            with open(self.query_f, 'r') as IN, open(self.split_f, 'w') as OUT:
                for record in SeqIO.parse(IN, 'fasta'):
                    seq = record.seq
                    seq_len = len(seq)
                
                    if seq_len < self.split_len:
                        continue

                    # calculate split_len for this contig this should ensure all splits are close to 
                    # one another in length
                    leftovers = seq_len % self.split_len
                    extra_per = int(leftovers / int(seq_len / self.split_len)) + 1
                    split_len = self.split_len + extra_per

                    # split the sequence into chunks and write each to the split_f 
                    # prepending the split index to the header
                    for seq_indx in range(0, seq_len, split_len):
                        # split out the seq
                        sp_seq = seq[seq_indx:seq_indx+split_len]

                        # calc the actual length of the seq
                        sp_seq_len = len(sp_seq)
                        
                        # make an QueryRegion object and add it to the list
                        region = QueryRegion(self.q_basename, record.description, seq_indx, seq_indx + sp_seq_len, split_number)
                        self.regions.append(region)

                        # write the seq
                        seqobj = SeqRecord(sp_seq, id=region.to_header(), description="")
                        SeqIO.write(seqobj, OUT, 'fasta')
                       
                        split_number += 1


    @staticmethod
    def param_dict_to_header(params):
        """ Converts a param dict to a header """

        header = "genome={base}|index={indx}|contig={contig}|slice=[{start}:{end}]".format(
                base=params["genome"],
                indx=str(params["index"]),
                contig=params["contig"],
                start=str(params["start"]),
                end=str(params["end"])
                )
        return header

    @staticmethod
    def header_to_param_dict(header):
        """ Converts a header to a param dict 
        
        header = "genome={base}|index={indx}|contig={desc}|slice=[{start}:{end}]"

        Uses brute force to ensure correct results
        """


        rest, slice = header.rsplit("|slice=", 1)
        rest, contig = rest.rsplit("|contig=", 1)
        rest, index = rest.rsplit("|index=", 1)
        genome = rest.split("genome=")[-1]

        start, end = slice[1:-1].split(":")

        return {"start": int(start), "end": int(end), "index": int(index), "contig": contig}
      
    def run_bbmap(self, reference, out):
        """ 
        Runs bbmap and returns a path of the output sam.
        Prints all top alignments for ambiguously mapped reads.
        """
        cpus = 8
        if self.split_len <= 500:
            prog = "bbmap.sh"
        else:
            prog = "mapPacBio.sh"

        #cmd = "{} ref={ref} in={split_f} local=t ambig=all nodisk overwrite=t sam=1.4 threads={cpus} out={out}".format(
        cmd = "{} ref={ref} in={split_f} local=t sssr=.01 secondary=t ambig=all maxsites=30 minid=.20 nodisk overwrite=t sam=1.4 threads={cpus} out={out}".format(
                prog,
                ref=reference,
                split_f=self.split_f,
                cpus=cpus,
                out=out)

        print("Running:\n    {}".format(cmd), file=sys.stderr)
        #code = subprocess.call(cmd.split(" "))         # safe way
    
        code = subprocess.call(cmd + " 2>>bbmap.err", shell=True)         # dangerous way

        if code:
            raise Exception("The bbmap command failed")
        else:
            return out

    @staticmethod
    def parse_SAM_file(sam_f):
        """ Yields SamRecords that meet some minimum criteria"""
        with open(sam_f, 'r') as IN:
            for record in samparser.parse(IN, mapq=0, aligned_only=True):
                yield record


class AmpliconAligner(object):

    @classmethod
    def process(cls, fasta_f, out_dir):

        # make directory if it doesn't exist
        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        align_f = out_dir + "/" + os.path.splitext(os.path.basename(fasta_f))[0] + ".afa"

        if os.path.isfile(align_f):
            LOG.info("Found alignment file '{}'. Skipping alignment step.".format(align_f))
            return align_f

        cls.run_muscle(fasta_f, align_f)
    
        return align_f

    @staticmethod
    def run_muscle(fasta_f, output_f):
        muscle_cmd = "muscle -in {} -out {}".format(fasta_f, output_f)

        bsub_cmd = "bsub -o muscle.out -e muscle.err -J auto_aln -n 2 -q week"
        
        cmd = bsub_cmd + " " + muscle_cmd

        code = subprocess.call(cmd, shell=True)         # dangerous way

        if code:
            raise Exception("The muscle command failed")

    @classmethod
    def wait_for_job(cls, job_name="auto_aln"):
        """ waits for job to complete, checks every 10 seconds """
        while cls._job_running(job_name):
            time.sleep(10)

    @staticmethod
    def _job_running(job_name="auto_aln"):

        output = subprocess.check_output([
                    "bjobs",
                    "-J", job_name
                ])
        #print(output)
        if output:
            return True
        else:
            return False


class GenomeRegion(object):
    """ 
    Basic information for accessing a region in a genome 
    
    Start and end are intended to be indices for extracting the region from a seq (base 0, half closed)
    """
 
    def __init__(self, genome, contig, start, end, index):

        self.genome = str(genome)
        self.contig = str(contig)
        self.start = int(start)
        self.end = int(end)
        self.index = int(index)

    def __str__(self):
        return self.to_header()

    @property
    def length(self):
        return self.end - self.start

    @classmethod
    def from_header(cls, header):
        """ Convert from the header in a FASTA file to a GenomeRegion """

        rest, slice = header.rsplit("|slice=", 1)
        rest, contig = rest.rsplit("|contig=", 1)
        rest, index = rest.rsplit("|index=", 1)
        genome = rest.split("genome=")[-1]

        start, end = slice[1:-1].split(":")

        return cls(genome, contig, start, end, index)

    def to_header(self):
        """ Convert a GenomeRegion to a FASTA header """
        
        return "genome={self.genome}|index={self.index}|contig={self.contig}|slice=[{self.start}:{self.end}]".format(self=self)


class QueryRegion(GenomeRegion):
    """ A specific GenomeRegion of the Query. """

    def __init__(self, genome, contig, start, end, index):
        super().__init__(genome, contig, start, end, index)

        # this is a dict of GenomeRegions (indexed by genome) 
        self.alignments = {}

    @classmethod
    def from_header(cls, header):
        """ An override method to init the right object """

        rest, slice = header.rsplit("|slice=", 1)
        rest, contig = rest.rsplit("|contig=", 1)
        rest, index = rest.rsplit("|index=", 1)
        genome = rest.split("genome=")[-1]

        start, end = slice[1:-1].split(":")

        return cls(genome, contig, start, end, index)

    def report_hit(self, genome, contig, start, length, perc_id):
        """ Adds a new alignment to the QueryRegion """
        
        ref_region = GenomeRegion(genome, contig, start - 1, (start - 1) + length, self.index) 

        try:
            self.alignments[genome].append(ref_region)
        except KeyError:
            self.alignments[genome] = [ref_region]

    def has_hits_from(self, genomes):
        """ Checks if alignment has hits from a list of genomes. Returns True is so, else False """

        for genome in genomes:
            if genome not in self.alignments:
                return False

        return True


class Amplicon(object):
   
    ampliconIndex = 0       # this is for indexing amplicons

    def __init__(self, query1, query2, max_len, min_len=20):
        self.q1, self.q2 = self._order_queries(query1, query2)

        self.amp_index = self.ampliconIndex
        Amplicon.ampliconIndex += 1

        if self.q1.contig != self.q2.contig or self.q2.end - self.q1.start > max_len:
            raise ValueError("{} and {} are not spaced an appropriate distance apart.".format(self.q1, self.q2))

        self.min_len = min_len
        self.max_len = max_len

    @staticmethod
    def _order_queries(q1, q2):
        # three possible options q1 is before, after, or overlaps q2

        # first determine if they overlap keeping in mind the end is half closed so the real end is -1
        # so as a twist on the usual expression we use < rather <=
        # really, because of the way I built these seqs, they should never overlap...
        if q1.start < q2.end and q2.start < q1.end:
            raise ValueError("{} and {} overlap.".format(q1, q2))

        # now that we know there is no messy over lap, we need to see which is first
        elif q1.start < q2.start:
            # q1 is before q2
            return q1, q2
        else:
            # q2 is before q1
            return q2, q1

    def check_all_distance(self, return_failed=False):
        """ Checks if query1 and query2 are spaced an appropriate distance apart in all genomes. True if so, else false """

        failed = []
        for genome in self.q1.alignments:
            try:
                self.get_correctly_spaced(genome)
            except ValueError:
                if return_failed:
                    failed.append(genome)
                else:
                    return False
       
        if return_failed:
            if failed:
                return False, failed
            else:
                return True, failed
        else:
            return True

    def get_correctly_spaced(self, genome):
        """ Returns the first correctly spaced GenomeRegion pair or raises ValueError """

        for aln1 in self.q1.alignments[genome]:
            for aln2 in self.q2.alignments[genome]:
                if aln1.contig == aln2.contig:
                    if aln2.start - aln1.end >= self.min_len and aln2.start - aln1.end <= self.max_len:
                        # found one pair, return
                        return aln1, aln2

        raise ValueError("No alignments meet the criteria.")
 
    def get_seq_shopping_list(self):
        """ Returns a list of GenomeRegion objects that make up this amplicon """

        shop_list = []

        # add the query
        qry_shop = GenomeRegion(self.q1.genome, self.q1.contig, self.q1.start, self.q2.end, self.amp_index)

        shop_list.append(qry_shop)


        # add all the others
        for genome in self.q1.alignments:
            g1, g2 = self.get_correctly_spaced(genome)
            shop_list.append(GenomeRegion(g1.genome, g1.contig, g1.start, g2.end, self.amp_index))

        return shop_list


class AlignedAmplicon(object):
    """ NOTE: THIS CLASS HAS NO PROGRAM RELATION TO Amplicon. THE NAMES ARE SIMILAR BY BIOLOGY ONLY"""

    def __init__(self, msa_f, min_prim, name=None):
        self.msa_f = msa_f
        self.msa = self._read_msa(msa_f)

        self.min_prim = min_prim

        if name:
            self.name = name
        else:
            self.name = os.path.splitext(os.path.basename(msa_f))[0]

        self.consensus = Consensus(self.msa)

        # params relating to uniqueness
        self.max_uniqueness = None
        self.not_unique = None
        self.min_unique_start = None
        self.min_unique_end = None

    @classmethod
    def print_amplicon_summaries(cls, amplicons, out_dir):

        if not os.path.isdir(out_dir):
            os.mkdir(out_dir)

        # calc uniqueness stats
        for amplicon in amplicons:
            amplicon.set_minimum_unique_region()
        
        with open(out_dir + "/" + "amplicon_stats.txt", 'w') as STATS, open(out_dir + "/" + "amplicon_masks.txt", 'w') as MASKS:
            STATS.write("\t".join(["amplicon", "uniqueness", "num_not_unique", "not_unique_groups"]) + "\n")

            for amplicon in sorted(amplicons, key=lambda amp: (amp.max_uniqueness, -amp.num_not_unique), reverse=True):

                if amplicon.max_uniqueness:
                    stats_list = [amplicon.name, str(amplicon.max_uniqueness), str(amplicon.num_not_unique), "None"]
                    amplicon.write_mask(MASKS)
            
                else:
                    stats_list = [amplicon.name, str(amplicon.max_uniqueness), str(amplicon.num_not_unique), "\t".join([",".join(group) for group in amplicon.not_unique])]



                STATS.write("\t".join(stats_list) + "\n")


    @property
    def num_not_unique(self):
        if self.max_uniqueness is None:
            raise AttributeError("get_uniqueness() has not been ran.")

        if self.not_unique is None:
            return 0
        else:
            count = 0
            for group in self.not_unique:
                count += len(group)

            return count


    def get_uniqueness(self, start, end, return_nodist=False):
        """ Makes a distance matrix using the MSA and returns the lowest value """
        # arbitrarily high starting number
        min_distance = 99999

        # calculate the distance using lower triangular algorithm
        print()
        print()
        #print("Here1")
        seqs = self.msa[:, start:end]
        no_dist = []
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

                if min_distance == 0:
                    if return_nodist:
                        new_no_dist = []
                        s1_group = None
                        s2_group = None
                        for group in no_dist:
                            if seq1.description in group:
                                s1_group = group

                            if seq2.description in group:
                                s2_group = group
                           
                            if seq1.description not in group and seq2.description not in group:
                                new_no_dist.append(group)

                        # check exactly what the additional group should be
                        if s1_group and s2_group is None:
                            print("h1")
                            s1_group.append(seq2.description)
                            new_no_dist.append(s1_group)
                        elif s1_group is None and s2_group:
                            print("h2")
                            s2_group.append(seq1.description)
                            new_no_dist.append(s2_group)
                        elif s1_group and s2_group:
                            # check if it is the same group
                            if s1_group is s2_group:
                                continue
                            print("h3")
                            group = set(s1_group + s2_group)
                            new_no_dist.append(list(group))
                        else:       # neither is in a group
                            print("h4")
                            new_no_dist.append([seq1.description, seq2.description])
                        
                        no_dist = new_no_dist
                        
                    else:
                        # if min = 0 we know that the uniqueness is 0 so no need to process any more
                        return 0
        
        if return_nodist:
            return min_distance, no_dist
        else:
            return min_distance
        
    def write_mask(self, fh):

        lprim = self.consensus.seq[0:self.min_unique_start]
        rprim = self.consensus.seq[self.min_unique_start:self.min_unique_end]
        num_Ns = self.min_unique_end - self.min_unique_start

        mask = lprim + "N"*num_Ns + rprim

        header = ">{self.name} [{self.min_unique_start}:{self.min_unique_end}] max_unique={self.max_uniqueness}".format(self=self)

        fh.write(header + "\n" + mask + "\n")
        
    def set_minimum_unique_region(self):
        """ 
        Sets a minimum unique region to fh where minimum unique mask is flanking regions as long as possible with a center that is sufficient to differentiate each seq from one another
        
        The algorithm:

        1. Split the seq into two parts at the midpoint. Start of left is start + primer. End of left is midpoint -1.
            Start of right is midpoint end of right is end - primer - 1
        
        :: Begin binary search
        2. Get midpoint of both sides from start to end
        3. Check uniqueness in inter region
        4. If not unique, set midpoint to end of left and beginning of right and goto 2
        4. If unique, set beginning of left to midpoint 
        
        """

        seq_len = len(self.msa[0])
        midpoint = int(seq_len / 2)     # get an approximate midpoint

        lstart = self.min_prim
        lend = midpoint
        rstart = midpoint
        rend = seq_len - self.min_prim - 1      # - 1 for half open index

        # first check if this amplicon can possibly be unique
        uniqueness, no_dist = self.get_uniqueness(lstart, rend, return_nodist=True)
        
        self.max_uniqueness = uniqueness
        if uniqueness == 0:
            self.not_unique = no_dist
            return
       
        while True:
            # get midpoints
            lmid = int((lend - lstart) / 2) + lstart
            rmid = int((rend - rstart) / 2) + rstart

            # test if range is unique
            if self.get_uniqueness(lmid, rmid):
                # if we are within 5 (so within 10 for both front and back) stop
                if lmid - lstart <= 5:
                    self.min_unique_start = lmid
                    self.min_unique_end = rmid
                    return
                else:
                    # otherwise make the range more narrow
                    lstart = lmid
                    rend = rmid
            else:
                # make the range more broad
                lend = lmid
                rstart = rmid

    @staticmethod
    def _read_msa(msa_f, aln_format="fasta"):
        return AlignIO.read(msa_f, aln_format)


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


def write_seqs_from_shopping_list(ref_files, shopping_list, out_dir):
    """ Basically just writes a bunch of seqs to FASTA files. Does a lot of stuff to make this go as quickly as possible. Returns a list of written paths. """

    # make directory if it doesn't exist
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    # refactor shopping list by index (output file)
    file_map = {}
    for region in shopping_list:
        try:
            file_map[region.index].append(region)
        except KeyError:
            file_map[region.index] = [region] 

    
    # check if any files can skipped because results from a previous run were found
    written_fastas = []
    for basen in list(file_map.keys()):
        fasta_name = "{}/amplicon{}.fna".format(out_dir, basen)
        if os.path.isfile(fasta_name):
            del file_map[basen]
            written_fastas.append(fasta_name)

    # check if any fasta were found
    if written_fastas:
        # check if all fastas were found
        if not file_map:
            LOG.info("Found directory that contained all amplicons already written. Skipping writing...")
            return written_fastas
        else: 
            LOG.info("Found directory that contained already written potential amplicons. Reused these where possible.")

    #
    ## Calculate how many files should be stored per iteration
    #

    maximum_per_set = 100000
    
    # all files should have the same number of seqs
    seqs_per_file = len(list(file_map.values())[0])
    files_per_set = int(maximum_per_set / seqs_per_file)

    #
    ## get each iteration of seqs
    # 
    files = list(file_map.keys())
    for indx in range(0, len(files), files_per_set):

        # make a dict of what I need to get from each FASTA
        to_get = {}
        for amp_f in files[indx:indx + files_per_set]:
            for region in file_map[amp_f]:

                try:
                    to_get[region.genome][region.contig].append(region)
                except KeyError:
                    try:
                        to_get[region.genome][region.contig] = [region]
                    except KeyError:
                        to_get[region.genome] = {region.contig: [region]}


        # Actually get the seqs from each file
        seqs_to_write = {}
        for ref in ref_files:
            r_basename = os.path.splitext(os.path.basename(ref))[0]
           
            # make sure there are seqs to get from this ref
            if r_basename not in to_get:
                raise ValueError("No sequences to get for {}".format(ref))

            with open(ref, 'r') as IN:
                for record in SeqIO.parse(IN, "fasta"):
                    # check if this contig has seqeuences to be extracted
                    if record.description not in to_get[r_basename]:
                        continue
                   
                    # make the seqobj and add it to be written
                    for region in to_get[r_basename][record.description]:
                        seqobj = record[region.start:region.end]
                        seqobj.id = region.to_header()
                        seqobj.description = ""

                        try:
                            seqs_to_write[region.index].append(seqobj)
                        except KeyError:
                            seqs_to_write[region.index] = [seqobj]

                    # report all seqs from the contig written
                    to_get[r_basename][record.description] = True

        # check all seqs were successfully found
        for ref in to_get:
            for contig  in to_get[ref]:
                if to_get[ref][contig] is False:
                    LOG.warning("Not all regions extracted from '{}' contig '{}'".format(ref, contig))

        # write all seqs
        for basen in seqs_to_write:
            fasta_name = "{}/amplicon{}.fna".format(out_dir, basen)
            written_fastas.append(fasta_name)
            with open(fasta_name, 'w') as OUT:
                SeqIO.write(seqs_to_write[basen], OUT, "fasta")

    return written_fastas

def main(args):

    # see if the user has a specific file they want for the query
    if args.qry:
        LOG.debug("Using {} as the query FASTA".format(args.qry))
        qry = args.qry
    
    else:
        LOG.debug("Using {} as the query FASTA".format(args.i[0]))
        qry = args.i[0]

    references = [genome for genome in args.i if genome != qry]
    

    #
    ## Map the query to each reference
    #

    qm = QueryMapper(qry, split_len=args.split_len, min_id=args.min_id, min_len=args.min_len)
    for ref in references:
        qm.process(ref)

    qry_regions = qm.regions

    #
    ## find seqs that were mapped in all 
    #

    LOG.info("Looking for splits that had alignments in all references...")
    ref_names = [os.path.splitext(os.path.basename(ref))[0] for ref in references]
    mapped_all = []
    for region in qry_regions:
        if region.has_hits_from(ref_names):
            mapped_all.append(region)

    LOG.info("{} of {} splits had alignments in all references.".format(str(len(mapped_all)), str(len(qry_regions))))


    if not mapped_all:
        LOG.warning("No splits had alignments in all references. Aborting.")
        return

    #
    ## Make and check amplicons
    #
    
    """
    Strategy:

    1. Loop through all alignments by index. 
    2. To get the second alignment in the pair begin looping with index + 1
    3. Break at the first bad pair because the pairs are ordered

    This should allow me to capture each pair only once and not process unnecessary pairs.
    """

    LOG.info("Making and checking potential amplicons...")

    amplicons = []

    # these will be used for some reporting stats
    failed_amplicons = {}
    total_amplicons = 0
    
    num_aln = len(mapped_all)
    for indx1 in range(num_aln):

        # report status of this operation
        if indx1 % 1000 == 0:
            print("{}% completed.".format(indx1 / num_aln * 100), end="\r")

        qry1 = mapped_all[indx1]
        indx2 = indx1 + 1
       

        # set up infinite loop
        while True:
            # check if index is out of bounds
            if indx2 == num_aln:
                break
            
            total_amplicons += 1 

            qry2 = mapped_all[indx2]
        
            # make the Amplicon, exception if query is inappripriate distance apart
            try:
                amp = Amplicon(qry1, qry2, args.amp_len)
            except ValueError:
                # after the first one fails, because regions are ordered, we know the rest will fail
                break

            # add the amplicon if it checks out
            correct_dist, failed = amp.check_all_distance(return_failed=True)
            if correct_dist:
                amplicons.append(amp)
            else:
                # add another to the failed column
                for genome in failed:
                    try:
                        failed_amplicons[amp].append(genome)
                    except KeyError:
                        failed_amplicons[amp] = [genome]

            indx2 += 1

    LOG.info("{} potential amplicons found.".format(len(amplicons)))

    # check if we need to go into amplicon report mode
    if len(amplicons) < 15:
        LOG.info("Low number of potential amplicons detected (of {} total), printing some statistics.".format(total_amplicons))
        num_failed_per_genome = {}
        num_failing_genomes = {}
        for amp in failed_amplicons:
            for genome in failed_amplicons[amp]:
                num_failed_per_genome[genome] = num_failed_per_genome.get(genome, 0) + 1

            num_failing_genomes[len(failed_amplicons[amp])] = num_failing_genomes.get(len(failed_amplicons[amp]), 0) + 1

        # print the stats
        LOG.info("Failing Amplicon Histogram:")
        LOG.info("  # fail\tcount")
        for num in sorted(num_failing_genomes, reverse=True):
            LOG.info("  {num}\t{count}".format(num=num, count=num_failing_genomes[num]))

        LOG.info("Failing Genome Histogram:")
        LOG.info("  genome\t# fail")
        for genome in sorted(num_failed_per_genome, key=lambda k: num_failed_per_genome[k], reverse=True):
            LOG.info("  {genome}\t{count}".format(genome=genome, count=num_failed_per_genome[genome]))

    if not amplicons:
        LOG.warning("Aborting due to lack of amplicons.")
        return

    
    #
    ## Write the amplicons 
    #

    LOG.info("Writting potential amplicons to FASTA files...")
    shopping_list = []
    for amp in amplicons:
        shopping_list += amp.get_seq_shopping_list()

    written_fastas = write_seqs_from_shopping_list(args.i, shopping_list, out_dir=args.out_dir + "/conserved_regions")


    #
    ## Align the written fastas
    #

    LOG.info("Aligning written FASTAs using MUSCLE...")

    aln_dir = args.out_dir + "/conserved_regions_aligned"
    aligned_files = []
    for fasta in written_fastas:
        aln_file = AmpliconAligner.process(fasta, aln_dir)
        aligned_files.append(aln_file)

    AmpliconAligner.wait_for_job()

    # a short sleep to give the system time to write the MSA
    time.sleep(20)

    #
    ## Write the final amplicon info
    #
    aligned_amps = []
    for aln_f in aligned_files:
        aligned_amps.append(AlignedAmplicon(aln_f, args.min_prim))

    AlignedAmplicon.print_amplicon_summaries(aligned_amps, args.out_dir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates ANI (average nucleotide identity) and AF (alignment fraction) between two genomes. This is done by breaking up each query genome into pieces and mapping those pieces back to the other genome in the pair. AF is the fraction of pieces that map back. ANI is the average identity of the pieces that map back. This script outputs a matrix where each cell is a tuple (ANI, AF). Currently, ANI is a percent and AF is a decimal but I should change them to make it a consistent measure before I put this in the scripts for others directory.")

    parser.add_argument("-i", help="bunch of fasta files to all-by-all compare", nargs="+", required=True)
    parser.add_argument("-qry", help="the file to be broken up and used as the query (if not the first one in '-i'")
    parser.add_argument("-split_len", help="target length to split query into", type=int, default=50)
    parser.add_argument("-min_id", help="minumum id for an alignment", type=float, default=.90)
    parser.add_argument("-min_len", help="minumum length for an alignment", type=int, default=40)
    parser.add_argument("-out_dir", help="directory to store output", default=os.getcwd())
    parser.add_argument("-amp_len", help="maximum amplicon length [%(default)s]", type=int, default=400)
    parser.add_argument("-min_prim", help="minumum primer length", default=20, type=int)
    args = parser.parse_args()

    main(args)
