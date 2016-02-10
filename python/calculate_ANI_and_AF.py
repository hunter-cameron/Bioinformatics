
import sys
import argparse
import logging
import os
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
import subprocess

from mypyli import samparser, parallelizer

logging.basicConfig()
LOG = logging.getLogger()
LOG.setLevel("INFO")


class QueryMapper(object):

    def __init__(self, query_f, out_dir=None, split_len=1000, min_id=.70, min_cov=.70):
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
        self.min_cov = min_cov

   
        self.split_f = self.out_dir + self.q_basename + ".split.fasta"
        self.aligns = None


        self.parallel_id = None
        self.mapped_files = {}

        self._setup()

    def _setup(self):
        """ Sets up for processing by either splitting the query or loading the splits from a previous run"""
        # either split the query or set the list from the last time to None
        if os.path.isfile(self.split_f):
            self.aligns = []
            self._load_splits()
        else:
            self.aligns = []
            self._split_query()

    def _load_splits(self):
        """ Loads splits from a previous run """
        LOG.info("Found existing split file. Loading that...")
        with open(self.split_f, 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                split_data = record.description.split("_", 1)[0]
                index, length = split_data.split("-len")
                self.aligns.append(Alignment(int(index), int(length)))

    def _split_query(self):
        """ 
        Splits query into fragments of around split_len

        TODO: Add a small algorithm to calculate the actual split_len based on the
        length of the contig to get it as close as possible to the supplied split_len

        Returns the total length of the sequences.
        """

        LOG.info("Splitting {} into pieces ~{}bp long.\nStoring as {}".format(self.query_f, self.split_len, self.split_f))

        split_number = 0 
        with open(self.query_f, 'r') as IN, open(self.split_f, 'w') as OUT:
            for seq_obj in SeqIO.parse(IN, 'fasta'):
                header = seq_obj.id
                seq = seq_obj.seq
                seq_len = len(seq)
            
                if seq_len < self.split_len:
                    continue


                # calculate split_len for this contig this should ensure all splits are close to 
                # one another in length
                leftovers = seq_len % self.split_len
                extra_per = int(leftovers / int(seq_len / self.split_len)) + 1

                # split the sequence into chunks and write each to the split_f 
                # prepending the split index to the header
                seq_indx = 0
                while seq_indx < seq_len:
                    # see if there are any leftovers remaining
                    if leftovers:
                        split_len = self.split_len + extra_per
                        leftovers -= extra_per
                    else:
                        split_len = self.split_len
     
                    # split out the seq
                    sp_seq = seq[seq_indx:seq_indx+split_len]

                    # calc the actual length of the seq
                    sp_seq_len = len(sp_seq)

                    # write a new header
                    sp_header = str(split_number) + "-len{}".format(sp_seq_len) + "_" + header
                    
                    # write the seq
                    sp_seq_obj = SeqRecord(sp_seq, id=sp_header, description='')
                    SeqIO.write(sp_seq_obj, OUT, 'fasta')
                    
                    # make an alignment object
                    self.aligns.append(Alignment(split_number, sp_seq_len))

                    split_number += 1
                    seq_indx += split_len

    def map_to_references(self, references, p_launcher=None, cpus=8):
        """ Maps the query to all references, optionally with the use of a parallelizer. """
        
        if p_launcher:
            args = {
                    'qry_split_f': self.split_f, 
                    'split_length': self.split_len,
                    'references': references, 
                    'cpus': cpus, 
                    'output_dir': self.out_dir
                    }

            self.parallel_id = p_launcher.run(args)
           
        else:
            self.mapped_files = map_to_references(self.split_f, self.split_len, references, cpus, self.out_dir)

    def get_ani_and_af(self, p_launcher=None):
        """ Returns two dicts, ani and af which are ani and af calculations keyed by reference name """

        # get the job from the parallelizer if necessary
        if self.parallel_id:
            self.mapped_files = p_launcher.get_results(self.parallel_id, wait=True)

        ani_per_ref = {}
        af_per_ref = {}
        if self.mapped_files:
            for ref, mapped_file in self.mapped_files.items():

                self.parse_SAM_file(mapped_file)

                # calculate gANI and AF
                total_length = 0
                aln_total_length = 0
                aln_aln_length = 0
                id_x_aln_length = 0
                num_aligned = 0
                for aln in self.aligns:
                    total_length += aln.length

                    # add up all alignments
                    if aln.aligned:
                        num_aligned += 1
                        aln_total_length += aln.length
                        aln_aln_length += aln.aln_len
                        id_x_aln_length += aln.aln_id * aln.aln_len

                try:
                    # the id is always divided by total length; not aligned length
                    gANI = id_x_aln_length / aln_total_length
                except ZeroDivisionError:
                    gANI = 0

                try:
                    # AF is scaled to length, not number of fragments
                    AF = aln_total_length / total_length
                except ZeroDivisionError:
                    AF = 0

                # store the result in the dict
                ani_per_ref[ref] = gANI
                af_per_ref[ref] = AF

                # write some output stats if in debug mode
                LOG.debug("For '{}' by '{}':".format(self.q_basename, ref))
                LOG.debug("\t{} of {} aligned.".format(num_aligned, len(self.aligns))) 
                LOG.debug("\t{} of {} length aligned.".format(aln_total_length, total_length))
                LOG.debug("\tANI={}\tAF={}\n".format(gANI, AF))

        return ani_per_ref, af_per_ref

    def parse_SAM_file(self, sam_f):
        """ Finds the hits in SAM file and sets the aligns to found as appropriate"""

        # reset all the aligns 
        for aln in self.aligns:
            aln.reset()

        with open(sam_f, 'r') as IN:
            for record in samparser.parse(IN, mapq=1, aligned_only=True):
                if record.mapped and record.perc_id >= self.min_id:
                    split_number, length = record.qname.split("_", 1)[0].split("-")
                    split_number = int(split_number)
                    if (record.length / self.aligns[split_number].length) >= self.min_cov:
                        self.aligns[split_number].report_found(record.perc_id, record.length)
  

class Alignment(object):
    """ Class used to record the splits of a query and if they align to a reference. """

    def __init__(self, index, length):
        self.index = index
        self.length = length

        self.aligned = False
        self.aln_len = 0
        self.aln_id = 0

    def __str__(self):
        if self.aligned:
            return "{}; len={}; aligned; aln_id={}; aln_len={}".format(self.index, self.length, self.aln_id, self.aln_len)
        else:
            return "{}; len={}; unaligned".format(self.index, self.length)

    def report_found(self, aln_id, aln_len):
        """ Sets the alignment to a found state"""
        self.aligned = True
        self.aln_id = aln_id
        self.aln_len = aln_len

    def reset(self):
        """ Resets the alignment to an unfound state """
        self.aligned = False
        self.aln_id = 0
        self.aln_len = 0


def map_to_references(qry_split_f, split_length, references, cpus, output_dir):
    """ Worker function to map a batch of files, returns a dict of SAM files from the mapping """

    q_basename = os.path.splitext(os.path.basename(qry_split_f))[0]
    q_basename = q_basename.rsplit(".split", 1)[0]      # second step to get the split suffix out

    sam_files = {}
    for ref in references:
        
        r_basename = os.path.splitext(os.path.basename(ref))[0]

        bbmap_f = output_dir + q_basename + "__--__" + r_basename + ".sam"

        # skip if this is done
        if os.path.isfile(bbmap_f):
            print("Found SAM file for '{}'. Skipping.".format(r_basename))
            sam_files[r_basename] = bbmap_f
            continue

        if split_length <= 500:
            prog = "bbmap.sh"
        else:
            prog = "mapPacBio.sh"

        cmd = "{} ref={ref} in={split_f} local=t ssao=f secondary=f nodisk overwrite=t sam=1.4 threads={cpus} out={out}".format(
                prog,
                ref=ref,
                split_f=qry_split_f,
                cpus=cpus,
                out=bbmap_f)

        print("Running:\n    {}".format(cmd))
        #code = subprocess.call(cmd.split(" "))         # safe way
    
        code = subprocess.call(cmd, shell=True)         # dangerous way
        code = 0

        if code:
            raise Exception("The bbmap command failed")
        else:
            sam_files[r_basename] = bbmap_f

    return sam_files

def main(args):

    # switch between all-by-all and reference based 
    if args.ref:
        LOG.debug("Using all in fastas as the reference group.")
        refs = args.ref
    else:
        LOG.debug("Using fastas supplied with -ref as the reference group.")
        refs = args.i
    
    # start a parallelizer if needed 
    if args.nodes > 0:
        p_launcher = parallelizer.Parallelizer(map_to_references, args.nodes, args.cpus, imports=["subprocess, os"], job_prefix="ani_af") 
    else:
        p_launcher = None
    
    # start all the mapping jobs
    LOG.info("Beginning mapping all by all...")
    mappers = []
    for qry in args.i:
        qm = QueryMapper(qry, split_len=args.split_len)
        qm.map_to_references(refs, p_launcher, args.cpus)
        mappers.append(qm)


    # process the results
    LOG.info("Parsing SAM files...")
    ani_and_af_matr = {}
    for qm in mappers:
        ani, af = qm.get_ani_and_af(p_launcher)
        ani_and_af_matr[qm.q_basename] = {}

        for ref, value in ani.items():
            ani_and_af_matr[qm.q_basename][ref] = (ani[ref], af[ref])

    qrys = sorted(ani_and_af_matr.keys())
    refs = sorted(ani_and_af_matr[qrys[0]].keys())

    # write the correlation table with tuples and also the classification table
    ani_out = args.prefix + ".ANI.tab"
    af_out = args.prefix + ".AF.tab"
    with open(ani_out, 'w') as ANI, open(af_out, 'w') as AF:
        ANI.write("\t".join(["query"] + refs) + "\n")
        AF.write("\t".join(["query"] + refs) + "\n")
        for q in qrys:
            ani_arr = [q]
            af_arr = [q]
            for r in refs:
                # add the val to be written to the quantified table
                ani_arr.append(str(ani_and_af_matr[q][r][0]))
                af_arr.append(str(ani_and_af_matr[q][r][1]))


            ANI.write("\t".join(ani_arr) + "\n")
            AF.write("\t".join(af_arr) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates ANI (average nucleotide identity) and AF (alignment fraction) between two genomes. This is done by breaking up each query genome into pieces and mapping those pieces back to the other genome in the pair. AF is the length of the fraction of pieces that map back. ANI is the average identity of the length that maps back. This script outputs two matrices where each cell is either ANI or AF. Optionally splits up the mapping across multiple nodes of an LSF cluster. This method is based on the paper 'Microbial species delineation using whole genome sequences' by Varghese et al http://nar.oxfordjournals.org/content/early/2015/07/06/nar.gkv657.full") 

    parser.add_argument("-i", help="bunch of fasta files to all-by-all compare", nargs="+", required=True)
    parser.add_argument("-ref", help="bunch of fasta files to use as the reference group (if you don't want all-by-all)", nargs="+")
    parser.add_argument("-split_len", help="approx length of fragment [%(default)s]", type=int, default=1000)
    parser.add_argument("-prefix", help="prefix for the tables [%(default)s]", default="genome_comparison")
    parser.add_argument("-nodes", help="the number of nodes to parallelize to, 0 for local only [%(default)s]", default=0, type=int)
    parser.add_argument("-cpus", help="the number of cpus to use for mapping [%(default)s]", default=4, type=int)
    args = parser.parse_args()

    main(args)
