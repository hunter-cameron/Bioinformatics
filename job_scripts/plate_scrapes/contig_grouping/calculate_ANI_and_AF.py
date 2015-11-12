
import sys
import argparse
import logging
import os
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
import subprocess

from mypyli import samparser


LOG = logging.getLogger()


class QueryMapper(object):

    def __init__(self, query_f, out_dir=None, split_len=1000, min_id=70, min_cov=.70):
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

    def process(self, ref_f):

        # either split the query or set the list from the last time to None
        if self.aligns is None:
            self.aligns = []
            self.split_query()
        else:
            for aln in self.aligns:
                aln.reset()

        r_basename = os.path.splitext(os.path.basename(ref_f))[0]

        bbmap_f = self.out_dir + self.q_basename + "__--__" + r_basename + ".sam"

        # skipper if this is done
        if not os.path.isfile(bbmap_f):
            self.run_bbmap(ref_f, bbmap_f)

        self.parse_SAM_file(bbmap_f)

        # calculate gANI and AF
        total_length = 0
        aln_total_length = 0
        aln_aln_length = 0
        id_x_aln_length = 0
        for aln in self.aligns:
            total_length += aln.length

            if aln.aligned:
                aln_total_length += aln.length
                aln_aln_length += aln.aln_len
                id_x_aln_length += aln.aln_id * aln.aln_len

        try:
            gANI = id_x_aln_length / aln_total_length
        except ZeroDivisionError:
            gANI = 0

        try:
            AF = aln_total_length / total_length
        except ZeroDivisionError:
            AF = 0

        return gANI, AF

    def split_query(self):
        """ 
        Splits query into fragments of around split_len

        TODO: Add a small algorithm to calculate the actual split_len based on the
        length of the contig to get it as close as possible to the supplied split_len

        Returns the total length of the sequences.
        """

        LOG.info("Splitting {} into pieces ~{}\nStoring as {}".format(self.query_f, self.split_len, self.split_f))

        split_number = 0 
        with open(self.query_f, 'r') as IN, open(self.split_f, 'w') as OUT:
            for seq_obj in SeqIO.parse(IN, 'fasta'):
                header = seq_obj.id
                seq = seq_obj.seq
                seq_len = len(seq)
            
                if seq_len < self.split_len:
                    continue

                # calculate split_len for this contig this should ensure all contigs are close to 
                # one another
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

        cmd = "{} ref={ref} in={split_f} local=t ssao=f secondary=f nodisk overwrite=t sam=1.4 threads={cpus} out={out}".format(
                prog,
                ref=reference,
                split_f=self.split_f,
                cpus=cpus,
                out=out)

        print("Running:\n    {}".format(cmd), file=sys.stderr)
        #code = subprocess.call(cmd.split(" "))         # safe way
    
        code = subprocess.call(cmd, shell=True)         # dangerous way
        code = 0

        if code:
            raise Exception("The bbmap command failed")
        else:
            return out

    def parse_SAM_file(self, sam_f):
        """ Generates hits to indices """
        with open(sam_f, 'r') as IN:
            for record in samparser.parse(IN, mapq=1, aligned_only=True):
                if record.mapped and (record.perc_id * 100) >= self.min_id:
                    split_number, length = record.qname.split("_", 1)[0].split("-")
                    split_number = int(split_number)
                    if (record.length / self.aligns[split_number].length) >= self.min_cov:
                        self.aligns[split_number].report_found(record.perc_id * 100, record.length)
  


class Alignment(object):

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
        self.aligned = True
        self.aln_id = aln_id
        self.aln_len = aln_len

    def reset(self):
        self.aligned = False
        self.aln_id = 0
        self.aln_len = 0


def main(args):

    # switch between all-by-all and reference based 
    if args.ref:
        LOG.debug("Using all in fastas as the reference group.")
        refs = args.ref
    else:
        LOG.debug("Using fastas supplied with -ref as the reference group.")
        refs = args.i

    table = {}
    for qry in args.i:
        q_basename = os.path.splitext(os.path.basename(qry))[0]
        table[q_basename] = {}

        qm = QueryMapper(qry)
        for ref in refs:
            print("Processing {} by {}...".format(qry, ref))
            gANI, AF = qm.process(ref)
            table[q_basename][os.path.splitext(os.path.basename(ref))[0]] = (gANI, AF)

    qrys = sorted(table.keys())
    refs = sorted(table[qrys[0]].keys())

    # write the correlation table with tuples and also the classification table
    with open(args.out, 'w') as QUAN:
        print("\t".join(["query"] + refs), file=QUAN)
        for q in qrys:
            quant_arr = [q]
            clas_arr = [q]
            for r in refs:
                # add the val to be written to the quantified table
                quant_arr.append(str(table[q][r]))
                
                # do a quick classification 
                gANI, AF = table[q][r]
                code = classify_relationship(gANI, AF)
                clas_arr.append(str(code))
                

            print("\t".join(quant_arr), file=QUAN)

def classify_relationship(gANI, AF):
    """ Returns a classification index based on the following coding scheme:

    1 = A is part/all of B
    2 = A is part/all of B plus extras (either A is contaminated or B is incomplete)
    3 = A has few regions that are conserved in B

    4 = A and B are likely sister species/of the same genus (there are plenty of point mutations between the two but most of the content is still conserved)
    5 = A and B are somehow related, though it is hard to say how
    6 = A has a few regions that are conserved in B

    7 = A and B share many regions but there is little base-by-base conservation. This should be rare (if not impossible based on mapping params.
    8 = A and B have some regions conserved at low identity
    9 = A and B are not related

    """
    hi_ani = .80
    mid_ani = .50
    lo_ani = .30

    hi_af = .80
    mid_af = .50
    lo_af = .30

    if gANI >= hi_ani:
        if AF >= hi_af:
            return 1
        elif AF >= mid_af:
            return 2
        else:
            return 3

    elif gANI >= mid_ani:
        if AF >= hi_af:
            return 4
        elif AF >= mid_af:
            return 5
        else:
            return 6

    else:
        if AF >= hi_af:
            return 7
        elif AF >= mid_af:
            return 8
        else:
            return 9



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates ANI (average nucleotide identity) and AF (alignment fraction) between two genomes. This is done by breaking up each query genome into pieces and mapping those pieces back to the other genome in the pair. AF is the fraction of pieces that map back. ANI is the average identity of the pieces that map back. This script outputs a matrix where each cell is a tuple (ANI, AF). Currently, ANI is a percent and AF is a decimal but I should change them to make it a consistent measure before I put this in the scripts for others directory.")

    parser.add_argument("-i", help="bunch of fasta files to all-by-all compare", nargs="+", required=True)
    parser.add_argument("-ref", help="bunch of fasta files to use as the reference group (if you don't want all-by-all)", nargs="+")
    parser.add_argument("-out", help="filepath for the quantified table", default="quantified_table.tab")
    args = parser.parse_args()

    main(args)
