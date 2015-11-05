
import sys
import argparse
import logging
import os
from Bio import SeqIO 
from Bio.SeqRecord import SeqRecord
import subprocess

from mypyli import samparser

logging.basicConfig()
LOG = logging.getLogger()
LOG.setLevel("DEBUG")

class QueryMapper(object):

    def __init__(self, query_f, out_dir=None, split_len=1000, min_id=.90, min_len=500):
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
        for record in self.parse_SAM_file(bbmap_f):
            if record.mapped and record.perc_id >= self.min_id:
                qry_region = GenomeRegion.from_header(record.qname)
    
                if record.length >= self.min_len:
                    self.regions[qry_region.index].report_hit(r_basename, record.rname, record.pos, record.length, record.perc_id)


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

    @staticmethod
    def parse_SAM_file(sam_f):
        """ Yields SamRecords that meet some minimum criteria"""
        with open(sam_f, 'r') as IN:
            for record in samparser.parse(IN, mapq=1, aligned_only=True):
                yield record


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

    def __str__(self):
        return self.to_header()

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
        
        if genome in self.alignments:
            raise ValueError("There has already been an alignment found in genome {} for {}".format(genome, str(self)))

        ref_region = GenomeRegion(genome, contig, start - 1, (start - 1) + length, self.index) 
        
        self.alignments[genome] = ref_region

    def has_hits_from(self, genomes):
        """ Checks if alignment has hits from a list of genomes. Returns True is so, else False """

        for genome in genomes:
            if genome not in self.alignments:
                return False

        return True


class Alignment(object):

    def __init__(self, index, header, start, length):
        self.index = index
        self.header = header
        self.start = start
        self.length = length

        self.aligned = False
        self.aln_id = 0
        self.aln_contig = ""
        self.aln_len = 0
        self.aln_start = 0

    def __str__(self):
        if self.aligned:
            return "{}; len={}; aligned; aln_id={}; aln_len={}".format(self.index, self.length, self.aln_id, self.aln_len)
        else:
            return "{}; len={}; unaligned".format(self.index, self.length)

    @property
    def aln_end(self):
        # aln_start is base 0, len is a distance, aln_end is closed so real end = (end -1)
        # so there is no need to +- 1 anywhere here
        return self.aln_start + self.aln_len

    def report_found(self, aln_id, aln_contig, aln_start, aln_len):
        self.aligned = True
        self.aln_id = aln_id
        self.aln_contig = aln_contig
        self.aln_start = aln_start - 1
        self.aln_len = aln_len

        #print(str(self))

    def reset(self):
        self.aligned = False
        self.aln_id = 0
        self.aln_len = 0


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

    for region in mapped_all[:10]:
        print(region)
    sys.exit()

    #
    ## build a dict with all the seqs to write (allows a single read/write phase)
    #

    LOG.info("Collating sequences to prepare for writing...")

    to_write = {}
    if mapped_all:
        
        # process the qry separately and first
        with open(qm.split_f, 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                params = QueryMapper.header_to_param_dict(record.description)

                if params["index"] in mapped_all:
                    # add the new subseq to the dict of seqs to write
                    try:
                        to_write[params["index"]].append(record)
                    except KeyError:
                        to_write[params["index"]] = [record]
 

        for ref in args.i:
            # skip the query here
            if ref == qry:
               continue 

            ref_basename = os.path.splitext(os.path.basename(ref))[0]
            contigs = set([conserved_regions[ref][indx].aln_contig for indx in mapped_all])
            with open(ref, 'r') as IN:
                for record in SeqIO.parse(IN, "fasta"):

                    # this is a check to avoid looping through all alns for contigs that have no aln
                    if record.description in contigs:

                        # get each region that corresponds to this contig
                        for indx in mapped_all:
                            alignment = conserved_regions[ref][indx]
                            if alignment.aln_contig == record.description:

                                # make a new seqobj with the alignment params
                                seq = record.seq[alignment.aln_start:alignment.aln_end]
                                header = QueryMapper.param_dict_to_header({"contig": record.description, "start": alignment.aln_start, "end": alignment.aln_end, "genome": ref_basename, "index": indx})
                                seqobj = SeqRecord(seq, id=header, description="")

                                # add the new subseq to the dict of seqs to write
                                try:
                                    to_write[indx].append(seqobj)
                                except KeyError:
                                    to_write[indx] = [seqobj]
    else:
        LOG.warning("No alignments found in all references. Nothing further to be done.") 

    #
    ## write all the seqs
    #

    LOG.info("Writing {} regions...".format(len(mapped_all)))

    for indx in to_write:
        with open("{prefix}_indx{indx}.fna".format(prefix=args.prefix, indx=str(indx)), 'w') as OUT:
            try:
                SeqIO.write(to_write[indx], OUT, "fasta")
            except TypeError:
                print(to_write[indx])
                raise





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculates ANI (average nucleotide identity) and AF (alignment fraction) between two genomes. This is done by breaking up each query genome into pieces and mapping those pieces back to the other genome in the pair. AF is the fraction of pieces that map back. ANI is the average identity of the pieces that map back. This script outputs a matrix where each cell is a tuple (ANI, AF). Currently, ANI is a percent and AF is a decimal but I should change them to make it a consistent measure before I put this in the scripts for others directory.")

    parser.add_argument("-i", help="bunch of fasta files to all-by-all compare", nargs="+", required=True)
    parser.add_argument("-qry", help="the file to be broken up and used as the query (if not the first one in '-i'")
    parser.add_argument("-split_len", help="target length to split query into", type=int, default=500)
    parser.add_argument("-min_id", help="minumum id for an alignment", type=float, default=.90)
    parser.add_argument("-min_len", help="minumum length for an alignment", type=int, default=500)
    parser.add_argument("-prefix", help="prefix for new fastas", default="quantified_table.tab")
    args = parser.parse_args()

    main(args)
