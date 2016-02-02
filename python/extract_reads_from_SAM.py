
import argparse
import sys
import os
import logging

logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel("INFO")

class Alignment(object):
    """ An abbreviated SAM parser that stores only the fields necessary for the purposes here """

    def __init__(self, name, reference, seq, qual):
        self.name = name
        self.reference = reference
        self.seq = seq
        self.qual = qual

    def __str__(self):
        return self.name + "\t" + self.reference + "\t" +  self.seq

    @classmethod
    def from_sam_line(cls, sam_line):
        """ 
        Creates an Alignment Object from a SAM file line 
        
        We assume that mapped reads have reference field set and unmapped reads will have it set to '*' though this may not always be the case.
        """
        elems = sam_line.rstrip().split("\t")

        name = elems[0]
        reference = elems[2]
        seq = elems[9]
        qual = elems[10]


        if reference == "*":
            reference = None

        return cls(name, reference, seq, qual)

class FastqWriter(object):
    """ Manages writing to a fastq file in a buffered way """

    writers = []

    def __init__(self, fastq_f, mode='w', buffer_lines=4000):
        self.fastq_f = fastq_f
        self.mode = mode
        self.buffer_lines = buffer_lines

        self.lines_to_write = []

        # register this instance
        FastqWriter.writers.append(self)

    @classmethod
    def finish_writing_all(cls):
        """ Writes anything left in the buffer for each Fastq Writer """
        for writer in cls.writers:
            writer.write_to_fastq()

    def add_alignment(self, alignment):
        """ Converts an alignment to its fastq representation and adds the lines to the array """
        self.lines_to_write.append("@" + alignment.name)
        self.lines_to_write.append(alignment.seq)
        self.lines_to_write.append("+")
        self.lines_to_write.append(alignment.qual)

        # check if it is time to write
        if len(self.lines_to_write) >= self.buffer_lines:
            self.write_to_fastq()

    def write_to_fastq(self):
        """ Writes the lines to the fastq """

        LOG.debug("Writting to fastq '{}'".format(self.fastq_f))
        
        # don't do anything if there are no lines to write
        if not self.lines_to_write:
            return

        with open(self.fastq_f, self.mode) as OUT:
            OUT.write("\n".join(self.lines_to_write) + "\n")

        # reset the lines
        self.lines_to_write = []

        # switch to append mode after the first write
        if self.mode == 'w':
            self.mode = 'a'


def yield_alignments(sam_fh):
    """ 
    Yields all the alignments in a SAM file
    """
    for line in sam_fh:

        # skip all header lines 
        if line.startswith("@"):
            continue
        
        alignment = Alignment.from_sam_line(line)
        yield alignment


def parse_sam_file(sam_f, contigs_to_writers):
    """
    The general plan here is to write each pair of reads, where either maps to a contig in the list of contigs, only once.

    However, it is possible for either or both of the reads to have multiple alignments. 
    If any of the alignments are to the contig, we want to write that pair.

    My strategy is to process the reads one at a time and store the R1 and use a flag system to know whether or not to write the pair when processing the R2 read.

    """

    LOG.info("Parsing sam file '{}'".format(sam_f))
    with open(sam_f, 'r') as IN:

        aln_iter = yield_alignments(IN)

        # flags
        read1 = None            # store R1 reads
        read2 = None            # store R2 reads
        write = False           # reference was on list, add the pair to be written

        while True:
            try:
                aln = next(aln_iter)
            except StopIteration:
                # check if there is a R1 waiting for an R2
                if read1 and not read2:
                    print(read1)
                    print(read2)
                    raise ValueError("Reached the end of the file and haven't found the R2 for the stored read.")
                break
            
            # check if the current aln is another alignment of read2
            if read2:
                if aln.name == read2.name:
                    # read pair has already been written, skip all remaining alignments
                    if write:
                        continue

                # this is a new R1 -- reset
                else:
                    read1 = None
                    read2 = None
                    write = False


            # check if the alignment is on the list
            if aln.reference in contigs_to_writers:
                write = True

            # store this aln if there isn't already one stored (this is R1)
            if read1 is None:
                read1 = aln
                continue


            # check if alignment is the same as the previous
            # The names are not unique because pairs were renamed to have the same name
            # Just found a case where name, seq combo wasn't unique
            # Now, I'll include quality too...
            if aln.name == read1.name and aln.seq == read1.seq and aln.qual == read1.qual:
                continue
            else:
                # check if the aln is the proper pair of the previous aln
                if check_proper_pair(aln, read1):
                    read2 = aln

                    # write the alns and set to skip any further alignments
                    if write:
                        LOG.debug("Adding read pair to write")

                        # map to all contigs -- for most pairs this will be one contig
                        for contig in set((read1.reference, aln.reference)):
                            # skip unmapped contigs
                            if contig is None:
                                continue

                            try:
                                writers = contigs_to_writers[contig]
                            except KeyError:
                                # one of the contigs was not binned
                                continue

                            # add pair to each writer (more than 1 if contig was placed in 2 bins)
                            for writer in writers:
                                writer.add_alignment(read1)
                                writer.add_alignment(aln)

                else:
                    # if we are supposed to write but this isn't a proper pair something bad has happened
                    if write:
                        print(read1)
                        print(aln)
                        raise ValueError("Discovered a read that mapped to one of the specified references but did not have a pair. Either the SAM file is not sorted by name, the reads were not PE, or the SAM file is corrupted.")

def check_proper_pair(aln1, aln2):
    """ 
    Returns true if two alignments are a pair 
    
    This function may change for different read sets...
    """

    aln1_base = aln1.name.split(" ")
    aln2_base = aln2.name.split(" ")

    if aln1_base == aln2_base:
        return True
    else:
        return False

def parse_contigs_files(contigs_files, output_dir):
    """ Parses a list of contigs files and returns a dict linking contigs to a list of fastq writer objects """

    LOG.info("Parsing contigs files...")

    contigs_to_writers = {}

    # loop through all contigs files
    for contigs_file in contigs_files:
        basename = os.path.splitext(os.path.basename(contigs_file))[0]
        fastq_file = output_dir + "/" + basename + ".fastq"

        fastq_writer = FastqWriter(fastq_file)

        # loop through all contigs in the contigs file
        with open(contigs_file, 'r') as IN:
            for line in IN:
                contig = line.rstrip()
                
                try:
                    contigs_to_writers[contig].append(fastq_writer)
                except KeyError:
                    contigs_to_writers[contig] = [fastq_writer]

    return contigs_to_writers

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Utility used to extract all reads (and their mates) that map to a set of contigs. Input should be sorted by read name.")
    parser.add_argument("-samfiles", help="one or more samfiles", nargs="+", required=True)
    parser.add_argument("-contigs", help="file(s) with names of contigs to extract reads for reads will be extracted to a file using the basename of the contigs files", nargs="+", required=True)
    parser.add_argument("-out_dir", help="an output directory to put the extracted fastq files in [%(default)s]", default=os.getcwd())

    args = parser.parse_args()

    # unlikely race condition here..
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    # map contigs to fastq writers
    contigs_to_writers = parse_contigs_files(args.contigs, args.out_dir)

    # parse each sam file
    for sam_file in args.samfiles:
        parse_sam_file(sam_file, contigs_to_writers)

    # dump all the rest of the writers
    FastqWriter.finish_writing_all()
