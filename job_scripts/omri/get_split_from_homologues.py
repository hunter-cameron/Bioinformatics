
import argparse
from Bio import AlignIO, SeqIO
import random

def get_primers_from_amplicons(amplicon_f):

    primers = {}
    with open(amplicon_f, 'r') as IN:
        for record in SeqIO.parse(IN, "fasta"):
            _, cluster, region = record.description.split(" ")[:3]

            seq_start = int(region.split("..")[0])
            seq_end = int(region.split("..")[1])

            # now, we need to find where the primer starts and ends
            first_N = None
            last_N = None
            for indx, nt in enumerate(record.seq):
                if first_N is None:
                    if nt == "N":
                        first_N = indx
                else:
                    if nt != "N":
                        last_N = indx - 1

            # make some primer objects
            fwd_primer = Primer(cluster, seq_start, seq_start + first_N)
            rev_primer = Primer(cluster, last_N, seq_end)

            # see if the primer has already been found, add it if not
            try:
                already_found = primers[cluster]
            except KeyError:
                primers[cluster] = [fwd_primer, rev_primer]
                continue

            if fwd_primer not in already_found:
                primers[cluster].append(fwd_primer)

            if rev_primer not in already_found:
                primers[cluster].append(rev_primer)
    return primers

class Primer(object):
    def __init__(self, cluster, start, end):
        self.cluster = cluster
        self.start = start
        self.end = end

        # to be added later
        self.seq = None
        self.contig = None
        self.cstart = None
        self.cend = None

    def __eq__(self, other):
        return self.cluster == other.cluster and self.start == other.start and self.end == other.end

    def __cmp__(self, other):
        if self.contig == other.contig:
            if self.cstart == other.cstart and self.cend == other.cend:
                return 0
            else:
                if self.cstart < other.cstart:
                    return -1
                else:
                    return 1
        else:
            first = sorted([self.contig, other.contig])
            if self.contig == first:
                return -1
            else:
                return 1

    def __lt__(self, other):

        if self.contig == other.contig:
            if self.cstart == other.cstart and self.cend == other.cend:
                return False
            else:
                if self.cstart < other.cstart:
                    return True
                else:
                    return False
        else:
            first = sorted([self.contig, other.contig])
            if self.contig == first:
                return True
            else:
                return False


    def generate_header(self, index):
        header = ">genome=13|index={indx}|contig={self.contig}|slice=[{self.cstart}:{self.cend}]".format(
                indx=str(index),
                self=self,
                )

        return header


def write_seqs_from_primers(primers, aligned_dir, out_f):
    for cluster in primers:
        print("Processing {}".format(cluster))
        with open(aligned_dir + "/" + cluster + ".faln", 'r') as IN:
            msa = AlignIO.read(IN, "fasta")
            
            # we only want to get the seq from genome 13
            target_seq = None
            for seq in msa:
                if seq.description.startswith("__Methylobacterium_13"):
                    target_seq = seq
                else:
                    continue

            
            # extract information from the header
            elems = target_seq.description.split("^")
            print(seq.description)

            # get start and end
            location = elems[0].split("__")[-1]     # remove first part of header
            location = location.split("_")[0]       # remove strand
            try:
                seq_start = int(location.split("-")[0]) - 1
                seq_end = int(location.split("-")[0])
            except ValueError:
                continue

            # get contig
            if "Contig" in elems[-1]:
                contig = elems[-1].split("_Contig_")[-1]
            else:
                continue

            gapless_seq = str(target_seq.seq).replace("-", "")

            for primer in primers[cluster]:
                # we need to count gaps before this position to get the real start
                gaps = target_seq.seq[:primer.start].count("_")

                real_start = primer.start - gaps

                primer_len = primer.end - primer.start
                # find out how much extra we need to get
                if primer_len >= 50:
                    primer_len = 50
                    front = 0
                    end = 0
                else:
                    extra = 50 - primer_len

                    # get a 50bp fragment that is within bounds
                    while True:
                        print("infinite loop?")
                        front = random.randint(0, extra)
                        end = 50 - front
                        if real_start - front >= 0 and real_start + primer_len + end <= len(gapless_seq):
                            break

                    
                primer.cstart = seq_start + (real_start - front)
                primer.cend = primer.cstart + primer_len + end
                primer.contig = contig
                primer.seq = gapless_seq[real_start - front:real_start + primer_len + end]

    with open(out_f, 'w') as OUT:

        # filter for good primers only
        prims = []
        for prim_list in primers.values():
            for prim in prim_list:
                if prim.contig is not None:
                    prims.append(prim)

        for indx, primer in enumerate(sorted(prims)):
            OUT.write(primer.generate_header(indx) + "\n" + primer.seq + "\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="this script is solely to write the primer regions from GH potential amplicons to use for testing. This seems like a futile task because often the headers are cut off so each sequence has no identifying information")
    parser.add_argument("-amplicons", help="potential amplicons from the GH run", required=True)
    parser.add_argument("-aligned_dir", help="the directory of aligned clusters", required=True)
    parser.add_argument("-out", help="the output FASTA", default="13.split.fasta")
    args = parser.parse_args()

    primers = get_primers_from_amplicons(args.amplicons)

    write_seqs_from_primers(primers, args.aligned_dir, args.out)
