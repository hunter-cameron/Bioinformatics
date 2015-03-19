
import argparse
import sys
import subprocess
import multiprocessing
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from mypyli import samparser
import os

class King(object):
    def __str__(self):
        string = "King: {}\n\n  ".format(self.name)
        string += "\n\n  ".join([str(contig) for contig in self.contigs.values()])
        return string

    def __init__(self, name):
        self.name = name
        self.contigs = {}

    def add_contig(self, contig, length=None):
        """ Adds a contig to the king """
        contig_obj = KingContig(contig, length)
        if contig in self.contigs:
            raise ValueError("Contig {} already exists!".format(str(index)))
        else:
            self.contigs[contig] = contig_obj

        return contig_obj

    def lookup_contig(self, contig):
        return self.contigs[contig]

    def lookup_piece(self, contig, index):
        """ Utility function to go right to a piece """
        return self.lookup_contig(contig).lookup_piece(index)

    def get_pieces(self):
        """ Convience function to iterate over pieces """
        for contig in self.contigs.values():
            for piece in contig.get_pieces():
                yield piece

    def get_fraction_found(self, minid=.90):
        """ 
        This will underrepresent the percentage found because it will
        count the fragment at the end of each contig as missing unless it is
        90% present. This may not be so bad though, because it should never be
        continued on another contig or it would have been a single contig
        in the first place.
        """
        found = 0
        total = 0
        for piece in self.get_pieces():
            total += 1
            if piece.found:
                if piece.f_identity >= minid:
                    found += 1

        return found / total


class KingContig(object):
    def __str__(self):
        string = "Contig: {}. Length: {}\n    ".format(self.name, self.length)
        string += "\n    ".join([str(piece) for piece in self.pieces.values()])
        return string

    def __init__(self, name, length=None):
        self.name = name
        self.length = length

        self.pieces = {}

    def add_piece(self, index, start=None, end=None):
        piece = KingPiece(self, index, start, end)
        if index in self.pieces:
            raise ValueError("Index {} already exists!".format(str(index)))
        else:
            self.pieces[str(index)] = piece
        
        return piece

    def lookup_piece(self, index):
        return self.pieces[str(index)]

    def get_pieces(self):
        for piece in self.pieces.values():
            yield piece

    def reconstruct(self):
        """ 
        Converts the king contig to a list of mapped ranges,
        Unmapped pieces are labeled unmapped

        Esentially returns an assembly of the King contig in terms of Cake contigs.
        """
        mapped_ranges = []
        ordered_pieces = [self.pieces[indx] for indx in sorted(self.pieces, key=lambda k: int(key))]
        

        for indx, piece in enumerate(ordered_pieces):
            pass   

        return mapped_ranges


class KingPiece(object):
    def __str__(self):

        return "Piece {}. Start: {} End: {}".format(str(self.index), str(self.start), str(self.end))

    def __init__(self, contig, index, start=None, end=None):
        # identity metrics
        self.contig = contig
        self.index = index
        self.start = start
        if self.start:
            self.start += 1         # SAM uses 1 based coords
        self.end = end              # SAM uses 1 based coords, but no need to add 1 b/c 
                                    # SAM is inclusive and the top bound in python is not

        self.found = False
        self.f_contig = None
        self.f_start = None
        self.f_end = None
        self.f_identity = None

    def is_adjacent_to(self, piece):
        """ Determines if one piece is adjacent to another (ascending or descending) """
        return self.contig == piece.contig and self.index in [piece.index - 1, piece.index + 1]


class CakeContig(object):
    contigs = {}

    @classmethod
    def lookup_contig(cls, name):
        return cls.contigs[name]

    @classmethod
    def sort_all_pieces(cls):
        # sort all king_pieces
        [contig._sort_king_pieces() for contig in cls.contigs.values()]

    def __init__(self, name, length):
        self.name = name
        self.length = {}

        self.king_pieces = []

        self.contigs[name] = self

    def add_king_piece(self, piece):
        self.king_pieces.append(piece)
    

    def _sort_king_pieces(self):
        """ This should be called before trying to make the assembly """
        self.king_pieces = [piece for piece in sorted(self.king_pieces, key=lambda piece: piece.f_start)]

    def get_global_alignment(self):
        """ 
        Returns an alignment for the whole contig with gaps in appropriate places. 
        """
        alignment = []
        for indx, piece in enumerate(self.king_pieces):
            # if first piece, need to check if it aligns at beginning
            if not alignment:
                # if not at beginning, need to add a gap
                if piece.f_start != 1:
                    alignment.append((1, piece.f_start - 1))

                alignment.append(piece)

            else:
                # if start of this one isn't one after the end of the previous
                # we need a (negative?) gap
                if piece.f_start != alignment[-1].f_end + 1:
                    alignment.append((f_end + 1, piece.f_start -1))




def build_king(fasta_f, regex):
    """ 
    Builds the king structure 

    Needs a regex with two named capture groups, contig and index
    The contig group may be any string but the index group must be an integer
    and should correspond to the ordering of the sequences (1 comes before 2)

    In the future, there will be support for matching 'start' and 'end' and
    possibly even contig length.
    """
    king = King(name="king")
    with open(fasta_f) as IN:
        for seq_obj in SeqIO.parse(IN, 'fasta'):
            header = seq_obj.description
            seq = seq_obj.seq

            #print(header)
            m = re.match(regex, header)

            if m:
                contig = m.group('contig')
                index = int(m.group('index'))
                #print(("contig", index))
            else:
                raise ValueError("Regex didn't match!")

            # try to lookup the contig or create it
            try:
                contig_obj = king.lookup_contig(contig)
            except LookupError:
                contig_obj = king.add_contig(contig)

            contig_obj.add_piece(index)

    return king

def split_king_fasta(king_f, sp_len, split_f="split_king.fasta"):
    """ This function splits the fasta into pieces and makes the king object """
    king = King(name="king")
    with open(king_f, 'r') as IN, open(split_f, 'w') as OUT:
        for seq_obj in SeqIO.parse(IN, 'fasta'):
            header = seq_obj.id
            seq = seq_obj.seq
            seq_len = len(seq)

            contig = king.add_contig(header, length=seq_len)
            
            # split the sequence into chunks and write each to the split_f 
            # prepending the split index to the header
            for sp_indx, seq_indx in enumerate(range(0, len(seq), sp_len)):
                # the min portion of this is for the last slice when the last indx
                # may be longer than the sequence
                contig.add_piece(sp_indx, seq_indx, min([seq_indx+sp_len, seq_len-1]))
                
                sp_seq = seq[seq_indx:seq_indx+sp_len]
                sp_header = str(sp_indx) + "_" + header

                sp_seq_obj = SeqRecord(sp_seq, id=sp_header, description='')
                SeqIO.write(sp_seq_obj, OUT, 'fasta')

    return king, split_f

def run_bbmap(cake_f, split_f, max_len, out="split_king.sam"):
    """ Runs bbmap and returns a path of the output sam """
    cpus = 4
    if max_len <= 500:
        prog = "bbmap.sh"
    else:
        prog = "mapPacBio.sh"

    cmd = "{} ref={cake_f} in={split_f} local=f ordered=t nodisk sam=1.4 threads={cpus} out={out}".format(
            prog,
            cake_f=cake_f,
            split_f=split_f,
            cpus=cpus,
            out=out)

    print("Running:\n    {}".format(cmd), file=sys.stderr)
    #code = subprocess.call(cmd.split(" "))         # safe way
    
    code = subprocess.call(cmd, shell=True)         # dangerous way

    if code:
        raise Exception("The bbmap command failed")
    else:
        return out

def get_cake_lengths(cake_f):
    with open(cake_f, 'r') as IN:
        headers = samparser.parse_headers(IN)

    for contig, length in headers['seqs'].items():
        CakeContig(contig, length)

def find_pieces(sam_f, king, regex):
    """ Finds pieces of the king from the SAM file """
    with open(sam_f, 'r') as IN:
        for record in samparser.parse(IN, mapq=0):
            if record.mapped:
                m = re.match(regex, record.qname)
                if m:
                    index = int(m.group('index'))
                    contig = m.group('contig')
                else:
                    raise ValueError("Regex didn't match SAM reference!")
                
                piece = king.lookup_piece(contig, index)

                # register the piece as found
                piece.found = True
                piece.f_contig = record.rname
                piece.f_start = record.pos
                piece.f_end = record.length
                piece.f_identity = record.perc_id
      
                # add the king piece to the cake
                #cake = CakeContig.lookup_contig(record.rname)
                #cake.add_king_piece(piece)


def reconstruct_king(king):
    king_assembly = []
    for cake_contig in CakeContig.contigs:
        if self.king_pieces:
            pass
            

def gene_pipeline(args):
    # this regex is pretty complicated but the end matches a ( and then skips
    # until a second opening ( to capture the contig
    regex='[^ ]+\ [^_]+_(?P<index>\d+)\ [^(]+\([^(]+\((?P<contig>[^)]+).*'
    king = build_king(args.s, regex=regex)

    sam_f = run_bbmap(args.a, args.s, max_len=6000, out="king.sam")

    find_pieces(sam_f, king, regex=regex)
    
    # use the name to print the results
    fasn = os.path.splitext(args.a)[0]

    print("{}\t{}".format(fasn, str(king.get_fraction_found(minid=.95) * 100)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Looks for the information from a trusted genome assembly in a metagenome by chopping the trusted genome into pieces of length -l and searching for them in the metagenome.")
    parser.add_argument("-s", help="fasta from the spiked-in genome (king)", required=True)
    parser.add_argument("-a", help="fasta from the assembly (cake)", required=True)
    parser.add_argument("-l", help="length of the piece to chop king genome into", type=int, default=1000)
    parser.add_argument("-g", action='store_true', help="flag to indicate the -s fasta is a fasta of genes from JGI and should not be spilt")
    args = parser.parse_args()

    if args.g:
        gene_pipeline(args)
        sys.exit()

    king, split_f = split_king_fasta(args.s, args.l)
    
    sam_f = run_bbmap(args.a, split_f, args.l)
    
    get_cake_lengths(sam_f)
   
    find_pieces(sam_f, king)
    print("Found {} percent of the spike in pieces present at > 90% identity.".format(str(king.get_fraction_found() * 100)))
    
