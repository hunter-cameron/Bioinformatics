
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
                if piece.identity >= minid:
                    found += 1

        return found / total

    def get_fraction_in_order(self):
        for contig in self.contigs.values():
            contig.find_index_ordered()
        
        ordered = 0
        total = 0
        for piece in self.get_pieces():
            if piece.f_ordered:
                ordered += 1
            total += 1

        return ordered / total


class KingContig(object):
    def __str__(self):
        string = "Contig: {}. Length: {}\n    ".format(self.name, self.length)
        string += "\n    ".join([str(piece) for piece in self.get_pieces()])
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
        """ Yields sorted pieces """
        for piece in sorted(self.pieces, key=lambda k: self.pieces[k].index):
            yield self.pieces[piece]

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

    def find_index_ordered(self):
        """ 
        This method struggles to resolve repetitions, it is is mapped to
        a repetition rather than the piece, it may be placed as out of order. 
        """
        pieces = [piece for piece in self.get_pieces()]
        
        prev = None
        curr = None
        nfound = False
        next = -1   # set to -1 so on the first loop it will be set to 0

        while True:

            # skip to the next found piece
            while nfound is False:
                next += 1
                try:
                    nfound = pieces[next].found
                except IndexError:      # this is the signal to exit
                    nfound = True
                    next = None

            # in general, I want to check if the previous piece and the
            # next piece have the correct relation to the current piece
            # to begin with I'll assume they do not
            p = False
            n = False
        
            if curr is not None:
                #print(pieces[curr])
                #print(("curr", curr))
                #print("cnode, " + str(pieces[curr].cnode))
                #print("indx, " + str(CakeNode.contigs[pieces[curr].cnode.contig].index(pieces[curr].cnode)) + "/" + str(len(CakeNode.contigs[pieces[curr].cnode.contig]) - 1))
                ### check prev
                if prev is not None:
                    if pieces[curr].cnode.is_beside(pieces[prev].cnode):
                        p = True
                        #print("pmatch")
                    else:   # check for contig break in cake
                        if pieces[prev].cnode.is_terminal() and pieces[curr].cnode.is_terminal():
                            p = True
                            #print("pbreak")
                        else:   # otherwise it doesn't match
                            #print("pno")
                            pass

                        
                else:       # is there is no prev, that is ok
                    p = True
                    #print("pnone")

                ### check next
                if next is not None:
                    if pieces[curr].cnode.is_beside(pieces[next].cnode):
                        n = True
                        #print("nmatch")
                    else:
                        if pieces[curr].cnode.is_terminal() and pieces[next].cnode.is_terminal():
                            n = True
                            #print("nbreak")
                        else:   # otherwise it doesn't match
                            #print("nno")
                            pass
                else:       # if there is no next, that is ok
                    n = True
                    #print("nnone")

                pieces[curr].f_ordered = p and n

                #print(pieces[curr].f_ordered)
                #print()
            # increment/reset counters
            if next is None:
                break
            prev = curr
            curr = next
            nfound = False
                        
    def set_optimal_cnodes(self):
        """ 
        Sets the primary cnode for each KingPiece to be the one
        that maximizes the number of KingPieces that are properly ordered.
        
        This method is a rough draft for this task. There are probably multiple
        ways to optimize this and also multiple bugs that I have looked over. 

        Basically, this gets the job done now, but don't expect it to always
        work.
        """
        pieces = [piece for piece in self.get_pieces()]
      
        multiple_paths = []
        # I really only care about the nodes at which there are multiple
        for indx in range(len(pieces)):
            if len(pieces[indx].alignments) > 1:
                multiple_paths.append(indx)
            else:
                if pieces[indx].found:
                    pieces[indx].cnode = pieces[indx].alignments[0]

        # if there are no multiples, return here
        if not multiple_paths:
            return

        # get tuples of bounds for multiple regions
        ranges = []
        low_bound = 0
        for indx in range(1, len(multiple_paths)):
            if multiple_paths[indx] > 1 + multiple_paths[indx-1]:
                ranges.append((multiple_paths[low_bound], multiple_paths[indx-1]))
                low_bound = indx
        #print((low_bound, len(multiple_paths)))
        ranges.append((multiple_paths[low_bound], multiple_paths[-1]))

        #print(ranges)

        
        def place_alignment(d, tag, pieces):
            """ 
            Function to recurse over a dict and find all possible
            spots for a given alignment. 

            This is an in-method function because I wanted to use
            recursion but saw no point in creating a static method for
            the class.
            """
            # set the piece and alignment indices for the tag
            pti, ati = [int(val) for val in tag.split("_")]

            for k, v, in d.items():
                if isinstance(v, dict):
                    d[k] = place_alignment(v, tag, pieces)
                else:
                    # get piece and alignment index
                    pi, ai = [int(val) for val in k.split("_")]
                    # skip alignments from the same index
                    if pti == pi:
                        continue
                    if pieces[pti].alignments[ati].is_logically_cons(pieces[pi].alignments[ai]):
                        d[k] = {tag: 1}


            return d
        


        # try to resolve the best node for the multiples
        for irange in ranges:
            #print(("irange", irange))
            imin, imax = irange
            
            paths = {}
            for indx in range(imin, imax+1):
                for align_indx in range(len(pieces[indx].alignments)):
                    # recurse through all paths looking for all possible
                    # placements
                    tag = "{}_{}".format(indx, align_indx)
                    #print(("tag", tag))

                    paths =  place_alignment(paths, tag, pieces)
                    paths[tag] = 1

                    #[print((k, paths[k])) for k in sorted(paths)]
                    #print()
                    #print()


            def get_all_paths(d, prefix=[]):
                """ 
                Another function to recurse through a dict that didn't
                need to be a static class method.

                Returns a list of all paths.
                """
                paths = []

                for k, v in d.items():
                    if isinstance(v, dict):
                        paths += get_all_paths(v, prefix + [k])
                    else:
                        paths.append(prefix + [k])
                return paths

            paths = get_all_paths(paths)

            #[print(path) for path in paths]
            

            top_paths = [[]]
            for path in paths:
                if len(path) > len(top_paths[0]):
                    top_paths = [path]
                elif len(path) == len(top_paths[0]):
                    top_paths.append(path)

            #print(top_paths)

            # try to score the paths based on whether they match the nodes
            # before and after them. 
            scores = [0]*len(top_paths)
            
            # score previous node
            if imin > 0:
                prev = pieces[imin-1]
                if not prev.cnode:
                    if prev.found:
                        print("WARNING! cannot use node with multiple alignments as a reference.\nScore will be 0.", file=sys.stderr)

                else:
                    for indx in range(len(top_paths)):
                        #print(("path", top_paths[indx][0]))
                        pi, ai = [int(val) for val in top_paths[indx][0].split("_")]

                        if prev.cnode.is_beside(pieces[pi].alignments[ai]):
                            scores[indx] += 3
                        elif prev.cnode.is_logically_cons(pieces[pi].alignments[ai]):
                            scores[indx] += 1

            # score next node
            if imax < len(pieces) - 1:
                next = pieces[imax+1]
                if not next.cnode:
                    if next.found:
                        print("WARNING! cannot use node with multiple alignments as a reference.\nScore will be 0.", file=sys.stderr)
                else:
                    for indx in range(len(top_paths)):
                        #print(("path", top_paths[indx][0]))
                        pi, ai = [int(val) for val in top_paths[indx][-1].split("_")]
                        if next.cnode.is_beside(pieces[pi].alignments[ai]):
                            scores[indx] += 3
                        elif next.cnode.is_logically_cons(pieces[pi].alignments[ai]):
                            scores[indx] += 1


            # get the top score (if two are tied it will take the first)
            top = max(scores)
            top_indx = scores.index(top)
            #print(("top", scores, top_indx))
    
            # assgn cnodes 
            for tag in top_paths[top_indx]:
                pi, ai = [int(val) for val in tag.split("_")]
                pieces[pi].cnode = pieces[pi].alignments[ai]

            #sys.exit()
  
class KingPiece(object):
    def __str__(self):
        if self.found:
            if self.cnode:
                return "Piece {}. Location {} {}-{} | Found {} {}-{} {}% id, Ordered: {}".format(
                    str(self.index), self.contig.name, str(self.start), str(self.end),
                    self.cnode.contig, str(self.cnode.start), str(self.cnode.end),
                    str(self.identity * 100), str(self.f_ordered))
            else:
                return "Piece {}. Location {} {}-{} | Found ambig node, {}% id, Ordered: {}".format(
                    str(self.index), self.contig.name, str(self.start), str(self.end),
                    str(self.identity * 100), str(self.f_ordered))

        else:
            return "Piece {}. Start: {} End: {} | Not found".format(str(self.index), str(self.start), str(self.end))

        print(self.index)

    def __init__(self, contig, index, start=None, end=None):
        # identity metrics
        self.contig = contig
        self.index = index
        self.start = start
        if self.start:
            self.start += 1         # SAM uses 1 based coords
        self.end = end              # SAM uses 1 based coords, but no need to add 1 b/c 
                                    # SAM is inclusive and the top bound in python is not

        # attributes reserved for when (if) the node is found
        self.found = False
        self.cnode = None
        self.f_identity = None
        self.f_ordered = None

        # this holds all the alignments, the best should be determined by 
        # set optimal cnodes
        self.alignments = []
 
    def is_beside(self, piece):
        """ Determines if one piece is adjacent to another (ascending or descending) """
        return self.contig == piece.contig and self.index in [piece.index - 1, piece.index + 1]

    def register_found(self, cnode, perc_id):
        """ 
        Registers a piece as found, if piece has been found before, 
        this will update the repeats section
        """

        if self.found:
            self.alignments.append(cnode)
        else:
            self.alignments.append(cnode)
            # it should be fine to set the perc id for the first one only
            # to be ambiguous, they should be exactly the same in terms of
            # % id
            self.identity = perc_id
            self.found = True
    

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

class CakeNode(object):
    
    contigs = {}
    indexed = False
    
    @classmethod
    def add_node(cls, contig, start, end, piece):
        """ 
        Add the node if it doesn't exist.
        
        I want to add nodes like this rather than with the __init__ traditional
        way because I don't want to add duplicate nodes but I want to 
        return the node that the index was assigned to. __init__ would just
        return the new node and not care that the node wasn't actually added
        to the dict.
        
        Perhaps I can use add node to set self in the __init__?
        """
        node = cls(contig, start, end, piece)

        cls.contigs[node.contig] = cls.contigs.get(node.contig, [])

        # if node has been found before, add a new piece
        if node in cls.contigs[node.contig]:
            #print("Node already found {}".format(str(node.start)))
            indx = cls.contigs[node.contig].index(node)
            cls.contigs[node.contig][indx].pieces += node.pieces
            return cls.contigs[node.contig][indx]
        # otherwise add the node
        else:
            cls.contigs[node.contig].append(node)
            return node

    @classmethod
    def index_nodes(cls):
        """ Sort the nodes and set the prev and next attributes of each """
        for contig in cls.contigs:
            cls.contigs[contig].sort()
            for indx, node in enumerate(cls.contigs[contig]):
                if indx > 0:
                    node.prev = cls.contigs[contig][indx-1]
                if indx < len(cls.contigs[contig]) - 1:
                    node.next = cls.contigs[contig][indx+1]

            #[print(node) for node in cls.contigs[contig]]

        cls.indexed = True

    def __str__(self):
        return "\t".join([self.contig, self.start, self.end])

    def __init__(self, contig, start, end, piece):
        self.pieces = [piece]
        self.contig = contig
        self.start = start
        self.end = end

        self.prev = None
        self.next = None

        #self._add_node(self)


    def __str__(self):
        return "CakeNode - {}, {}-{}".format(self.contig, str(self.start), str(self.end))

    # comparison operator overloads
    def __eq__(self, other):
        return (self.contig == other.contig) and (self.start == other.start) and (self.end == other.end)

    def __ne__(self, other):
        return not self == other
    
    def __lt__(self, other):
        if self.contig == other.contig:
            return self.start < other.start
        else:
            return self.contig < other.contig

    def __le__(self, other):
        return not other < self

    def __gt__(self, other):
        if self.contig == other.contig:
            return self.start > other.start
        else:
            return self.contig > other.contig
    
    def __ge__(self, other):
        return not self < other

    def is_beside(self, other):
        beside = []
        if self.prev:
            beside.append(self.prev)
        if self.next:
            beside.append(self.next)

        #print((str(self), str(other), [str(item) for item in beside]))

        #print((other, beside))
        return other in beside

    def is_logically_cons(self, other):
        """
        This is a bit of a convenience method. It returns a boolean of whether
        one node "is logically consistent" with another in terms of order.

        In order to be logically consistent, the nodes must be beside one
        another or both nodes must be terminal on a contig (ie. there is a 
        contig break in the cake assembly that wasn't present in the king asm.)
        """

        if self.is_beside(other):
            return True
        else:
            return self.is_terminal() and other.is_terminal()

    def is_at_end(self):
        return self.contigs[self.contig][-1] == self

    def is_at_beginning(self):
        return self.contigs[self.contig][0] == self

    def is_terminal(self):
        return self.is_at_beginning() or self.is_at_end()

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
    """ 
    Runs bbmap and returns a path of the output sam.
    Prints all top alignments for ambiguously mapped reads.
    """
    cpus = 4
    if max_len <= 500:
        prog = "bbmap.sh"
    else:
        prog = "mapPacBio.sh"

    #cmd = "{} ref={cake_f} in={split_f} local=f ordered=t ssao=t nodisk overwrite=t sam=1.4 threads={cpus} out={out}".format(
    cmd = "{} ref={cake_f} in={split_f} local=f ordered=t ssao=f secondary=f nodisk overwrite=t sam=1.4 threads={cpus} out={out}".format(
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
                    raise ValueError("Regex didn't match SAM reference! {}".format(record.qname))
                
                piece = king.lookup_piece(contig, index)

                # make the CakeNode
                # NOTE: Check this!
                # end position is -1 to offset the length -- is this right??
                cnode = CakeNode.add_node(record.rname, record.pos, record.pos + record.length - 1, piece)

                # register the piece as found
                piece.register_found(cnode, record.perc_id)
      

def gene_pipeline(args):
    # this regex is pretty messy but the end matches a ( and then skips
    # until a second opening ( to capture the contig
    
    # this beginning anchored regex fails for some JGI names (they have some
    # sort of bug in their python code with the positions)
    #regex='[^ ]+\ [^_]+_(?P<index>\d+)\ [^(]+\([^(]+\((?P<contig>[^)]+).*'
    
    # This is a begin and end anchored regex that will hopefull be better
    # this takes advantage of the fact(?) that names seem to only have one
    # of brackets
    regex='[^ ]+\ [^_]+_(?P<index>\d+).*\((?P<contig>[^)]+)\) \[.*$'
    king = build_king(args.s, regex=regex)

    sam_f = run_bbmap(args.a, args.s, max_len=6000, out="king.sam")
    #sam_f = 'king.sam'
    find_pieces(sam_f, king, regex=regex)
    
    CakeNode.index_nodes()

    # use the name to print the results
    fasn = os.path.splitext(args.a)[0]

    #[print(node.start, node.end) for node in CakeNode.contigs['N515DRAFT_scaffold00013.13']]

    #[print(str(node)) for node in CakeNode.contigs['N515DRAFT_scaffold00001.1'][:10]]
    #king.contigs['N515DRAFT_scaffold00001.1'].set_optimal_cnodes()
    #king.contigs['N515DRAFT_scaffold00001.1'].find_index_ordered()
    #print(king.contigs['N515DRAFT_scaffold00001.1'])

    for contig in king.contigs.values():
        contig.set_optimal_cnodes()

    print("{}\t{}\t{}".format(fasn, str(king.get_fraction_found(minid=.95) * 100), str(king.get_fraction_in_order() * 100)))

    #print(str(king))

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
    
    #sam_f = run_bbmap(args.a, split_f, args.l)
    sam_f = "king.sam" 
    get_cake_lengths(sam_f)
   
    find_pieces(sam_f, king)
    print("Found {} percent of the spike in pieces present at > 90% identity.".format(str(king.get_fraction_found() * 100)))
    
