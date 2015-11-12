
import re
import sys

def parse(sam_fh, aligned_only=False, mapq=1, local=False):
    for line in sam_fh:
        if line.startswith("@"):
            continue
        # store all the basic sam attributes
        aln_dict = {k: v for k, v in zip(['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'], line[:-1].split('\t')[:11])}

        aln_dict['pos'] = int(aln_dict['pos'])
        aln_dict['mapq'] = int(aln_dict['mapq'])
        aln_dict['flag'] = int(aln_dict['flag'])

        # store the optional attribs
        attrib_dict = {}
        for attrib in line[:-1].split("\t")[11:]:
            name, type, value = attrib.split(":")

            if type == "f":
                attrib_dict[name] = float(value)
            elif type == "i":
                attrib_dict[name] = int(value)
            else:
                attrib_dict[name] = value

        aln_dict["attributes"] = attrib_dict

        if aligned_only:
            if aln_dict['flag'] in [4, 516]:
                continue
        
        if aln_dict['mapq'] < mapq:
            continue

        yield SamRecord(aln_dict, local)

def parse_headers(sam_fh):
    headers = {'seqs': {}}
    for line in sam_fh:
        if line.startswith("@"):

            if line.startswith("@SQ"):
                sn, ln = line.split("\t")[1:]
                name = sn[3:]
                length = ln[3:]
                headers['seqs'][name] = int(length)
        else:
            return headers


class SamRecord(object):
    """ 
    A single read in a SAM file
    
    Attributes are added as needed to completely minimize computation (nothing
    unnecessary computed, no computations performed twice)
    """

    def __init__(self, attrib, local=False):
        for k, v in attrib.items():
            setattr(self, k, v)
        
        self.local = local


        if self.flag == 4:
            self.mapped = False
        else:
            self.mapped = True

        # hidden attributes to hold calculated values
        self._length = None
        self._matches = None
        self._perc_id = None
        self._mismatches = None
        self._soft_clipped = None

    @property
    def length(self):
        if self._length is None:
            self._parse_cigar()

        if self.local:
            return self._length - self.soft_clipped
        else:
            return self._length
    @length.setter
    def length(self, value):
        self._length = value

    @property
    def matches(self):
        """ 
        Right now this just looks for matches in the Cigar string of
        SAM cigar format 1.4
        """
        if self._matches is None:
            self._parse_cigar()

        return self._matches
    @matches.setter
    def matches(self, value):
        self._matches = value
    
    @property
    def mismatches(self):
        if self._mismatches is None:
            self._parse_cigar()

        return self._mismatches
    @mismatches.setter
    def mismatches(self, value):
        self._mismatches = value

    @property
    def soft_clipped(self):
        if self._soft_clipped is None:
            self._parse_cigar()

        return self._soft_clipped
    @soft_clipped.setter
    def soft_clipped(self, value):
        self._soft_clipped = value

    @property
    def perc_id(self):
        if self._perc_id is None:
            self.perc_id = (self.matches / self.length)

        return self._perc_id
    @perc_id.setter
    def perc_id(self, value):
        self._perc_id = value

    def _parse_cigar(self):
        cigar_stats = {k: 0 for k in ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X']}
        length = 0
        m = re.findall("(\d+[MIDNSHP=X])", self.cigar)
        if m:
                
            for group in m:
                cigar_stats[group[-1]] = cigar_stats.get(group[-1], 0) + int(group[:-1])
                length += int(group[:-1])
                
            self.length = length
            self.soft_clipped = cigar_stats['S']
                    
        else:
            raise ValueError("Cigar string {} could not be parsed. It may be malformed.".format(self.cigar))
        
        # check for format 1.4
        if "=" in self.cigar:
            self.matches = cigar_stats['=']
            self.mismatches = cigar_stats['X']
        else:   # old format
            # we can get something like mismatches from the edit distance
            try:
                edit_dist = self.attributes["NM"]
            except KeyError:
                raise ValueError("Entry doesn't have the edit distance field 'NM' and is in the old style SAM format. It is impossible to calculate mismatches.")

            # mismatches is the edit distance - insertions and deletions
            self.mismatches = edit_dist - (cigar_stats["I"] + cigar_stats["D"])
            self.matches = self.length - edit_dist
