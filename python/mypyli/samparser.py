
import re
import sys

def parse(sam_fh, aligned_only=False, mapq=1):
    for line in sam_fh:
        if line.startswith("@"):
            continue
        # store all the basic sam attributes
        aln_dict = {k: v for k, v in zip(['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'], line[:-1].split('\t')[:11])}

        aln_dict['pos'] = int(aln_dict['pos'])
        aln_dict['mapq'] = int(aln_dict['mapq'])
        aln_dict['flag'] = int(aln_dict['flag'])

        # store the optional attribs
        aln_dict.update({'op_fields': line[:-1].split("\t")[11:]})

        if aligned_only:
            if aln_dict['flag'] == 4:
                continue
        
        if aln_dict['mapq'] < mapq:
            continue

        yield SamRecord(aln_dict)

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

    def __init__(self, attrib):
        for k, v in attrib.items():
            setattr(self, k, v)

        if self.flag == 4:
            self.mapped = False
        else:
            self.mapped = True

        # hidden attributes to hold calculated values
        self._length = None
        self._matches = None
        self._perc_id = None

    @property
    def length(self):
        if self._length is None:
            self.length = len(self.seq)

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
            if "=" in self.cigar:       # sam 1.4
                m = re.findall("(\d+[MIDNSHP=X])", self.cigar)
                
                if m:
                    matches = 0
                    for group in m:
                        if group[-1] == "=":
                            matches += int(group[:-1])

                    self.matches = matches
                else:
                    raise ValueError("Cigar string {} could not be parsed. It may be malformed.".format(self.cigar))
            else:
                raise NotImplementedError("Right now, this requires SAM format 1.4")

        return self._matches
    @matches.setter
    def matches(self, value):
        self._matches = value
    
    @property
    def perc_id(self):
        if self._perc_id is None:
            self.perc_id = (self.matches / self.length)

        return self._perc_id
    @perc_id.setter
    def perc_id(self, value):
        self._perc_id = value

    

