import regex

def parse(blast_fh, outfmt):
    """ Function to provide similar use to biopython's SeqIO """
    if str(outfmt) == "6":
        return BlastParser.parse_outfmt_6(blast_fh=blast_fh)
    



class BlastParser(object):
    """ 
    Just a shell of a BLAST parser so I don't have to rewrite this. 
    Wow, this is a really, really terrible class...
    """
    
    def __init__(self, blast_f, outfmt):
        with open(blast_f, 'r') as IN:
            self.parse(IN, outfmt)
    
    @classmethod
    def parse_outfmt_6(cls, blast_fh):
        """ 
        Parses standard outfmt 6 BLAST results
        Returns: iterate Blastrecord objects
        Excepts: AssertionError 
        """

        for line in blast_fh:
            elems = line[:-1].split("\t")

            # attempt to make sure it is in standard format
            assert len(elems) == 12
            
            param_dict = {k: v for k, v in 
                    zip(["query", "subject", "perc_id", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"], elems)}

            # assign correct types to the data
            for param in ['length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send']:
                param_dict[param] = int(param_dict[param])

            for param in ['perc_id', 'evalue', 'bitscore']:
                param_dict[param] = float(param_dict[param])

            yield BlastRecord(param_dict)




class BlastRecord(object):
    """ An individual hit? Or a query? Right now just a single hit"""

    def __init__(self, param_dict):
        """ Needs param & error checking """
        for k, v in param_dict.items():
            setattr(self, k, v)

    def get_subj_gi(self):
        """
        Returns: subject's GI as int
        Excepts: ValueError
        """

        try:
            match = regex.search("gi\|(\d*)\|", self.subject)
        except KeyError:
            raise ValueError("Record has no 'subject' attribute")
        else:
            if match:
                return(int(match.group(1)))
            else:
                raise RuntimeError("GI wasn't found for subject: {}".format(self.subject))
