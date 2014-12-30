

class SamParser(object):
    """ A class for parsing SAM files """
    def __init__(self, sam_f, aligned_only=False, mapq=1, mismatches=1):
        self.sam_f = sam_f
        self.qc = {'aligned_only': aligned_only, 'mapq': mapq, 'tag_filter': {'NM:i': mismatches}}
   
    def _keep_alignment(self, aln):
        # eliminate unmapped reads
        if self.qc['aligned_only']:
            if aln['flag'] == '4':
                return False

        # eliminate reads that are not unique enough based on mapping quality
        if self.qc['mapq']:
            if aln['mapq'] <= self.qc['mapq']:
                return False

        # eliminate reads based on their optional tags; assumes desire is to be below the value
        for tag, value in self.qc['tag_filter'].items():
            if aln.get(tag, 0):
                if aln[tag] <= value:
                    return False

        return True
                
             

    def parse_sam_file(self, iterate=True):
        """
        SAM generator that makes a SamAlignment for desired alignments in a file.
        Could be a class method in SamAlignment instead but that might seem counter-intuitive to call
        """

        with open(self.sam_f, 'r') as IN:
            for line in IN:
                if line.startswith("@"):
                    continue
                # store all the basic sam attributes
                aln_dict = {k: v for k, v in zip(['qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'], line[:-1].split('\t')[:10])}

                aln_dict['pos'] = int(aln_dict['pos'])
                aln_dict['mapq'] = int(aln_dict['mapq'])

                if self._keep_alignment(aln_dict):

                    if iterate:
                        yield aln_dict
                    else:
                        alignments.append(aln_dict)


