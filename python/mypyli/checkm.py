
import sys
import argparse
import json

def parse(checkm_fh):
    """ This method is for new-style checkm files. For old files, use legacy_parse. """

    for line in checkm_fh:
        bin_name, param_json = line[:-1].split("\t")

        params = json.loads(param_json.replace("'", '"'))

        yield CheckmRecord(name=bin_name, params=params)

def legacy_parse(checkm_fh):
    # replace any single quotes from CheckM with double quotes for JSON
    all_dict = json.loads(checkm_fh.read().replace("'", '"'))

    for entry in all_dict:
        yield CheckmRecord(name=entry, params=all_dict[entry])

def write_table(records, out="checkm_table.tsv"):
    """ This mostly shouldn't be used. CheckM offers an option that does this """
    with open(out, 'w') as OUT:
        OUT.write("\t".join(["name", "genome_size", "N50", "num_contigs", "longest_contig", "marker_lineage", "num_markers", "completeness", "contamination"]) + "\n")
        for record in records:
            OUT.write("\t".join([   str(record.name),
                                    str(record.genome_size),
                                    str(record.N50),
                                    str(record.num_contigs),
                                    str(record.longest_contig),
                                    str(record.marker_lineage),
                                    str(record.num_markers),
                                    str(record.completeness),
                                    str(record.contamination)
                                    ]) +"\n")

def write_record_classifications(records, out="checkm_record_classifications.txt"):
    classifications = {}
    for record in records:
        try:
            classifications[record.classify()].append(record.name)
        except KeyError:
            classifications[record.classify()] = [record.name]

    with open(out, 'w') as OUT:
        for category in ["quality", "sibling", "contaminated", "incomplete", "junk"]:
            OUT.write(category)
            try:
                records = classifications[category]
                OUT.write("({})\t".format(len(records)))
                OUT.write(" ".join(records))
            except KeyError:
                OUT.write("(0)\tNone")

            OUT.write("\n")

class CheckmRecord(object):
    """ Class to parse a Checkm Output file. Used when classifying CheckM bins """

    GOOD_COMP = 80
    MIN_COMP = 40
    MAX_CONTAM = 15
    SIB_CONTAM = 85

    def __init__(self, name, params):
        
        self.name = name

        self.genome_size = int(params.get('Genome size', 0))
        self.longest_contig = int(params.get('Longest contig', 0))
        self.marker_lineage = params.get('marker lineage', "")
        self.gc_content = float(params.get('GC', 0))
        self.num_scaffolds = int(params.get('# scaffolds', 0))
        self.num_contigs = int(params.get('# contigs', 0))
        self.completeness = float(params.get('Completeness', 0))
        self.contamination = float(params.get('Contamination', 0))
        self.num_markers = int(params.get('# markers', 0))
        self.num_marker_sets = int(params.get('# marker sets', 0))
        self.coding_density = float(params.get('Coding density', 0))
        self.N50 = float(params.get('N50 (contigs)', 0))

    def classify(self, good_comp=None, min_comp=None, max_contam=None, sib_contam=None):
        """ Classifies a record as "quality", "sibling", "contaminated", "incomplete", or "junk"

        Where:
            quality(the good) -> above completeness threshhold, below contamination threshhold
            sibling(the related) -> above completeness, high, high contam
            contaminated(the bad) -> above completeness threshhold, above contamination threshhold
            incomplete(the fragments) -> below completeness threshhold but still above some baseline and below contam thresh
            junk(the ugly) -> completeness below baseline or incomplete status but high contam

        The idea is that below some minimal completeness, there will be no differentiation between
        contigs that belong and those that do not. Chimeric bins will complement each other 
        and look more complete but really be a bag of junk.

        So, contamination threshholds will be a percent of completion.

        """
        
        # set omitted args to defaults
        if good_comp is None:
            good_comp = self.GOOD_COMP
        if min_comp is None:
            min_comp = self.MIN_COMP
        if max_contam is None:
            max_contam = self.MAX_CONTAM
        if sib_contam is None:
            sib_contam = self.SIB_CONTAM

        try: 
            contam_ratio = (self.contamination / self.completeness) * 100
        except ZeroDivisionError:
            return "junk"
        
        if self.completeness >= good_comp:
            if contam_ratio <= max_contam:
                return "quality"
            elif contam_ratio > max_contam and contam_ratio < sib_contam:
                return "contaminated"
            else:
                return "sibling"

        elif self.completeness >= min_comp:
            if contam_ratio <= max_contam:
                return "incomplete"
            else:
                return "junk"

        elif self.completeness < min_comp:
            return "junk"


class CheckmBin(object):
    """ Class to represent a bin from CheckM """

    def __init__(self, name):
        self.name = name
    
        self.completeness = None
        self.contamination = None

        self.contigs = []

    @staticmethod
    def parse_checkm_data(checkm_f):
        """ Parses CheckM's json files, yields tuples of bin_name, params """
        with open(checkm_f, 'r') as IN:
            for line in IN:
                bin_name, param_json = line[:-1].split("\t")

                params = json.loads(param_json.replace("'", '"'))

                yield bin_name, params

    @classmethod
    def from_marker_gene_stats(cls, marker_gene_stats_f):
        """ Yields CheckM bin objects for all the bins in a marker_gene_stats file """

        for bin_name, params in cls.parse_checkm_data(marker_gene_stats_f):
            checkm_bin = CheckmBin(bin_name)
            
            # go through the dict of contigs
            for contig, contig_data in params.items():
                checkm_contig = CheckmContig(contig)

                for key, value in contig_data.items():
                    pass

class CheckmContig(object):
    """ Represents a contig that is part of a CheckmBin """

    def __init__(self, name):
        self.name = name

        self.length = None

        self.markers = []


class CheckmMarker(object):
    """ Represents a marker and, optionally, its location """

    def __init__(self, name):
        self.name = name

        self.location = None







if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bin_stats_ext.tsv file", required=True)
    parser.add_argument("-o", help="output file", default="checkm_table.tsv")
    args = parser.parse_args()
    with open(args.i, 'r') as IN:
        records = parse(IN)
        write_record_classifications(records, args.o)
