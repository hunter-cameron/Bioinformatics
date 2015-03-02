
import sys
import argparse
import json

def parse(checkm_fh):
    # replace any single quotes from CheckM with double quotes for JSON
    all_dict = json.loads(checkm_fh.read().replace("'", '"'))

    for entry in all_dict:
        yield CheckmRecord(name=entry, params=all_dict[entry])


def write_table(checkm_fh, out="checkm_table.tsv"):
    with open(out, 'w') as OUT:
        OUT.write("\t".join(["name", "genome_size", "N50", "num_contigs", "longest_contig", "marker_lineage", "num_markers", "completeness", "contamination"]) + "\n")
        for record in parse(checkm_fh):
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


class CheckmRecord(object):
    """ Class to parse a Checkm Output file """

    def __init__(self, name, params):
        
        self.name = name

        self.genome_size = params.get('Genome size', "")
        self.longest_contig = params.get('Longest contig', "")
        self.marker_lineage = params.get('marker lineage', "")
        self.gc_content = params.get('GC', "")
        self.num_scaffolds = params.get('# scaffolds', "")
        self.num_contigs = params.get('# contigs', "")
        self.completeness = params.get('Completeness', "")
        self.contamination = params.get('Contamination', "")
        self.num_markers = params.get('# markers', "")
        self.num_marker_sets = params.get('# marker sets', "")
        self.coding_density = params.get('Coding density', "")
        self.N50 = params.get('N50 (contigs)', "")
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input bin_stats_ext.tsv file", required=True)
    parser.add_argument("-o", help="output file", default="checkm_table.tsv")
    args = parser.parse_args()
    with open(args.i, 'r') as IN:
        write_table(IN, args.o)
