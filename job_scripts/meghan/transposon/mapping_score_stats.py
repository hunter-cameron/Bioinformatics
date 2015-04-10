
"""
I want to read in all the SAM files and figure out if length
(of alignment) or mismatches is the main contributor in reads with a low
quality score.
"""

from mypyli import samparser
import sys

def make_count_table(sam_files):
    counts_table = {}
    for sam_f in sam_files:
        with open(sam_f) as IN:
            for record in samparser.parse(IN, aligned_only=True):
                if not record.flag in [0, 16]:      # primary alignments only
                    continue

                mismatches = record.mismatches
                # true length shouldn't count soft clipped bases
                length = record.length - record.soft_clipped
                mapq = record.mapq


                try:
                    counts_table[mapq][length][mismatches] += 1
                except KeyError:
                    try:
                        counts_table[mapq][length][mismatches] = 1
                    except KeyError:
                        try:
                            counts_table[mapq][length] = {mismatches: 1}
                        except KeyError:
                            counts_table[mapq] = {length: {mismatches: 1}}
                    
    return counts_table




if __name__ == "__main__":
    sam_files = sys.argv[1:]

    counts_table = make_count_table(sam_files)

    for qual in sorted(counts_table):
        print("Qual: {}".format(str(qual)))
        for length in sorted(counts_table[qual]):
            print("\tLen: {}".format( str(length)))
            for mismatch in sorted(counts_table[qual][length]):
                print("\t\tMismatch {}: {}".format(str(mismatch), str(counts_table[qual][length][mismatch])))
            
        print()
