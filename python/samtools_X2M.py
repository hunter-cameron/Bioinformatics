
import sys

import argparse
import re

def X2M(sam_f, out_f):
    with open(sam_f, 'r') as IN, open(out_f, 'w') as OUT:
        # compile the regex since we are going to use it a bunch
        reg = re.compile("\d+\D")
        for line in IN:
            
            # write any header lines
            if line.startswith("@"):
                OUT.write(line)
                continue
            

            fields = line[:-1].split("\t")

            cigar = fields[5]
            cigar = cigar.replace("=", "M")
            cigar = cigar.replace("X", "M")

            #print(cigar)
            m = reg.findall(cigar)
            # capture all the cigar groups (each with some number of digits and a single non-digit)
            #print(m)
            if m:
                cigar_tups = []
                for group in m:
                    #print(group)
                    op = group[-1]
                    #print(op)

                    tup = (int(group[:-1]), op)
                    #print(tup)
                    # add the groups if needed
                    if cigar_tups:
                        if cigar_tups[-1][1] == tup[1]:
                            cigar_tups[-1] = ((cigar_tups[-1][0] + tup[0], cigar_tups[-1][1]))
                    else:
                        cigar_tups.append(tup)
                
                fields[5] = "".join([str(count) + op for count, op in cigar_tups])

            else:
                fields[5] = "*"

            OUT.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts SAM cigar formats that use =/X instead of M into M format")
    parser.add_argument("-sam", help="the SAM file to convert", required=True)
    parser.add_argument("-out", help="the filename for the converted SAM file", default="M_format.sam")

    args = parser.parse_args()
    
    X2M(args.sam, args.out)
