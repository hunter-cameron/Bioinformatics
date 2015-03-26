

from Bio import SeqIO
from Bio.Seq import Seq
import difflib
import sys
import regex
import argparse

def mine_seqs(in_f, out_f):
    with open(in_f, 'r') as IN, open(out_f, 'w') as OUT:
        rec = 1
        for record in SeqIO.parse(IN, 'fasta'):
            
            genomic_frag = _fuzzy_match(record.seq)
            if genomic_frag:
                OUT.write(">{}; genomic_len={}".format(record.id, len(genomic_frag)) + "\n")
                OUT.write(genomic_frag + "\n")

def _fuzzy_match(seq, transposon, errors=2):

    # matches 2-4 T, then 10-30 nucs, then the transposon allowing errors
    reg = "[NT]{3,4}([ATGC]{10,30})(:?" + transposon + "){e<=" + str(errors) + "}"
    
    #reg = ".*([ATGC]{14})(:?" + transposon + "){e<=" + str(errors) + "}"
    match = regex.match(reg, str(seq))

    #print(("match", match))

    if match:
        return match.group(1)
    else:
        return ""

def single_miner(read_f, out_f, reg):
    pass

def paired_miner(r1_f, r2_f, out_f, reg1, reg2):
    """ Parses the fwd and rev reads to yield consensus seqs """
    with open(r1_f, 'r') as R1, open(r2_f, 'r') as R2, open (out_f, 'w') as OUT:
        missing_r1 = 0
        missing_r2 = 0
        missing_both = 0
        short_conseq = 0
        r1_len = {}
        r2_len = {}
        subseq_len = {}
        for rec1, rec2 in zip(SeqIO.parse(R1, 'fasta'), SeqIO.parse(R2, 'fasta')):
            r1_genom = get_genomic_frag(rec1.seq, reg1)
            r2_genom = get_genomic_frag(rec2.seq, reg2)

            if (not r1_genom) and (not r2_genom):
                missing_both += 1
                continue
            elif not r1_genom:
                missing_r1 += 1
                continue
            elif not r2_genom:
                missing_r2 += 1
                continue

            #print(("s2ge", r2_genom))
            r1_len[len(r1_genom)] = r1_len.get(len(r1_genom), 0) + 1
            r2_len[len(r2_genom)] = r2_len.get(len(r2_genom), 0) + 1

            r2_genom_rc = Seq(r2_genom).reverse_complement()

            #print(("seq1", r1_genom))
            #print(("seq2", str(r2_genom_rc)))


            subseq = longest_common_substring(r1_genom, r2_genom_rc)
            subseq_len[len(subseq)] = subseq_len.get(len(subseq), 0) + 1
            #print(("subs", subseq))
            #print()
            #sys.exit()
            
            #if len(subseq) < 5:
                #print(("s2ge", r2_genom))
                #print(("seq1", r1_genom))
                #print(("seq2", str(r2_genom_rc)))
                #print(("subs", subseq))
                #print()


            if len(subseq) < 15 or len(subseq) > 17:
                short_conseq += 1
                continue
            OUT.write(">{}; genomic_len={}".format(rec1.id, len(subseq)) + "\n")
            OUT.write(subseq + "\n")

    print("Both bad: {}".format(str(missing_both)))
    print("R1 bad: {}".format(str(missing_r1)))
    print("R2 bad: {}".format(str(missing_r2)))
    print("Subseq too short/long: {}".format(str(short_conseq)))

    names = ["R1 length", "R2 length", "Subseq Length"]
    for indx, seq_type in enumerate([r1_len, r2_len, subseq_len]):
        print(names[indx])
        for key in sorted(seq_type):
            if int(key) >= 15 and int(key) <= 17:
                print("  {}\t{}".format(str(key), str(seq_type[key])))

def longest_common_substring(seq1, seq2):
    """ Longest common substr dynamic programming style """

    # make a matrix of seq2 rows and seq1 columns
    matr = [[0]*(len(seq1) + 1) for i in range(len(seq2) + 1)]

    length = 0
    end = 0
    for row in range(1, len(seq2) + 1):
        for col in range(1, len(seq1) + 1):
            if seq2[row-1] == seq1[col-1]:
                matr[row][col] = matr[row-1][col-1] + 1

                # if there is a non-0 score, treat this as the optimal
                # solution (this can be changed in later iterations)
                if matr[row][col] > length:
                    length = matr[row][col]
                    end = row

            else:
                matr[row][col] = 0

    #for row in matr:
    #    print(row)

    # return the subseq derriving the start as the end minus the length
    return str(seq2)[end-length:end]

def merge_with_n(seq1, seq2):
    pass


def get_genomic_frag(seq, reg):
    m = regex.match(reg, str(seq))

    if m:
        return m.group(1)
    else:
        return ""


def _n_minus_T(seq, n=25):
    """ Matches the first n characters and then strips off the starting "T"s """
    match = regex.match("T*([ACTG]*)", str(seq)[:n])

    if match:
        return match.group(1)
    else:
        return ""

if __name__ == "__main__":
    #mine_seqs(sys.argv[1], sys.argv[2])
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-r1", help="fasta of r1 reads")
    parser.add_argument("-r2", help="fasta of r2 reads")
    parser.add_argument("-e", help="number of errors to allow in anchor sequence", default="0")
    parser.add_argument("-out", help="out file for sequences", required=True)
    args = parser.parse_args()

    #TRANSPOSON="TAACAGGTTGGATGATAAGTCCCCGGTCT"
    #TRANSPOSON="TAACAGGTT"
    
    
    #reg1 = "C{0,1}TGTTA([ACTG]{10,30})(AAAAGATCGGAAGAGCACACGTCTGAACTCC){e<=" + args.e + "}"
    #reg1 = "C{0,1}TGTTA([ACTG]{10,30})(AAAAGAT){e<=" + args.e + "}"
    
    # insertion may leave TA in genome so take that off filtering
    TRANSPOSON="ACAGGTT"
    reg1 = "C{0,1}TGT([ACTG]{10,30})(AAAAGAT){e<=" + args.e + "}"
    reg2 = "[NT]{3,4}([ATGC]{10,30})(:?" + TRANSPOSON + "){e<=" + args.e + "}"

    if args.r1 and args.r1:
        paired_miner(args.r1, args.r2, args.out, reg1, reg2)

    elif args.r1:
        single_miner(args.r1, args.out, reg1)
    elif args.r2:
        single_miner(args.r2, args.out, reg2)

