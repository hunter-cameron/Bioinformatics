
import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def find_good_cols(msa, frac_trash):
    """ 
    Returns an array of ranges (as tuples) to keep by counting gaps and Ns at each 
    position and comparing their fraction to the maximum amount of trash allowed.
    """
    
    num_seqs = len(msa)
    
    # iterate over each column of the MSA
    to_keep = []
    start = None
    for i in range(len(msa[0])):
        col = msa[:, i]

        # count gaps and Ns
        gaps = col.count("-")
        Ns = col.count("N")

        # check against frac_trash
        if (gaps + Ns) / num_seqs < frac_trash:
            # check if we need to open a new range
            if start is None:
                start = i

        else:
            # check if we need to end an existing range
            if start:
                # this position is 1 past where we want to keep but that's good because 
                # indexing is half open
                to_keep.append((start, i))
                start = None
    
    # end the last range if there is one
    if start:
        # add an end anchor 1 past the last index.
        # technically out of bounds but it will let us get the end
        to_keep.append((start, len(msa[0])))

    return to_keep

def extract_cols(old_msa, to_keep):
    """ Extracts to_keep regions from the old MSA and returns the new one"""

    # build up a new MSA by concatenating all the ranges from to_keep
    new_msa = None
    for start, end in to_keep:
        if new_msa is None:
            new_msa = old_msa[:, start:end]

        else:
            new_msa += old_msa[:, start:end]

    return new_msa
    
def filter_seqs(msa, frac_trash):
    """ Returns a new MSA with sequences above the frac_trash threshhold removed """

    seq_len = len(msa[0])

    to_keep = [] 
    for index, record in enumerate(msa):
        gaps = record.seq.count("-")
        Ns = record.seq.count("N")

        if (gaps + Ns) / seq_len < frac_trash:
            to_keep.append(record)
        else:
            print("Removing seq {seqid} trash={trash}%".format(seqid=record.id, trash=(gaps + Ns) * 100 / seq_len))
    return MultipleSeqAlignment(to_keep)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Removes columns from a multiple sequence alignment with greater than a given fraction of gaps or 'N's (trash)")
    parser.add_argument("-msa", help="multiple sequence alignment in FASTA format", required=True)
    parser.add_argument("-col_filt", help="remove columns with this decimal proportion or greater of trash", required=True, type=float)
    parser.add_argument("-seq_filt", help="afer col_filt, remove seqs with this decimal proportion or greater of trash [%(default)s]", default=.1, type=float)
    parser.add_argument("-out", help="path to output the trimmed MSA file [%(default)s]", default="trimmed_MSA.fa")

    args = parser.parse_args()

    # read the MSA
    with open(args.msa, 'r') as IN:
        msa = AlignIO.read(IN, "fasta")

    to_keep = find_good_cols(msa, args.col_filt)
    new_msa = extract_cols(msa, to_keep)
    new_msa = filter_seqs(new_msa, args.seq_filt)

    # write the new MSA
    with open(args.out, 'w') as OUT:
        AlignIO.write(new_msa, OUT, "fasta")


