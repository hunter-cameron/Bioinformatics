

from Bio import SeqIO
import sys
import argparse
import regex

def parse_fastx(fastx_f, reg_obj):
    """ Parses a fastx and yields only sequences that don't match the regex """
    # check if there is a "q" in the ext
    if "q" in fastx_f.split(".")[-1]:
        f_type = "fastq"
    else:
        f_type = "fasta"

    with open(fastx_f, 'r') as IN:
        for seq_obj in SeqIO.parse(IN, f_type):
            if not reg_obj.search(str(seq_obj.seq)):
                yield str(seq_obj.seq)


def find_mismatch(seq, reg_obj):
    
    match_obj = reg_obj.search(seq, regex.BESTMATCH)


def main(fastx_f, reg, anchor=0):
    """
    Current algorithm:

    Progressively add portions to the regex until it no longer matches. Wherever
    it stops matching is where the error is. 

    This has the potential to be slow.
    """
    # split regex into parts
    # confusing regex to return all matches of ( and then anything until ) 
    matches = regex.findall("(\([^)]*\))", reg)
    print("Groups found:", file=sys.stderr)
    print(matches, file=sys.stderr)
    
    # convert negative indexes into positive ones
    if anchor < 0:
        tmp = [x for x in range(len(matches))]
        anchor = tmp[anchor]
        print("Anchor converted: " + str(anchor))
    
    print("Anchored to group: {}".format(matches[anchor]), file=sys.stderr)

    # start with the anchor index
    cur_regex = matches[anchor]
    regexes = [regex.compile(cur_regex)]
    regex_to_group = {matches[anchor]: 0}     # will be used to relate the compiled regexes back to the group
    # first look to the left and build those regexes
    if anchor > 0:
        for indx in reversed(range(anchor)):
            cur_regex = matches[indx] + cur_regex
            regexes.append(regex.compile(cur_regex))
            regex_to_group[matches[indx]] = len(regexes) - 1

    # now look to the right and build those regexes
    if anchor < len(matches) - 1:
        cur_regex = matches[anchor]
        for indx in range(anchor + 1, len(matches)):
            cur_regex += matches[indx]
            regexes.append(regex.compile(cur_regex))
            regex_to_group[matches[indx]] = len(regexes) - 1

    # begin parsing seqs, submit the last regex for the initial check (the full thing)
    elem_counts = [0 for elem in regexes]   # sort of a hack to initialize the array to be as large as regexes
    for seq in parse_fastx(fastx_f, regexes[-1]):
        for indx, regx in enumerate(regexes):
            # increment the count of the index that doesn't match
            if not regx.search(seq):
                elem_counts[indx] += 1
                break

    print("group\terror_count")
    for group in matches:
        print("{}\t{}".format(group, elem_counts[regex_to_group[group]]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Checks a regex against the sequences in a fastx file and returns counts of the point at which the sequence stops matching the regex. NOTE: My experience has been that, to run this using bsub, you need to enclose the python ... in double quotes and the regex in single quotes.
            
            Example: bsub -o out -e err "python script.py -f fasta -r 'my_regex'"
            
In particular, I had trouble with a regex with a '|' character.""")
    parser.add_argument("-f", help="fasta or fastq file to parse (looks for a 'q' in the file extension for fastq)", required=True)
    parser.add_argument("-r", help="regex, be sure to capture all the groups you want to check using '()' Example: (elem1)(elem2)(...)(elemn)", required=True)
    parser.add_argument("-a", help="the group number from the regex to use to anchor the search (0 based)", default=0, type=int)

    args = parser.parse_args()

    main(args.f, args.r, args.a)
