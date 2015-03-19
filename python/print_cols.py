
import sys
import argparse




parser = argparse.ArgumentParser("Prints columns from a file as specified by the indices")
parser.add_argument("-file", help="initial file", required=True)
parser.add_argument("-indx", help="indices to print separated by commas. Ex: 0,2,4")
parser.add_argument("-sep", help="delimiter in the file", default="\t")
args = parser.parse_args()

indices = args.indx.split(",")
with open(args.file, 'r') as IN:
    for line in IN:
        elems = line[:-1].split(args.sep)

        to_print = []
        for indx in indices:
            to_print.append(elems[int(indx)])

        print(args.sep.join(to_print))




