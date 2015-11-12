#!/usr/bin/env python

import argparse
import sys
import os

description = """
Author: Hunter Cameron
Date: Oct 19, 2015
Description:
Small utility for comparing two lists. It outputs the elements that 
are only present in one of the lists and elements that are present in
both. Allows each list of elements to be output to a different file
handle or standard out/error.
"""



## READ ARGS
parser = argparse.ArgumentParser(description=description)
parser.add_argument("l1")
parser.add_argument("l2")
parser.add_argument("-l1_out", help="path or stdout/stderr to write elements that are only in l1")
parser.add_argument("-l2_out", help="path or 'stdout'/'stderr' to write elements that are only in l2")
parser.add_argument("-both_out", help="path or 'stdout'/'stderr' to write elements that are in both lists")
args = parser.parse_args()

itm_dict = {}
for indx, lst in enumerate([args.l1, args.l2]):
    indx += 1
    with open(lst, 'r') as IN:
        for line in IN:
            itm = line.rstrip()

            # check if item already in dict
            try:
                # if it was already found in this list, raise error
                if itm_dict[itm] >= indx:
                    raise ValueError("Duplicate name(s) in list {}. Example: '{}'".format(str(indx), itm))
                else:
                    itm_dict[itm] += indx

            # add new item to list
            except KeyError:
                itm_dict[itm] = indx

"""
parse the lists; 
in l1 but not l2 = 1
in l2 but not l1 = 2
in both = 3
"""
l1 = []
l2 = []
both = []
for k, v in itm_dict.items():
    if v == 1:
        l1.append(k)
    elif v == 2:
        l2.append(k)
    elif v == 3:
        both.append(k)
    else:
        print((k, v))
        raise Exception("This should have never happened...")

def open_output_fh(arg):
    """ Opens and returns a FH depending on the argument supplied """

    if arg:
        if arg == "stdout":
            return sys.stdout
        elif arg == "stderr":
            return sys.stderr
        else:
            try:
                fd = os.open(arg, os.O_WRONLY | os.O_CREAT | os.O_EXCL)

            except OSError:
                print("Cowardly refusing to overwrite path {}. Exiting.".format(arg))
                sys.exit("Exited with error status.")

            return os.fdopen(fd, 'w')

    else:
        return sys.stdout

## WRITE THE OUTPUTS TO THE APPROPRIATE PLACES

# write both
both_fh = open_output_fh(args.both_out)
both_fh.write("{} items in both lists:\n".format(str(len(both))))
for i in both:
    both_fh.write(i + "\n")

# write l1 only
l1_fh = open_output_fh(args.l1_out)
# print a blank line if this is a fh we have already written to
if args.l1_out == args.both_out:
    l1_fh.write("\n")
l1_fh.write("{} items in List 1 only:\n".format(str(len(l1))))
for i in l1:
    l1_fh.write(i + "\n")

# write l2 only
l2_fh = open_output_fh(args.l2_out)
# print a blank line if this is a fh we have already written to
if args.l2 == args.both_out or args.l2 == args.l1_out:
    l2_fh.write("\n")
l2_fh.write("{} items in List 2 only:\n".format(str(len(l2))))
for i in l2:
    l2_fh.write(i + "\n")


## NOTE -- all FH are still open, they will all close when the program ends.
## This is done this way as a shortcut to avoid closing stdout/stderr.
