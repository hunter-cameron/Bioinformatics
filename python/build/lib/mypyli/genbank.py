#!/usr/bin/env python3

"""
Author: Hunter Cameron
Date: 11/31/2014
Updated: 11/31/2014
Description:
    A module for working with GenBank files that uses BioPython under the hood.

"""

from Bio import SeqIO
from taxstring import TaxString

class GenBank(object):
    """ 
    Class for broad processing of GenBank files.
    USAGE: GenBank(gbk_file=my_file)
    """

    def __init__(self, gbk_file=""):
        pass

    def get_tax_from_gbk(gbk_file):
        """ Returns a TaxString object from the info in a gbk file (So far just tested with ones from JGI """
        tax_lines = []
        append = False
        with open(gbk_file, 'r') as IN:
            for line in IN:
                if line.startswith("  ORGANISM"):
                    append = True
                if line.startswith("FEATURES"):
                    break

                if append:
                    tax_lines.append(line)

        # remove the first line because it's junk and remove leading and trailing spaces from lines
        clean_lines = [line.strip() for line in tax_lines[1:]]

        # concatenate the lines and remove that period they put at the end
        tax_string = " ".join(clean_lines)
        tax_string = tax_string.replace(".", "")

        #rank_tax = []
        #for rank, tax in zip(["domain", "phylum", "class", "order", "family", "genus", "species"], tax_string.split("; ")):
        #    rank_tax.append((rank, tax))

        
        return TaxString(tax_string)

