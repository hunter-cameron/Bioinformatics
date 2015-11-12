#!/usr/bin/env python3

"""
Author: Hunter Cameron
Date: 11/31/2014
Updated: 11/31/2014
Description:
    A module for working with GenBank files that uses BioPython under the hood.

"""

from Bio import SeqIO

class GenBank(object):
    """ 
    Class for broad processing of GenBank files.
    USAGE: GenBank(gbk_file=my_file)
    """

    def __init__(self, gbk_f):
        self.gbk_f = gbk_f

        # get some basic info
        with open(self.gbk_f, 'r') as IN:
            for record in SeqIO.parse(IN, "genbank"):
                self.organism = record.annotations["organism"]
                self.taxonomy = record.annotations["taxonomy"]
                break
