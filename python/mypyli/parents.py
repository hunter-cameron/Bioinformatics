#!/usr/bin/env python3

"""
Author: Hunter Cameron
Date: 11/31/2014
Updated: 11/31/2014
Description:
    A module for containing parent classes that are never used but kept only for the purpose of inheritance.

"""

class Sequence(object):
    def __init__(self, header="", seq=""):

        for param, value in {'header': header, 'seq': seq}:
            if value == "":
                raise AssertionError("{} must not be empty.".format(param))

        self.header = header
        self.seq = seq

    def get_header(self):
        return self.header

    def get_seq(self):
        return self.seq
