__author__ = 'Hunter Cameron'


class Parent(object):

    def __init__(self, name, length):
        self.name = name
        self.length = length
        self.alignments = []

    def add_alignment(self, alignment):
        self.alignments.append(alignment)

    def sort_alignments_by_length(self):
        """  Sorts alignments by length and stores it as the alignments attribute.
        :return list:
        """
        self.alignments = sorted(self.alignments, key=lambda x: x.length, reverse=True)      # use a lambda b/c key expects a function
        return self.alignments      # also return the list in case that is what people expected
