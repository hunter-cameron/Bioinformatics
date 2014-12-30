__author__ = 'Hunter Cameron'

# TODO Consider some sort of filtering mechanism to select only best subjects when BLASTing to a (large) database
# TODO Add option to filter redundant query hits. (So redundant hits don't stack up if the user mostly wants to see how well they have the subject genome covered)

# TODO Add 'real' alignment mode where the query length and orientation is taken into acoount when creating scaffolds. (Useful to align assemblies where you would expect certain HSPs to be in a particular relationship with one another ex. on the same contig.)



from lucidParent import Parent
import lucidArtist

class Subject(Parent):

    def __init__(self, args):
        Parent.__init__(self, args["name"], args["length"])
        self.scaffold = []

    def make_scaffold(self):
        """  Places the alignments on the subject sequence with the largest alignments closest to the reference.
        :return:
        """

        alignments = super(Subject, self).sort_alignments_by_length()

        scaffold = [self]       # subject is at 0 for scaffold

        for alignment in alignments:
            _check_overlap_by_tier(scaffold, 1, alignment)

        self.scaffold = scaffold


# TODO check this method again. As a first go at recursion, it could probably be done better.
def _check_overlap_by_tier(scaffold, tier, alignment):
    """
    The workhorse or the make scaffold function. Recursively makes tiers of non-overlapping alignments.
    :param scaffold:
    :param tier:
    :param alignment:
    :return:
    """
    if tier < len(scaffold):
        for indx, existing in enumerate(scaffold[tier]):
            # print(tier, alignment.query.name, "\t", existing.query.name)
            if alignment.query_overlaps(existing):
                # print("overlap")
                _check_overlap_by_tier(scaffold, tier + 1, alignment)
                return      # return required to get out of recursion

            # if doesn't match the last alignment for the tier, add it in
            elif not alignment.query_overlaps(existing) and indx == len(scaffold[tier]) - 1:
                # print("No overlap, adding to tier {}".format(tier))
                scaffold[tier].append(alignment)        # if it gets here, there is no overlap in the tier
                return      # return required to get out of recursion

    else:       # alignment overlapped on all tiers made so far, need to make a new one
        scaffold.append([alignment])

