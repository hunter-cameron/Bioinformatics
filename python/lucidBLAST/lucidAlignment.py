__author__ = 'Hunter Cameron'


class Alignment(object):

    def __init__(self, args):
        self.name = "{}_{}".format(args["subject"], args["query"])
        self.subject = args["subject"]
        self.query = args["query"]
        self.length = args["length"]
        self.qry_start = args["qry_start"]
        self.qry_end = args["qry_end"]
        self.subj_start = args["subj_start"]
        self.subj_end = args["subj_end"]
        self.perc_identity = args["perc_identity"]
        self.score = args["score"]
        self.evalue = args["evalue"]
        self.reverse = args["reverse"]

    def query_overlaps(self, alignment2):
        """ Checks for overlap between two alignments
        :param self
        :param alignment2:
        :return:
        """

        return (self.subj_end >= alignment2.subj_start) and (self.subj_start <= alignment2.subj_end)

