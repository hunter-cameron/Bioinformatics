#!/usr/bin/env python3

import argparse
from lxml import etree
from pprint import pprint

from lucidSubject import Subject
from lucidQuery import Query
from lucidAlignment import Alignment
import lucidArtist


def parse_BLAST(args):
    """
    Selects which parsing method to use. Currently only format supported is XML.
    :param args:
    :return:
    """
    if args["input_format"] == "XML":
        hit_hash = parse_XML(args["input"], args["qc_params"])
        # pprint(hit_hash)

    return dict_to_subject_list(hit_hash)

# TODO consider rewriting this with biopython. Also, add QC measures to remove redundancy.
def parse_XML(xml_file, qc_params):
    """ Takes a BLAST output XML file and parses it to return a nested dictionary of hits.
    Dictionary ordered: subject->"queries"->query->"hsps"->hsp->{stats}

    Need to add some QC to shrink the potential dictionary size and remove useless hits.

    """

    print("Parsing XML...")

    master_hash = {}
    for event, element in etree.iterparse(xml_file, tag="Iteration"):
        query = element.findtext("Iteration_query-def")
        qry_length = int(element.findtext("Iteration_query-len"))

        iteration_hits = element.find("Iteration_hits")
        for hit in iteration_hits.findall("Hit"):
            subject = hit.findtext("Hit_def")
            subj_length = int(hit.findtext("Hit_len"))

            if subject not in master_hash:
                master_hash[subject] = {"length": subj_length, "queries": {}}

            if query not in master_hash[subject]["queries"]:
                master_hash[subject]["queries"][query] = {"length": qry_length, "hsps": []}

            hit_hsps = hit.find("Hit_hsps")
            for hsp in hit_hsps.findall("Hsp"):
                hsp_hash = {"score": float(hsp.findtext("Hsp_score")),      # floats may not be precise enough, may want to import decimal.Decimal (probably more storage required though)
                            "evalue": float(hsp.findtext("Hsp_evalue")),
                            "qry_start": int(hsp.findtext("Hsp_query-from")),
                            "qry_end": int(hsp.findtext("Hsp_query-to")),
                            "subj_start": int(hsp.findtext("Hsp_hit-from")),
                            "subj_end": int(hsp.findtext("Hsp_hit-to")),
                            #"identity": hsp.findtext("Hsp_identity"),      #not percent ID
                            "length": int(hsp.findtext("Hsp_align-len")),
                            "midline": hsp.findtext("Hsp_midline"),
                            }

                qc_test = hit_qc(hsp_hash, subj_length, qc_params)
                if qc_test:
                    master_hash[subject]["queries"][query]["hsps"].append(qc_test)

        element.clear()

    return remove_hitless_queries(master_hash)


def hit_qc(hit, subj_length, qc_params):
    """ Some simple quality control to determine if the hit should be stored.
    Also checks orientation of the hit and assigns a flag.

    Returns None to delete the hit and returns the edited hit if the hit should be stored.

    Better way to do this, i.e. without passing s_length?
    """

    #check length
    if qc_params["min_length"] < 1:         # if < 1 use min_length as percent of subject -- need option to use % of query for BLASTing reads or short seqs?
        if hit["length"] < subj_length * qc_params["min_length"]:
            return None
    else:
        if hit["length"] < qc_params["min_length"]:
            return None

    #check score
    if hit["score"] < qc_params["min_score"]:
        return None

    #check evalue
    if hit["evalue"] > qc_params["max_evalue"]:
        return None

    #check identity -- this might be a costly check if the sequences get large enough
    # also may not work if there are gaps
    hit["perc_identity"] = hit["midline"].count("|") / hit["length"]
    if hit["perc_identity"] < qc_params["min_perc_identity"]:
        return None
    else:
        del hit["midline"]

    ###
    # End QC checking
    ###

    # check for reverse hit
    if hit["subj_start"] > hit["subj_end"]:
        hit["reverse"] = True

        # reverse the subject start and end -- should I reverse the query start and end too? Not currently at least.
        end = hit["subj_start"]
        start = hit["subj_end"]
        hit["subj_start"] = start
        hit["subj_end"] = end

    else:
        hit["reverse"] = False

    return hit


def remove_hitless_queries(dictionary):
    """ Deletes queries that had no hits. Could make it work for subjects but might want to display blank ones.

    :param dictionary: A dictionary
    :return:
    """

    queries_to_del = {}
    for subject in dictionary:
        for query in dictionary[subject]["queries"]:
            if not dictionary[subject]["queries"][query]["hsps"]:
                if subject in queries_to_del:
                    queries_to_del[subject].append(query)
                else:
                    queries_to_del[subject] = [query]

    for subject in queries_to_del:
        for query in queries_to_del[subject]:
            del dictionary[subject]["queries"][query]

    return dictionary


def dict_to_subject_list(dictionary):
    """
    Takes the dictionary generated from parse XML and converts it to a list of lucidSubject objects complete with alignments.
    :param dictionary:
    :return:
    """
    subjects = []
    for subject in dictionary:
        subj = Subject({'name': subject, 'length': dictionary[subject]['length']})

        for query in dictionary[subject]['queries']:
            queries = dictionary[subject]['queries']

            qry = Query({'name': query, 'length': queries[query]['length']})

            for hsp in queries[query]['hsps']:
                hsp['subject'] = subj
                hsp['query'] = qry

                align = Alignment(hsp)

                # I should do some sort of check for original alignments
                # ex. if alignment overlaps 95%+ with an existing alignment is isn't original
                # ex2. if starts and ends are both within 100bp of each other

                subj.add_alignment(align)

        subjects.append(subj)

    return subjects


def plot_subjects(subjects):
    scaffolds = {}
    [scaffolds.setdefault(subject.name, subject.scaffold) for subject in subjects]

    artist = lucidArtist.MultiPlotScaffold(scaffolds=scaffolds, num_panels=3)
    artist.plot_scaffolds()




def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", help="Input a file of BLAST results.", required=True)
    parser.add_argument("-fmt", help="The format of the ", choices=["XML", ], default="XML")
    parser.add_argument("-min_perc_id", help="The minimum percent identity (decimal) to retain a high-scoring pair.", type=float, default=.90)
    parser.add_argument("-min_score", help="The minimum score to retain a high-scoring pair.", type=float, default=0.0)
    parser.add_argument("-min_length", help="The minimum length to retain a high-scoring pair. Values below 1 will be treated as a percentage of the subject length while values 1 and above will be treated as a number of base pairs.", default=1000.0, type=float)
    parser.add_argument("-max_evalue", help="The maximum evalue to retain a high-scoring pair.", type=float, default=.05)

    args = parser.parse_args()

    # reformat the args into a dictionary
    d_args = {"input": args.i,
              "input_format": args.fmt,
              "qc_params": {"min_perc_identity": args.min_perc_id,
                            "min_score": args.min_score,
                            "min_length": args.min_length,
                            "max_evalue": args.max_evalue,
                            },
              }

    return d_args

if __name__ == "__main__":

    args = parse_args()

    subjects = parse_BLAST(args)

    for subject in subjects:
        subject.make_scaffold()

    plot_subjects(subjects)