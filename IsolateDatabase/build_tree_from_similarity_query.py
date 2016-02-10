
import argparse
import sys
from Bio.Phylo import TreeConstruction
from Bio import Phylo

import logging
logging.basicConfig()
LOG = logging.getLogger(__name__)
LOG.setLevel("INFO")



def build_tree(distance_matrix, alg="nj"):

    # make the calculator and get the tree constructor. 
    # truly the calculator serves no purpose, it is demanded as an argument,
    # however, the methods that build the tree and instance methods that should
    # really be class methods, so I will call those directly
    calculator = TreeConstruction.DistanceCalculator()
    tc = TreeConstruction.DistanceTreeConstructor(calculator, method="nj")

    if alg == "nj":
        return tc.nj(distance_matrix)

    elif aln == "upgma":
        return tc.upgma(distance_matrix)

    else:
        raise ValueError("alg must be either 'nj' or 'upgma'")

def build_distance_matrix(nested_dict):
    """ Builds a _DistanceMatrix object from a nested dict representation of a matrix """

    names = list(nested_dict.keys())

    # build a lower-triangular matrix
    matrix = []
    for n1 in names:
        row = []
        for n2 in names:
            if n1 == n2:
                row.append(0)
                matrix.append(row)
                break

            # check the matrix is symmetric, otherwise make it symmetric
            v1 = nested_dict[n1][n2]
            v2 = nested_dict[n2][n1]
            if v1 == v2:
                value = v1
            else:
                # average values to make symmetric
                value = (v1 + v2) / 2.0

            row.append(value)

    return TreeConstruction._DistanceMatrix(names, matrix)

def parse_input(input_fh):
    """ Parses the input list and returns a distance matrix in nested dict format  """

    nested_dict = {}
    for line in input_fh:

        line = line.rstrip()

        query, reference, similarity = line.split("\t")
            
        # try to convert similarity to a float
        # this may fail if there is a header line or if columns are in the wrong order
        try:
            similarity = float(similarity)
        except ValueError:
            LOG.warning("Value in similarity field '{}' could not be converted to float. Skipping entry.".format(similarity))
            continue

        distance = 1 - similarity
        
        # add the value to the nested dict
        try:
            nested_dict[query][reference] = distance
        except KeyError:
            nested_dict[query] = {reference: distance}

    # check that there is an entry for each combination
    queries = set(nested_dict.keys())
    tmp_set = set()
    references = set()
    for value in nested_dict.values():
        references.update(value.keys())

    if references != queries:
        raise ValueError("Queries are not the same as references. Some pairwise combinations were forgotten.\nQueries: {}\n\nReferences: {}".format(str(queries), str(references)))

    # now that we know the two sets are the same we want to make sure each combination has been made
    matrix_is_invalid = False
    for q1 in queries:
        for q2 in queries:
            # any key errors in lookups mean missing fields
            try:
                nested_dict[q1][q2]
            except KeyError:
                # set matrix as invalid and then do an additional step to give more informative warnings
                matrix_is_invalid = True
                try:
                    nested_dict[q1]
                except KeyError:
                    LOG.warning("Genome '{}' is not present as a query.")
                    break

                LOG.warning("Genome '{}' is not present as a reference for '{}'.".format(q2, q1))

    # fail if invalid, otherwise return
    if matrix_is_invalid:
        raise ValueError("One or more warnings were encountered while checking the input. Aborting.")
    else:
        return nested_dict



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Builds a tree from the output of a query from the GenomeSimilarity table of the isolate database. Reads input from stdin and expects 3 columns 'query, reference, similarity'. The query and reference fields will be used as the names in the tree.")
    parser.add_argument("input", help="results from the MySql query", nargs="?", type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument("-alg", help="the tree building algorithm", choices=['nj', 'upgma'], default="nj")
    parser.add_argument("-format", help="the format to write the tree in", choices=["newick"], default="newick")
    args = parser.parse_args()

    nested_dict = parse_input(args.input)

    distance_matrix = build_distance_matrix(nested_dict)

    # build the tree and write it to stdout
    tree = build_tree(distance_matrix, alg=args.alg)
    Phylo.write(tree, sys.stdout, args.format)
