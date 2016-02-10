
import argparse

from Bio.Phylo import TreeConstruction
from Bio import Phylo



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

def build_distance_matrix(nested_dict, names_dict):
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

            # check the matrix is symmetric
            v1 = nested_dict[n1][n2]
            v2 = nested_dict[n2][n1]
            if v1 == v2:
                value = v1
            else:
                # average values to make symmetric
                value = (v1 + v2) / 2.0

            row.append(value)

    # update names if required
    new_names = []
    if names_dict:
        for name in names:
            new_names.append(names_dict[name])

        names = new_names

    return TreeConstruction._DistanceMatrix(names, matrix)

def parse_input_matrix(input_f):
    """ Parses the input matrix and returns a nested dict representation of the matrix """

    with open(input_f, 'r') as IN:
        headers = []
        nested_dict = {}
        for line in IN:

            line = line.rstrip()

            # read in the header line if this is the first line
            if headers == []:
                headers = line.split("\t")[1:]
                continue

            elems = line.split("\t")
            
            # initialize this row in the dict
            genome = str(elems[0])
            nested_dict[genome] = {}

            # fill out all the columns
            for header, value in zip(headers, elems[1:]):
                # get distance from similarity score
                nested_dict[genome][str(header)] = 1 - float(value)

        return nested_dict

def create_names_dict(rename, names):
    """ 
    Function that changes names. 
    
    Takes a list of names and returns a dict that maps each names to a new name.

    This function can be edited by the user.
    """
    names_dict = {}
    with open(rename, 'r') as IN:
        for line in IN:
            old_name, new_name = line.rstrip().split("\t")
            names_dict[old_name] = new_name

    for name in names:
        try:
            names_dict[name]
        except KeyError:
            names_dict[name] = name

    return names_dict



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Builds a tree from the output of the ANI/AF program. It is highly recommended to use the AF matrix for the tree.")
    parser.add_argument("-matr", help="the matrix to use", required=True)
    parser.add_argument("-alg", help="the tree building algorithm", choices=['nj', 'upgma'], default="nj")
    parser.add_argument("-out", help="the output file name for the tree", default="tree.tre")
    parser.add_argument("-format", help="the format to write the tree in", choices=["newick"], default="newick")
    parser.add_argument("-rename", help="file that has old_name<tab>new_name to rename a set of the nodes")
    args = parser.parse_args()


    nested_dict = parse_input_matrix(args.matr)

    if args.rename:
        names = list(nested_dict.keys())
        names_dict = create_names_dict(args.rename, names)
    else:
        names_dict = None

    distance_matrix = build_distance_matrix(nested_dict, names_dict)
    tree = build_tree(distance_matrix, alg=args.alg)
    Phylo.write(tree, args.out, args.format)
