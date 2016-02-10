#!/usr/bin/env python

description = """
Author: Hunter Cameron
Date: Aug 17, 2015
Description:

Small utility for manipulating a phylogenetic tree. It is possible to
keep (or remove) nodes based on a list supplied as a file or using a 
regular expression to match tree node names. There is also the option
to rename nodes on the list by providing a second column that is the 
desired new name.

As a convenience, it is also possible to list all the external nodes in
the tree. This can be used to generate a list to use to supply new names
or to check the results of prunings.

"""


import argparse
import re

from Bio import Phylo
from Bio.Phylo.PhyloXML import Taxonomy, Property


class TreeEditor(object):
    """ Wraps around BioPython to edit a tree """

    def __init__(self):
        self.tree = None

    def read_tree(self, tree_f, format="phyloxml"):
        """ Use biopython to read in a tree """

        # comments are confidence allows all digit node names to parse correctly
        # I think this is only needed for Newick
        if format == "newick":
            self.tree = Phylo.read(tree_f, format, comments_are_confidence=True)
        else:
            self.tree = Phylo.read(tree_f, format)
        
    def write_tree(self, format, path=None):
        """ Use biopython to write the tree """
        if not path:
            path = "edited_tree.tre"
        Phylo.write(self.tree, path, format)
  
    def list_nodes(self):
        """ Lists all of the external nodes in the tree """
        for node in self.tree.get_terminals():
            print(node.name)

    def edit(self, nodes, regex, remove):
        """ Just a wrapper for the recursive trim_tree function """

        self.trim_tree(self.tree.root, nodes, regex, remove)

    @classmethod
    def trim_tree(cls, root, nodes, regex, remove):
        """ Trims unwanted nodes and collapses unecessary internal nodes """
        # store a list of nodes that need to be pruned to avoid pruning while traversing
        to_prune = []

        for child in root.clades:
            '''
            This flag is here to detect internal nodes that have all their children removed
            and become terminal themselves. We want to delete such nodes.
            '''
            terminal_flag = child.is_terminal()

            # process the child node
            cls.trim_tree(child, nodes, regex, remove)
            
            # check if the child should be pruned
            if child.is_terminal():

                # check if child was terminal (before it was processed)
                if not terminal_flag:
                    to_prune.append(child)
                    continue
                
                # look for the child in the node list
                if child.name in nodes:
                    if remove:
                        to_prune.append(child)
                    else:
                        if nodes[child.name]:
                            child.name = nodes[child.name]

                # next check the name against the regex, if supplied
                elif regex:
                    if re.match(regex, child.name):
                        if remove:
                            to_prune.append(child)

                # node wasn't found in list, see if we want to remove it
                else:
                    if not remove:
                        to_prune.append(child)

        
        '''
        Pruning must be done separately because we can't change the list as we iterate over it
        while pruning, we want to make sure any nodes in the node list were found
        we also want to count how many nodes were removed
        '''
        for node in to_prune:
            root.clades.remove(node)


        # after the pruning, if there is only one child, there is no point in having the internal node
        # must make a copy because root.clades is modified with each call to collapse
        for child in root.clades[:]:
            if len(child.clades) == 1:
                root.collapse(child)


def parse_nodes(nodes_f):
    """ Parses nodes and (optionally) new names from a file; returns dict """
    nodes = {}
    with open(nodes_f, 'r') as IN:
        for line in IN:
            sp_line = line[:-1].split("\t")

            # store the optional rename if provided, none otherwise
            if len(sp_line) == 2:
                nodes[sp_line[0]] = sp_line[1]
            elif len(sp_line) == 1:
                nodes[sp_line[0]] = None
            else:
                raise ValueError("There were more columns than expected in line:\n  {}".format(line))

    return nodes

def main(args):

    # make the Editor object and read in the tree
    te = TreeEditor()
    te.read_tree(args.tree, args.format)

    # check if all we want to do is list nodes
    if args.list:
        te.list_nodes()
        return
    
    # process a nodes file if given
    if args.nodes:
        nodes = parse_nodes(args.nodes)
    else:
        nodes = {}

    # prune the tree
    te.edit(nodes, args.regex, args.remove)

    # write the tree
    if args.out_format:
        format = args.out_format
    else:
        format = args.format
    te.write_tree(format=format, path=args.out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-tree", help="the tree", required=True)
    parser.add_argument("-format", help="the format of the tree", choices=["newick", "nexus", "phyloxml"], required=True)
    parser.add_argument("-nodes", help="a list of nodes to keep, optional new names in the second column")
    parser.add_argument("-regex", help="a regular expression to match node names from the beginning")
    parser.add_argument("-remove", help="flag to indicate you want to remove the nodes in nodes/regex", action="store_true")
    parser.add_argument("-out", help="the path for the pruned tree")
    parser.add_argument("-out_format", help="the output format if it is different from the input", choices=["newick", "nexus", "phyloxml"])

    parser.add_argument("-list", help="writes all external node names to stdout; overwrites other options", action="store_true")

    args = parser.parse_args()

    main(args)
