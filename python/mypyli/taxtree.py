

import sys
import argparse
import re
import difflib
from Bio import Entrez
import pickle


# TODO: Incorporate taxstrings lookup functions to allow the lookup of unknown data.



class TaxTree(object):
    """ A class to manage taxonomic relationships by acting as a container for TaxNode objects """
    

    @classmethod
    def load_tree(cls, filename):
        """ Loads and returns a pickled TaxTree """
        with open(filename, 'rb') as IN:
            tree = pickle.load(IN)

        return tree


    def __init__(self):
        self.taxnodes = {}

    def save_tree(self, filename):
        """ Saves a TaxTree using pickle. Can be loaded with TaxTree.load_tree """
        with open(filename, 'wb') as OUT:
            pickle.dump(self, OUT)

    def add_node(self, taxnode):
        """ Adds a TaxNode object to the TaxTree, checks if TaxNode already exists and prints error if so """
        if self.taxnodes.get(taxnode.taxid, None):
            print("Node {} already exists.".format(taxnode.taxid), file=sys.stderr)
        else:
            self.taxnodes[taxnode.taxid] = taxnode

    def lookup_taxid(self, taxid):
        """ Looks up a TaxNode by it's id. Raises a KeyError if taxid is not in Tree """
        try:
            return self.taxnodes[taxid]
        except KeyError as e:
            print("No entry for taxid: {}".format(taxid))
            raise e


class TaxNode(object):
    """ An entry in the TaxTree """

    taxTree = None
    taxRanks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxRanksEx = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    
    
    
    @classmethod
    def set_default_tree(cls, tree):
        """ Sets the tree for all TaxNodes """
        cls.taxTree = tree

    def __init__(self, taxid, name, rank, parent, children=[], taxtree=""):

        # set the instance tree to the class TaxTree unless other is specified
        if taxtree:
            self.taxtree = taxtree
        else:
            self.taxtree = self.taxTree

        if rank != "root":

            if not all([taxid, name, rank, parent]):
                print("\n".join(["{}= {}".format(name, var) for name, var in zip(["taxid", "name", "rank", "parent"], [taxid, name, rank, parent])]))
                raise ValueError("Cannot create tax node: taxid, name, rank, and parent required.")

        self.taxid = taxid
        self.name = name
        self.rank = rank
        self.parent = parent
        self.children = children

        self.add_node_to_tree()

    def __repr__(self):
        return "TaxNode with taxid: {}; name: {}; rank: {}".format(self.taxid, self.name, self.rank)

    def add_node_to_tree(self):
        """ Adds a node to TaxTree (error checking happens in TaxTree method) """
        self.taxtree.add_node(self)

    def add_child(self, child):
        """ Adds a new child to a node """

        self.children.append(child)
        #if child.taxid in [c.taxid for c in self.children]:
        #    print("Ignoring attempt to add an existing child to the Node.", file=sys.stderr)
        #else:
        #    self.children.append(child)

    def get_taxonomy(self):
        """ Returns an array of tuples with each tuple containing the rank and name of a step in the lineage. 
        I picked array over dict because several of the entries don't have a rank and would have used the same key.
        """
        # populate an array with all the lineage information
        tax_ranks = [(self.rank, self.name)]
        parent = self.parent
        while parent:

            tax_ranks = [(parent.rank, parent.name)] + tax_ranks

            parent = parent.parent

        return tax_ranks

    def get_tax_string(self, extended=False, trim_to=None, truncate=True, include_root=False):
        """ 
        Returns a tax string for the node. 

            Options:
                extended - include sub species (Ex: suborder == so_suborder)
                trim_to - make the specified rank the lowest one
                truncate - trim off unknown taxonomic levels
                include_root - include the root (Ex: root; k_Bacteria...)

            For example, to display a tax string at the genus level do: truncate=False, trim_to="genus".
        """

        if extended:
            t_ranks = self.taxRanksEx
        else:
            t_ranks = self.taxRanks
  
        taxonomy = self.get_taxonomy()
        
        #
        ## Now, I need some way to convert from a list of tuples to a dict or something I can order
        ## I think I'll just convert it to a dict
        #

        tax_dict = {k: v for k, v in taxonomy}
        tax_lst = []
        for rank in t_ranks:

            # get the rank's abbreviation
            if rank == "subphylum":
                abbrev = "sp"
            elif rank == "subspecies":
                abbrev = "ss"
            elif rank == "superkingdom":
                abbrev = "k"
            elif rank in self.taxRanks:     # take the first letter if it is a standard rank
                abbrev = rank[0]
            else:
                raise ValueError("Unknown rank {} in lineage".format(rank))

            # get the name at the rank
            name = tax_dict.get(rank, "")
              
            if name:
                tax_lst.append(abbrev + "__" + name)
            else:
                # truncate here if unknown and truncate is on
                if truncate:
                    break
                else:
                    tax_lst.append(abbrev + "__")
           
            # check if it is to the trim
            if trim_to:
                if rank == trim_to:
                   break

        # check if need to include root
        if include_root:
            tax_lst = ["root"] + tax_lst

        return "; ".join(tax_lst)

    def is_parent_of(self, tax_node):
        """ Checks if node(self) is the parent of another node(tax_node) """
        parent = tax_node.parent

        while parent:
            if parent == self:
                return True
            parent = parent.parent
        return False

    def get_tax_at_rank(self, rank, null=None):
        """ 
        Used like dict.get(key, [default]). Returns the tax at the rank given.
        
        If the rank is not available for the instance (even if only because the rank is misspelled, ex: speces) the "null" value will be returned. If there is no null value, an AssertionError is raised.

        """

        for r, name in self.get_taxonomy():
            if r == rank:
                return name
        else:
            if null is None:
                raise AssertionError("No entry for rank: {} in taxonomy of {}".format(rank, self))
            else:
                return null


# These functions are used to build a TaxTree that can be imported later. 
# they use a self import so when the object is pickled, it gets pickled 
# with an absolute path (mypyli.taxtree.TaxTree) so when this module is 
# imported, it is not required to import * (just import mypyli.taxtree).
# There may be a better way to do this but this is the best I've found.
def build_tree_from_NCBI(names_f, nodes_f):
    """ Builds and returns a TaxTree built from the NCBI names.dmp and nodes.dmp files """
    
    from mypyli import taxtree

    
    data_dict = {}

    # get the names for each taxid
    with open(names_f, 'r') as IN:
        for line in IN:
            elements = line.split("\t|\t")
            # replace the line break of the last element
            elements[-1] = elements[-1].rstrip("\t|\n")

            if elements[-1] == "scientific name":
                data_dict[elements[0]] = {'name': elements[1]}
    
    # get parent for each taxid         
    with open(nodes_f, 'r') as IN:
        for line in IN:
            elements = line.split("\t|\t")
            elements[-1] = elements[-1].rstrip("\t|\n")

            id=elements[0]
            parent=elements[1]
            rank=elements[2]

            if not parent:
                print("Warning: taxid {} doesn't have parent".format(id))

            data_dict[id].update({'parent': parent, 'rank': rank})

    # build the TaxTree
    
    tree = taxtree.TaxTree()
    taxtree.TaxNode.set_default_tree(tree)
  
    # make a parent to id dict
    parent2child = {}
    for taxid, data in data_dict.items():
        parent2child[data['parent']] = parent2child.get(data['parent'], []) + [taxid]


    # add the root node
    taxtree.TaxNode(taxid="131567", name="cellular organisms", rank="root", parent=None)

    print("Adding nodes...", file=sys.stderr)
    tree = _add_nodes_recurs(["131567"], parent2child, data_dict, tree)

    return tree

def _add_nodes_recurs(to_add, parent2child, data_dict, tree):
    """ 
    Adds nodes to a TaxTree from a dict (used in the build_tree_from_NCBI method)
    
    Dictionary is expected to look like this  id_dict = {taxid: {'name': id_name; 'rank': id_rank, 'parent': parent_id}}

    Where each 'parent' in each entry is another entry.

    Initial call will look like _add_nodes_recurs([root_taxid], parent2child, id_dict, tree)
            
        Where parent2child is essentially the inverse of id_dict, so, {parent_id: [children]}
        and tree is a TaxTree with a single root node whose id is an entry in parent2child
        and to_add is a list containing only the taxid of the root node.

    """
    from mypyli import taxtree


    next_add = []
    for parent in to_add:
        pnode = tree.lookup_taxid(parent)
        for child in parent2child.get(parent, []):
            entry = data_dict.get(child, None)
            
            # check if superkingdom should be changed to kingdom
            # I prefer to do this because the Bacteria entry is listed as superkingdom
            # and I would like it to be listed as kingdom instead.
            # If you want it as it appears in the taxonomy, comment out this if statement
            if entry['rank'] == "superkingdom":
                if pnode.rank == "root":
                    entry['rank'] = "kingdom"

            cnode = taxtree.TaxNode(child, name=entry['name'], rank=entry['rank'], parent=pnode)

            pnode.add_child(cnode)
            next_add.append(child)

    print("Processed {} nodes.".format(len(tree.taxnodes)))
    if next_add:
        tree = _add_nodes_recurs(next_add, parent2child, data_dict, tree)

    return tree




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Module for working with taxonomies and tax strings. Access mosr functionality by importing this to your script. Run this as a standalone module to build a TaxTree object to load into your script.")
    parser.add_argument("-names", help="NCBI taxonomy names.dmp file", required=True)
    parser.add_argument("-nodes", help="NCBI taxonomy nodes.dmp file", required=True)
    parser.add_argument("-out", help="filename for pickled output file", default="TaxTree.pickle")
    args = parser.parse_args()

    tree = build_tree_from_NCBI(args.names, args.nodes)
    tree.save_tree(args.out)
