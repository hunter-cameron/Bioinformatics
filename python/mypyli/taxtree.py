
from __future__ import print_function   # needed this module to be py2 compatible

import sys
import argparse
import re
import difflib
import pickle
from Bio import Entrez
import re


class TaxTree(object):
    """ A class to manage taxonomic relationships by acting as a container for TaxNode objects """
   
    taxRanks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
    taxRanksEx = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
   
    def __init__(self, remote=False, email=None):
        self.taxnodes = {}

        
        if remote:
            if email is None:
                raise ValueError("When remote mode is enabled, you must supply an email (Eg. tree = TaxTree(remote=True, email=myemail@domain.com).")
            Entrez.email = email

        self.email = email
        self.remote = remote

    # TREE UTILTY METHODS
    @classmethod
    def load_tree(cls, filename):
        """ Loads and returns a pickled TaxTree """
        with open(filename, 'rb') as IN:
            tree = pickle.load(IN)
        
        Entrez.email = tree.email

        return tree

    def save_tree(self, filename):
        """ Saves a TaxTree using pickle. Can be loaded with TaxTree.load_tree """
        with open(filename, 'wb') as OUT:
            pickle.dump(self, OUT)

    def enable_remote(self, email=None):
        """ 
        Allows connection to NCBI taxonomy databases to look up information.
        Requires an email address.
        """

        if self.email is None:
            if email is None:
                raise ValueError("Email address must be provided to enable remote mode! (Eg. obj.enable_remote(myemail@domain.com)")
            else:
                self.email = email

        self.remote = True

    def disable_remote(self, remove_email=False):
        """ Disables remote lookup, optionally removes the email address associated with this tree"""

        if remove_email:
            self.email = None

        self.remote = False

    # TREE CONSTRUCTION METHODS

    def add_node(self, taxnode):
        """ Adds a TaxNode object to the TaxTree, checks if TaxNode already exists and prints error if so """
        if self.taxnodes.get(taxnode.taxid, None):
            print("Node {} already exists.".format(taxnode.taxid), file=sys.stderr)
        else:
            self.taxnodes[taxnode.taxid] = taxnode

    # TREE LOOKUP METHODS

    def lookup_taxid(self, taxid):
        """ Looks up a TaxNode by it's id. Raises a KeyError if taxid is not in Tree """
        try:
            return self.taxnodes[taxid]
        except KeyError as e:
            #print("No entry for taxid: {} in taxtree: {}".format(taxid, self))
            raise e


    def _taxstring2dict(self, taxstring):
        """
        Returns the dict form of a raw tax string. Not optimized for different types of input
        Expects a string of form: 'k_aaa; p_bbb; c_ccc'
        """
        
        delimiter = ';'     # hard-coded now to allow easy delimiter switch (as argument) later
        
        tax_levels = taxstring.split(delimiter)
        # remove blank entry at the end if the string ends with the delim.
        if not tax_levels[-1].strip():
            tax_levels.pop()

        # build the taxonomy dict
        taxonomy = {}
        for level in tax_levels:
            level = level.strip()
            match = re.match("([kpcofgs])_{1,2}(.*)$", level)
            if match:
                for rank in TaxTree.taxRanks:
                    if match.group(1) == rank[0]:
                        if rank in taxonomy:
                            # right now just warn and skip. May want to throw error
                            print("Warning! rank '{}' already assigned for {}".format(rank, taxstring), file=sys.stderr)
                            return None
                        else:
                            if match.group(2):
                                taxonomy[rank] = match.group(2)

            else:
                print("Warning! Rank tags (k__, s__, etc) were not found for {}. Assigning taxonomy in order.".format(taxstring), file=sys.stderr)
                rnks = TaxTree.taxRanks[:]
                for level in tax_levels:
                    try:
                        rank = rnks.pop(0)
                    except:
                        break
                    taxonomy[rank] = level

        return taxonomy


    def lookup_single_tax(self, tax):
        """ 
        I needed a way to look up a single tax term (with no rank info)
        so this function was born. It may be incorporated into one of the other
        lookup methods in the future.
        """
        matches = []
        for node in self.taxnodes.values():
            if node.name.lower() == tax.lower():
                matches.append(node)

        if len(matches) == 1:
            return matches[0]
        elif len(matches) < 1:
            raise LookupError("Tax: {} was not found in the tree.".format(tax))
        elif len(matches) > 1:
            raise LookupError("Tax: {} returned multiple matches in the tree.".format(tax))

    def lookup_taxstring(self, taxstring):
        """ 
        Attempts to lookup a TaxNode by a taxstring.
        Returns the taxnode or raises exception
        This is still in development. Needs more work.
        Algorithm is currently to just loop through all tax nodes to find one that matches.
        This is actually quite fast
        
        Really, only need to search for lowest node. If the lowest node isn't 
        found, search is meaningless. 
        """
        try:
            taxonomy = self._taxstring2dict(taxstring)
        except:
            raise 
        for rank in reversed(TaxTree.taxRanks):
            if rank in taxonomy:
                for node in self.taxnodes.values():
                    if node.name.lower() == taxonomy[rank].lower():
                        # match rank
                        if node.rank == rank:
                            return node
                        else:
                            continue
                else:
                    if self.remote:
                        try:
                            return self._remote_lookup(taxonomy)

                        except:
                            raise
                    else:
                        raise Exception("Taxstring {} not found in Tree".format(taxstring))

    def lookup_taxstring_fast(self, taxstring):

        try:
            tax_dict = self._taxstring2dict(taxstring)
        except:
            raise 
        print("here")
    
        # lookup the root node
        node = self.lookup_taxid("131567") 

        taxonomy = []
        flag = False
        for rank in TaxTree.taxRanks:
            try:
                tax = tax_dict[rank]
                if flag:
                    raise Exception("Malformed taxstring {}. Blank fields were followed by more data.".format(taxstring))
                else:
                    taxonomy.append((rank, tax.replace("_", " ")))
            except KeyError:
                flag = True
                

        print(("taxonomy", taxonomy))
        while taxonomy:
            rank, taxname = taxonomy[0]

            print(("Looking for", rank, taxname))
            for child in node.children:
                if child.name == taxname and child.rank == rank:
                    if len(taxonomy) == 1:
                        return child
                    else:
                        node = child
                        taxonomy = taxonomy[1:]
                        break
            else:
                raise Exception("Taxstring {} not found in Tree".format(taxstring))


 

    def _remote_lookup(self, taxonomy):
        """ Tries to look up a taxonomy on NCBI's servers. """

        if Entrez.email == self.email:
            Entrez.email = self.email

        for rank in reversed(TaxTree.taxRanks):
            if rank in taxonomy:
                name = taxonomy[rank]
            else:
                continue

            id_list = self._entrez_esearch(name)

            for taxid in id_list:
                record = self._entrez_efetch(taxid)


                if record["Rank"] == ["superkingdom"]:
                    rec_rank = "kingdom"
                else:
                    rec_rank = record["Rank"]
                
                if rec_rank == rank:
                    if record["ScientificName"] == name:
                        # nodes match, need to add a new node?
                        try:
                            return self.lookup_taxid(str(taxid))
                        except KeyError:
                            raise Exception("Taxid {} was not found in the tree. Method for adding new taxids to the tree is still under development.".format(taxid))
                            parent = None
                            # add in the entry and necessary parent nodes
                            ref_taxonomy = record["LineageEx"]
                            if ref_taxonomy == "superkingdom":
                                ref_taxonomy = "kingdom"
                            for rank in TaxTree.taxRanks:
                                for taxa in ref_taxonomy:
                                    ref_rank = taxa["Rank"]
                                    if ref_rank == rank:
                                        ref_id = taxa["TaxId"] 
                                        try:    # lookup the node and set it as the parent for the next
                                            parent_node = self.lookup_taxid(ref_taxonomy)
                                        except:
                                            # need to add ancestor nodes
                                            parent_id = record.get("ParentTaxId", None)
                                            rec_taxonomy = record.get("LineageEx", None)
            else:
                raise Exception("Found no results that matches {}: {}".format(rank, name))

                
                
        else:
            print(taxonomy)
            raise Exception("Something crazy happened!")

 

    def _entrez_esearch(self, term):
        """ Searches for a term and returns a taxid/list of taxids for it. """
        search = Entrez.esearch(term=term, db="taxonomy", retmode="XML")
        id_list = Entrez.read(search)["IdList"]

        # give warning if more than one id was found
        if len(id_list) > 1:
            print("Warning: search returned more than one result for {}.".format(term), file=sys.stderr)

        return id_list

    def _entrez_efetch(self, taxid):
        fetch = Entrez.efetch(id=taxid, db="taxonomy", retmode="XML")
        return Entrez.read(fetch)[0]

         
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

    def is_ancestor_of(self, tax_node):
        """ Checks if node(self) is the parent of another node(tax_node) at any level in the ancestry """
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
        # swap synonyms 
        if rank == "domain":
            rank = "kingdom"

        for r, name in self.get_taxonomy():
            if r == rank:
                return name
        else:
            if null is None:
                raise AssertionError("No entry for rank: {} in taxonomy of {}".format(rank, self))
            else:
                return null


class TaxLookup(object):
    """ Handles the Entrez interface/lookup of taxids/taxstrings """

    def __init__(self, email, tree):
        self.email = email
        self.tree = tree

        Entrez.email = email


    def taxstring2dict(taxstring):
        """
        Returns the dict form of a raw tax string. Not optimized for different types of input
        Expects a string of form: 'k_aaa; p_bbb; c_ccc'
        """
        
        delimiter = ';'     # hard-coded now to allow easy delimiter switch (as argument) later
        tax_ranks = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

        
        tax_levels = taxstring.split(delimiter)
        # remove a blank entry at the end if the string ends with the delim.
        if not tax_levels[-1].strip():
            tax_levels.pop()

        # build the taxonomy dict
        taxonomy = {}
        for level in tax_levels:
            level = level.strip()
            match = re.match("([kpcofgs])_{1,2}(.*)$", level)
            if match:
                for rank in tax_ranks:
                    if match.group(1) == rank[0]:
                        if rank in taxonomy:
                            # right now just warn and skip. May want to throw error
                            print("Warning! rank '{}' already assigned for {}".format(rank, taxstring), file=sys.stderr)
                            return None
                        else:
                            taxonomy[rank] = group(2)

            else:
                print("Warning! Rank tags (k__, s__, etc) were not found for {}. Assigning taxonomy in order".format(taxstring), file=sys.stderr)
                rnks = tax_ranks.copy()
                for level in tax_levels:
                    try:
                        rank = rnks.pop(0)
                    except:
                        break
                    taxonomy[rank] = level

        return taxonomy 

    def lookup_taxonomy(self, taxonomy):
        """
        Looks up a taxonomy, [(rank1, name1), (rank2, name2)], in the NCBI database 
        I want the whole list so I can find the best match in the event of multiples"""
        
            
        (rank, name) = taxonomy[-1]

        id_list = self._entrez_esearch(name)

        

        for taxid in id_list:
            record = _entrez_efetch(taxid)

            parent_id = record.get("ParentTaxId", None)
            rec_taxonomy = record.get("Lineage")
        else:

           pass 

    def _entrez_esearch(self, term):
        """ Searches for a term and returns a taxid/list of taxids for it. """
        search = Entrez.esearch(term=term, db="taxonomy", retmode="XML")
        id_list = Entrez.read(search)["IdList"]

        # give warning if more than one id was found
        if len(id_list) > 1:
            print("Warning: search returned more than one result for {}.".format(term), file=sys.stderr)

        return id_list

    def _entrez_efetch(self, taxid):
        fetch = Entrez.efetch(id=taxid, db="taxonomy", retmode="XML")
        return Entrez.read(fetch)[0]


# These functions are used to build a TaxTree that can be imported later. 
# they use a self import so when the object is pickled, it gets pickled 
# with an absolute path (mypyli.taxtree.TaxTree) so when this module is 
# imported, it is not required to import * (just import mypyli.taxtree).
# There may be a better way to do this but this is the best I've found.
def build_tree_from_NCBI(names_f, nodes_f, remote=False, email=None):
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
    
    tree = taxtree.TaxTree(remote, email)
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

    parser = argparse.ArgumentParser(description="""Module for working with taxonomies and tax strings. Access mosr functionality by importing this to your script. Run this as a standalone module to build a TaxTree object to load into your script.
            
            NOTE: A tree you generate in this manner should be considered personal if the -remote option is enabled because the email address provided will be wrapped in the resulting TaxTree object (and will be avialiable to anyone with access to the object)!!!""")
    parser.add_argument("-names", help="NCBI taxonomy names.dmp file", required=True)
    parser.add_argument("-nodes", help="NCBI taxonomy nodes.dmp file", required=True)
    parser.add_argument("-out", help="filename for pickled output file", default="TaxTree.pickle")
    parser.add_argument("-remote", help="allow connection to NCBI servers to lookup missing info", action="store_true")
    parser.add_argument("-email", help="email address to use for Entrez (required if -remote is given)")
    args = parser.parse_args()

    if args.remote:
        if not args.email:
            raise ValueError("When remote mode is enabled, you must supply an email")

        print("="*70 + "\nWarning! A TaxTree is being generated using the -remote option. This tree will use the email address '{}' when executing Entrez searches. Therefore, this tree should not be shared with others whom should not execute searches using the provided email address.\n".format(args.email) + "="*70)
    
    
    tree = build_tree_from_NCBI(args.names, args.nodes, args.remote, args.email)
    tree.save_tree(args.out)
