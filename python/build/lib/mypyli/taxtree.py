

import sys
import argparse
import re
import difflib
from Bio import Entrez
import pickle



class TaxString(object):

    @classmethod
    def _tax_string_to_tax_dict(cls, tax, email='hjcamero@live.unc.edu'):
        
        # if tax is a full/partial string, get only lowest element
        if ";" in tax:
            tax = tax.split[";"][-1]
        
        # if tax string is in the format: k__Bacteria; p__abcdef, get rank to use later
        if "__" in tax:
            rank = tax[0]
            tax = tax[3:]
        else:
            rank = ""
       
        Entrez.email = email

        id_list = cls._entrez_search(tax)
        # try to resolve which id should be used if multiple results are returned
        # this section could use more work
        if len(id_list) > 1:
            print("Warning: search returned more than one result for {}.".format(tax))
            if rank:
                
                for tid in id_list:
                    record = cls._entrez_fetch(tid)
                    rec_rank = record["Rank"]
                    if rec_rank == "superkingdom":
                        abrev = "k"
                    else:
                        abrev = rec_rank[0]

                    if rank == abrev:
                        print("Found best result: taxid = {}".format(tid))
                        taxid = tid
                        #record = cls._entrez_fetch(taxid)[0]
                        break
                else:
                    print("No suitable result found.")
                    #record = {}     # setting record to this will act as if nothing was found
                    taxid = None

            else:
                print("  Warning! Could not identify rank, no suitable result found.")
                #record = {}
                taxid = None

        # if only one record, use that (still should do some checking here probably    
        else:
            #record = cls._entrez_fetch(id_list[0])
            taxid = id_list[0]
        
        if taxid:
            taxdict, parent_id = cls._taxid_to_taxdict(taxid, rank)
            return taxdict, taxid, parent_id
        else:
            return {}, taxid, None

    @classmethod
    def _taxid_to_taxdict(cls, taxid, rank="", email='hjcamero@live.unc.edu'):
       
        Entrez.email = email
        record = cls._entrez_fetch(taxid)
        #return(record)     # uncomment this to assign the record to the taxonomy
        tax_dict = {}

        parent_id = record.get("ParentTaxId", None)
        if record.get("LineageEx", 0):
            for entry in record["LineageEx"]:
                rank = entry["Rank"]
                if rank == "superkingdom":
                    rank = "kingdom"
                elif rank == "no rank":
                    continue

                name = entry["ScientificName"]
                # taxid = entry["TaxId"]        # may be used for looking up relatedness to build tree with accurate lengths?

                tax_dict[rank] = name

            # input entry for current record
            if record["Rank"] == "superkingdom":
                record["Rank"] = "kingdom"
            tax_dict[record["Rank"]] = record["ScientificName"]
            return tax_dict, parent_id
                
        else:
            print("Entrez lookup failed for {}".format("'" + taxid + "'"))
            return {}, None

    @staticmethod
    def _parse_input_string(tax):
        """ Parses tax strings and returns a dict with key=rank, value=tax 
        
        Since this doesn't know what kind of input it will get, this is the procedure:
            For each ; separated element - 
            1. Look for the rank prefix
            2. Assign to the highest rank available
        """
        elements = tax.split(";")
 
        #print(tax)

        tax_dict = {}
        for element in elements:
            element = element.strip()
           
            # check if root is included
            if element == "root":
                continue

            # check for the rank prefixes
            match = re.search("^([kpcofgs]_{1,2})", element)

            if match:
                # remove the prefix
                #print(match.groups(1))
                tax = element.replace(match.groups(1)[0], "")
                #print(tax, match.groups(1)[0])
                for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
                    if rank[0] == match.groups(1)[0][0]:
                        entry = tax_dict.get(rank, "")
                        
                        if entry:
                            raise AssertionError("Error: Rank '{}' has already been assigned for {}.".format(rank, tax))
                        else:
                            tax_dict[rank] = tax

            else:
                for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
                    entry = tax_dict.get(rank, "")
                    if entry:
                        pass
                    else:
                        tax_dict[rank] = element
                        break
                else:
                    raise AssertionError("All ranks have already been assigned. No place for '{}'.\n{}".format(element, tax))

        # make sure all elements have been assigned
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            entry = tax_dict.get(rank, "")
            if not entry:
                tax_dict[rank] = ""

        return tax_dict

    @staticmethod
    def _entrez_search(term):
        search = Entrez.esearch(term=term, db="taxonomy", retmode="XML")
        return Entrez.read(search)["IdList"]

    @staticmethod
    def _entrez_fetch(tax_id):
        fetch = Entrez.efetch(id=tax_id, db="taxonomy", retmode="XML")
        return Entrez.read(fetch)[0]

    @staticmethod
    def is_same_tax(tax1="", tax2="", rank=""):
        """
        Matches taxonomy at a given rank. Needs a better built-in fuzzy string matcher.
        
        Percent match is a bad idea because scientific names have conserved parts (Alphaproteobacteria, betaproteobacteria)
        that would make them appear to be better matches than they are. 
        Also, percent is dependent on string size 'bacteria' to 'bacteriia' is worse than a 90% match. (maybe)
        
        Ex: JGI calls 'sphingobacteria' => 'sphingobacteriia'
        """
        try:
            first = tax1.get_tax_at_rank(rank)
            second = tax2.get_tax_at_rank(rank)

        # Errors result from the tax not having the rank (or the rank being misspelled, etc)
        except AssertionError:
            return False

        if first == second:
            return True
        else:
            match = difflib.SequenceMatcher(None, first, second)
            
            return match.ratio() > .9
            
    @classmethod
    def lineage_divides_at(cls, tax1="", tax2=""):
        """ Returns the rank at which the lineage of two tax divide. """
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
        
            if not cls.is_same_tax(tax1, tax2, rank):
                return rank
        else:   # this will return "species" if the lineage never divides
            return "Never"

    def __init__(self, tax="", is_id=False, name="", lookup=False):

        if lookup:
            if is_id:
                self.taxid = tax
                self.taxonomy, self.parent = self._taxid_to_taxdict(tax)
            else:
                self.taxonomy, self.taxid, self.parent = self._tax_string_to_tax_dict(tax)
        else:
            self.taxonomy = self._parse_input_string(tax)
            self.taxid = None
            self.parent = None

        self.name = name

    def get_taxid(self):
        if self.taxid:
            return self.taxid
        else:
            print("Taxid not defined. Returning ''.")
            return ""

    def get_parent_id(self):
        if self.parent:
            return self.parent
        else:
            print("Parent id not defined. Returning ''.")
            return ""

    def get_taxonomy(self):
        if self.taxonomy:
            return self.taxonomy
        else:
            raise AssertionError("Taxonomy not defined!")

    def get_tax_string(self, trim_to=None, include_root=False, truncate=True):
        tax_string = ""
        if include_root:
            tax_string += "root"
        
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            if not self.taxonomy.get(rank, 0):
                if truncate:
                    continue
            tax_string += "; {}__{}".format(rank[0], self.taxonomy.get(rank, ""))
            
            # if trim is specified, stop once reached
            if rank == trim_to:
                break

        if not tax_string.startswith("root"):
            tax_string = tax_string[2:]
        if tax_string:
            return tax_string
        else:
            print("Tax string is null (taxonomy was never assigned).", file=sys.stderr)
            print(self.taxonomy, file=sys.stderr)
            return self.get_taxid()

    def get_lowest_rank(self):
        for rank in ["species", "genus", "family", "order", "class", "phylum", "kingdom"]:
            try:
                self.get_tax_at_rank(rank=rank)
                return rank
            except AssertionError:
                pass

        return ""

    def get_tax_at_rank(self, rank="", suppress=False):
        rank = rank.lower()

        if rank == "domain":
            rank = "kingdom"

        if self.taxonomy.get(rank, 0):
            return self.taxonomy[rank]
        else:
            if suppress:
                print("Warning: No entry for rank: '{}' in entry: '{}'\nSuppressing error...".format(rank, self.name))
                return ""
            else:
                raise AssertionError("No entry for rank: '{}' in entry: '{}'".format(rank, self.name))

    def get_LCA(self, other):
        """ Returns a tuple (rank, name) for the lowest common ancestor between two """
        tax1 = self.get_taxonomy()
        tax2 = other.get_taxonomy()

        for rank in ["species", "genus", "family", "class", "order", "phylum", "kingdom"]:
            first = tax1.get(rank, None)
            second = tax2.get(rank, 0)

            if first == second:
                return rank, first
        else:
            return None, None

class TaxTree(object):
    """ A class to manage taxonomic relationships """


    @classmethod
    def load_tree(cls, filename):
        with open(filename, 'rb') as IN:
            tree = pickle.load(IN)

        return tree



    def __init__(self):
        self.taxnodes = {}

    def save_tree(self, filename):
        with open(filename, 'wb') as OUT:
            pickle.dump(self, OUT)

    def add_node(self, taxnode):
        if self.taxnodes.get(taxnode.taxid, None):
            print("Node {} already exists.".format(taxnode.taxid), file=sys.stderr)
        else:
            self.taxnodes[taxnode.taxid] = taxnode

    def lookup_taxid(self, taxid):
        try:
            return self.taxnodes[taxid]
        except KeyError as e:
            print("No entry for taxid: {}".format(taxid))
            raise e


    def build_tree_from_file(self, tree_f):
        """ Builds a tax tree from a flat file in the format taxid<tab>tax string """
        tree_dict = {}
        with open(tree_f, 'r') as IN:
            for line in IN:
                taxid, taxonomy = line[:-1].split("\t")
                print((taxid, taxonomy))
                
                # separate path for the root node
                if taxid == "1":
                    tax = 'root'
                    lowest_rank = 'root'
                    name = 'root'
                    parent = ""
                    parent_rank = ""

                else:
                    tax = TaxString(taxonomy)
                    lowest_rank = tax.get_lowest_rank()
                    name = tax.get_tax_at_rank(lowest_rank, suppress=True)

                    # get the parent name and rank
                    set_rnk = False
                    for indx, rank in enumerate(["species", "genus", "family", "order", "class", "phylum", "kingdom"]):
                        if set_rnk:
                            parent_rank = rank
                            parent_name = tax.get_tax_at_rank(parent_rank, suppress=True)
                            break
                    
                        if rank == lowest_rank:
                            set_rnk = True

                    else:
                        parent_rank = "root"
                        parent_name = "root"

                TmpBuild(taxid, name, lowest_rank, parent, parent_rank)

        
        TmpBuild.build()
    
        return TmpBuild

class TmpBuild(object):
    """ Temporary object to hold a data structure to help build the tree from file. """

    objDict = {}
    
    @classmethod
    def build(cls):
        """ Replaces the rank + name notation with the object """

        for name in cls.objDict:
            parent = cls.objDict[name].parent
            if parent == "root":
                continue
            parent_obj = cls.objDict.get(parent, parent)
            
            if type(parent_obj) == "TmpBuild":
                parent_obj.add_child(cls.objDict[name])
                cls.objDict[name] = parent_obj
            else:
                raise AssertionError("Could not find parent {}".format(parent))

            
            
    def get_existing(self):

        return self.objDict.get(self.rank + self.name, None)


    def store_obj(self):
        self.objDict[self.rank + self.name] = self

    def __init__(self, taxid, name, rank, parent, p_rank):
        
        self.name = name
        self.rank = rank
        
        existing = self.get_existing()
        
        if not existing:
            self.taxid = taxid
            self.parent = p_rank + parent
            self.children = {}
            self.store_obj()


    def add_child(self, child):
        if not child.rank + child.name in self.children:
            self.children[child.rank + child.name] = child

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

        # set the taxtree to the default
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
       

    def get_tax_at_rank(self, rank, null=None):
        for r, name in self.get_taxonomy():
            if r == rank:
                return name
        else:
            if null is None:
                raise AssertionError("No entry for rank: {} in taxonomy of {}".format(rank, self))
            else:
                return null

if __name__ == "__main__":

    tree = TaxTree()
    tmpbuild = tree.build_tree_from_file(sys.argv[1])
