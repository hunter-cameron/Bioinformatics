#!/usr/bin/env python3

"""
Author: Hunter Cameron
Date: 11/31/2014
Updated: 11/31/2014
Description:
    A module for working with taxonomy strings.

    NOTE: A single taxonomy can have multiple taxon ids ex: '1747','1134454' are the exact same species, just different strains.
"""


import re
import difflib
from Bio import Entrez
import sys

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
                            raise AssertionError("Error: Rank '{}' has already been assigned".format(rank))
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
