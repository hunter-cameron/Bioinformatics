
import sys
from Bio import Entrez

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
                        record = cls._entrez_fetch(taxid)[0]
                        break
                else:
                    print("No suitable result found.")
                    record = {}     # setting record to this will act as if nothing was found


            else:
                print("  Warning! Could not identify rank, no suitable result found.")
                record = {}
        # if only one record, use that (still should do some checking here probably    
        else:
            record = cls._entrez_fetch(id_list[0])
            taxid = id_list[0]
        
        tax_dict = {}
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
            tax_dict[record["Rank"]] = record["ScientificName"]
                
        else:
            print("Entrez lookup failed for {}".format("'" + tax + "'"))
            return ({}, "Not assigned")

        return (tax_dict, taxid)

    @staticmethod
    def _entrez_search(term):
        search = Entrez.esearch(term=term, db="taxonomy", retmode="XML")
        return Entrez.read(search)["IdList"]

    @staticmethod
    def _entrez_fetch(tax_id):
        fetch = Entrez.efetch(id=tax_id, db="taxonomy", retmode="XML")
        return Entrez.read(fetch)[0]

    
    def __init__(self, gen_name, tax):
        (self.taxonomy, self.taxid) = self._tax_string_to_tax_dict(tax)
        self.name = gen_name

    def get_tax_string(self, include_root=False, truncate=True):
        tax_string = ""
        if include_root:
            tax_string += "root"
        
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            if not self.taxonomy.get(rank, 0):
                if truncate:
                    continue
            tax_string += "; {}__{}".format(rank[0], self.taxonomy.get(rank, ""))
        if not tax_string.startswith("root"):
            tax_string = tax_string[2:]
        return tax_string

class CheckmOut(object):
    """ Class to parse checkm output, create genome objects, and manage them. -- Mostly works as a container"""
    def __init__(self, checkm_ext):
        self.genomes = []
        with open(checkm_ext, 'r') as IN:
            data = eval(IN.read())
        for gen, info in data.items():
            info["name"] = gen
            self.genomes.append(CheckmGenome(info))

    def checkm_to_taxid(checkm_ext):
        pass


    def print_genome_table(self, fh=sys.stdout, lookup_tax=False):
        """ Prints a table of info for each genome. Import to R in mind """
        
        # write a header
        # fh.write("\t".join(["genome", "completeness", "contamination", "marker_lineage", "0", "1", "2", "3", "4", "5+"]) + "\n")

        for genome in self.genomes:
            if lookup_tax:
                taxstring = TaxString(genome.get_name(), genome.get_marker_lineage())
                new_lin = taxstring.get_tax_string()
                if new_lin:
                    genome.set_marker_lineage(new_lin.split("; ")[-1])

            fh.write("\t".join([genome.get_name(),
                                str(genome.get_completeness()),
                                str(genome.get_contamination()),
                                genome.get_marker_lineage(),
                                ] +
                                [str(genome.get_marker_counts()[key]) for key in ["0", "1", "2", "3", "4", "5+"]]
                                ) + "\n")

class CheckmGenome(object):
    def __init__(self, gen_hash):
        
        # essential parameters
        try:

            self.name = gen_hash["name"]
            self.completeness = gen_hash["Completeness"]
            self.contamination = gen_hash["Contamination"]
            self.marker_lineage = gen_hash["marker lineage"]
            self.marker_counts = {k: gen_hash[k] for k in ["0", "1", "2", "3", "4", "5+"]}

        except KeyError as err:
            print(err)
            print("CheckmGenome object not successfully created! Essential parameter missing.")
            sys.exit()
    
    def get_name(self):
        return self.name

    def get_completeness(self):
        return self.completeness

    def get_contamination(self):
        return self.contamination
    
    def get_marker_lineage(self):
        return self.marker_lineage

    def set_marker_lineage(self, new_lin):
        self.marker_lineage= new_lin

    def get_marker_counts(self):
        return self.marker_counts




# printing functions
def print_tax_id(genomes):
   
    tax_objects = []
    for (gen_name, tax) in genomes:
        print(gen_name + "\t" + tax)
        tax_objects.append(TaxString(gen_name=gen_name, tax=tax))


    print("\t".join(["Genome_name", "taxid", "taxstring"]))
    for tax in tax_objects:
        print("\t".join([tax.name, tax.taxid, tax.get_tax_string(include_root=True)]))

def print_tax_for_R(genomes):
    tax_objects = []
    for (gen_name, tax) in genomes:
        tax_string = tax.get_tax_string(include_root=True, truncate=False)

def print_info_for_bar_graph(checkm_ext):

    with open(checkm_ext, 'r') as IN:

        data = eval(IN.read())


if __name__ == "__main__":
    """
    Work done is to parse CheckM's extended output file for which marker bin was used and dump a list of tax ids
    """
    
    parser = CheckmOut(sys.argv[1])

    parser.print_genome_table(lookup_tax=True)
