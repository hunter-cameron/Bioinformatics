
import argparse

from Bio import Phylo
from Bio.Phylo.PhyloXML import Taxonomy, Property

class Genome(object):
    """ A genome in the CheckM tree """

    def __init__(self, name, type, marker_genes=None, taxonomy=None):
        self.name = name

        self.type = type
        self._marker_genes = marker_genes
        self.taxonomy = taxonomy

    def __str__(self):
        string = self.name
        
        if self.marker_genes is not None:
            string += "\tmarkers: {}".format(self.marker_genes)

        if self.taxonomy is not None:
            string += "\t{}".format(self.taxonomy)

        return string

    @property
    def marker_genes(self):
        return self._marker_genes 
    @marker_genes.setter
    def marker_genes(self, value):
        self._marker_genes = int(value)

def parse_tree_qa(tree_qa_f):
    """ Parses the tree_qa file and yields Genome objects """

    with open(tree_qa_f, 'r') as IN:
        IN.readline()

        for line in IN:
            genome_name, num_markers, dup_markers, taxonomy = line[:-1].split("\t")
            
            yield Genome(genome_name, "qry", num_markers, taxonomy)
           
def parse_checkm_ref(checkm_ref_f):

    with open(checkm_ref_f, 'r') as IN:
        for line in IN:
            name, taxonomy = line[:-1].split("\t")

            yield Genome(name, "ref", taxonomy=taxonomy)

class TreeEditor(object):
    """ Wraps around BioPython to edit a tree """

    def __init__(self):
        self.tree = None

    def read_tree(self, tree_f, tre_format="phyloxml"):
        """ Use biopython to read in a tree """
        self.tree = Phylo.read(tree_f, tre_format)
        
    def write_tree(self, path="annotated_tree.phyloxml"):
        Phylo.write(self.tree, path, "phyloxml")
  
    def annotate_taxonomy(self, genome_dict, trim=False):
        """ 
        Adds taxonomy to the nodes in the keys of tax_dict 
        
        Do I only want to annotate terminal nodes or do I want to annotate them all?
        I think I probably only want to anntate the terminal nodes...
        """

        # trim first to remove unnecessary nodes
        if trim:
            self.trim_tree(self.tree.root, genome_dict)

        # loop over all remaining clades to update name info
        for clade in self.tree.find_clades():

            name = clade.name

            # format of UID names are: UIDxx|tax(may be blank)|bootstrap

            # lookup the genome object for the name
            try:
                genome = genome_dict[name]
            except KeyError:
                # genome was not in the dict; do something
                continue

            clade.name += " | "  + genome_dict[name].taxonomy


            """ ALL this is worthless because you can't search properties in A """
            # make a new property with the genome type
            # genome_type = Property(genome.type, "Genome:type", "node", "xsd:string")

            # add the property... might be a better way
            # clade.properties.append(genome_type)



            """ ALL this is worthless because I can't convert those names into taxid because they aren't formatted right """
            # make the Taxonomy object
            #taxonomy = Taxonomy(id=ID(genome.taxid, provider="ncbi_taxonomy"), rank=genome.tax_rank)

            # set the taxonomy
            #clade._set_taxonomy(taxonomy)

            # set the color of the branch (and all descendents unless overridden) 
            #clade._set_color()


    @classmethod
    def trim_tree(cls, root, genome_dict):
        """ Trims unwanted nodes and collapses unecessary internal nodes """
        print("Processing Node {}...".format(root))
        to_prune = []
        for child in root.clades:
            cls.trim_tree(child, genome_dict)

            if child.is_terminal():
                if not child.name in genome_dict:
                    to_prune.append(child)
        
        # pruning must be done separately because we can't change the list as we iterate over it
        for node in to_prune:
            root.clades.remove(node)


        # after the pruning, if there is only one child, there is no point in having the internal node
        # must make a copy because root.clades is modified with each call to collapse
        for child in root.clades[:]:
            if len(child.clades) == 1:
                root.collapse(child)


def test(args):
    genome_dict = {genome.name: genome for genome in parse_tree_qa(args.tree_qa) if genome.marker_genes >= 10}

    #genome_dict = {name: Genome(name, "Isolate", None, "Bacteria") for name in ["psg_MFr1_23", "iso_2561511073", "iso_2574179747"]}

    #genome_dict.update({genome.name: genome for genome in parse_checkm_ref(args.checkm_ref)})

    te = TreeEditor()
    te.read_tree(args.tree)
    te.annotate_taxonomy(genome_dict, trim=True)
    te.write_tree()

    return te

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-tree", help="the tree from storage/tree/concatenated.tre", required=True)
    parser.add_argument("-tree_qa", help="output from checkm tree_qa --tab_table")
    parser.add_argument("-checkm_ref", help="genome_tree.taxonomy.csv file from CheckM data directory")

    args = parser.parse_args()

    te = test(args)
