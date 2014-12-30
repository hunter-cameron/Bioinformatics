

import sys
from mypyli import taxstring, taxtree

class TaxNode(object):

    taxNodeCollection = {}

    @classmethod
    def get_node(cls, id):
        return cls.taxNodeCollection.get(id, None)

    @classmethod
    def add_node(cls, id, taxnode):
        cls.taxNodeCollection[id] = taxnode

    @classmethod
    def get_all_nodes(cls):
        return cls.taxNodeCollection

    @classmethod
    def create_node_from_id(cls, id):
        """ 
        Uses the taxstring module to create a full TaxNode from only a taxid.
        The lookup is slow so this is best used as a supplement to fill in gaps.
        """
        tax = taxstring.TaxString(id, is_id=True, lookup=True)
        rank = tax.get_lowest_rank()
        name = tax.get_tax_at_rank(rank)
        parent = tax.get_parent_id()

        cls(id, name, rank, parent)



    def get_full_taxonomy(self, taxonomy=""):
        if not taxonomy:
            taxonomy = {}
        if self.parent:
            parent = self.get_parent_node()
        else:
            parent = None
        rank = self.get_rank()
        
        if parent:
            taxonomy = parent.get_full_taxonomy(taxonomy)

        if not self.get_name():
            print("looking up {}".format(self.id))
            taxstr = taxstring.TaxString(tax=self.id, is_id=True, lookup=True)
            name = taxstr.get_tax_at_rank(taxstr.get_lowest_rank(), suppress=True)
            print(name)
            self.set_name(name)
        if rank:
            if rank == "superkingdom":
                rank = "kingdom"
            taxonomy[rank] = self.get_name()

        return taxonomy

    def __init__(self, id, name="", rank="", parent=""):
        
        # check if this node has already been created
        if self.get_node(id):
            print("Warning: attempting to add duplicate node. Ignoring attempt.")
        else:       # add a new node
            self.id = id
            if self.id == "1":
                self.name = "root"
                self.parent = None
                self.rank = "root"
            else:
                #print("id=" + id)
                self.name = name
                self.rank = rank
                self.parent = parent
            self.add_node(id, self)


    def get_parent_node(self):
        if self.parent:
            if self.get_node(self.parent):
                return self.get_node(self.parent)
            else:   # we have a taxid but no parent entry so make a new entry
                print("Looking up parent {}".format(self.parent))
                self.create_node_from_id(self.parent)
                return self.get_node(self.parent)
        else:
            raise AssertionError("Node has no parent.")

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def get_rank(self):
        return self.rank

    def get_full_lineage(self, lineage=[]):
        """ Returns a list of the nodes in the lineage (ending with the lowest) """
        
        if self.parent:
            parent = self.get_parent_node()
        else:
            parent = None

        if not self.get_name():
            print("looking up {}".format(self.id))
            taxstr = taxstring.TaxString(tax=self.id, is_id=True, lookup=True)
            name = taxstr.get_tax_at_rank(taxstr.get_lowest_rank(), suppress=True)
            print(name)
            self.set_name(name)

        lineage = [self] + lineage

        if parent:
            parents = parent.get_full_lineage(lineage)


        return lineage




def get_names_for_ids(names_f):

    name_dict = {}
    with open(names_f, 'r') as IN:
        for line in IN:
            elements = line.split("\t|\t")
            # replace the line break of the last element
            elements[-1] = elements[-1].rstrip("\t|\n")

            if elements[-1] == "scientific name":
                name_dict[elements[0]] = elements[1]
             

    return name_dict

def build_tree(nodes_f, name_dict):
    
    with open(nodes_f, 'r') as IN:
        for line in IN:
            elements = line.split("\t|\t")
            elements[-1] = elements[-1].rstrip("\t|\n")

            TaxNode(id=elements[0], parent=elements[1], rank=elements[2],
                    name=name_dict[elements[0]])
    

    return TaxNode.get_all_nodes()


def convert_to_taxtree(nodes):
        
    # init the tree and set all nodes as belonging to the tree
    tree = taxtree.TaxTree()
    taxtree.TaxNode.set_default_tree(tree)

    for taxid, node in nodes.items():
        # start at the root and work down such that each node has a parent (except root)
        parent = None
        for elem in node.get_full_lineage():
            if parent is None:
                p_node = None
            else:
                p_node = tree.lookup_taxid(parent)
            parent = taxtree.TaxNode(taxid=taxid, name=elem.get_name(), rank = elem.get_rank(), parent=p_node)
            parent = taxid
    return tree



if __name__ == "__main__":
    
    names = get_names_for_ids(sys.argv[1])

    print(len(names))

    nodes = build_tree(sys.argv[2], names)
    
    print(len(nodes))

    for k in nodes:
        if not names.get(k, ""):
            print("No entry for {}".format(k))

    
    tree = convert_to_taxtree(nodes)


    stop = 0
    with open("tax_dump.txt", 'w') as OUT:
        for node in sorted(nodes, key=lambda node: int(node)):
            tax = nodes[node].get_full_taxonomy()

            tax_string = ["root"]
            for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
                tax_string.append("{}__{}".format(rank[0], tax.get(rank, "")))

            OUT.write(node + "\t" + "; ".join(tax_string) + "\n")
    
            #if stop > 10:
            #    sys.exit()
            #stop += 1
