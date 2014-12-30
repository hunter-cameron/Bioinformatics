

import sys
from mypyli import taxstring

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

    def get_full_taxonomy(self, taxonomy=""):
        if not taxonomy:
            taxonomy = {}
        
        parent = self.get_parent()
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
        
        # check if this node has been created as a parent
        if self.get_node(id):
            old = self.get_node(id)
            if not old.get_parent():
                old.set_parent(parent)
            if not old.get_name():
                old.set_name(name)

            self.add_node(id, old)
        else:
            self.id = id
            if self.id == "1":
                self.name = "root"
                self.parent = None
                self.rank = "root"
            else:
                #print("id=" + id)
                self.name = name
                self.rank = rank
                self.set_parent(parent)
            self.add_node(id, self)

    def get_parent(self):
        return self.parent

    def set_parent(self, parent):
        if parent:
            if self.get_node(parent):
                self.parent = self.get_node(parent)
            else:
                self.parent = TaxNode(id=parent)
        else:
            self.parent = ""

    def get_name(self):
        return self.name

    def set_name(self, name):
        self.name = name

    def get_rank(self):
        return self.rank

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

if __name__ == "__main__":
    
    names = get_names_for_ids(sys.argv[1])

    print(len(names))

    nodes = build_tree(sys.argv[2], names)
    
    print(len(nodes))

    for k in nodes:
        if not names.get(k, ""):
            print("No entry for {}".format(k))

    stop = 0
    for node in sorted(nodes, key=lambda node: int(node)):
        tax = nodes[node].get_full_taxonomy()

        tax_string = ["root"]
        for rank in ["kingdom", "phylum", "class", "order", "family", "genus", "species"]:
            tax_string.append("{}__{}".format(rank[0], tax.get(rank, "")))

        print("; ".join(tax_string))

        if stop > 10:
            sys.exit()
        stop += 1
