

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

    data_dict = {}
    with open(names_f, 'r') as IN:
        for line in IN:
            elements = line.split("\t|\t")
            # replace the line break of the last element
            elements[-1] = elements[-1].rstrip("\t|\n")

            if elements[-1] == "scientific name":
                data_dict[elements[0]] = {'name': elements[1]}
             

    return data_dict

def get_parents_for_ids(nodes_f, data_dict):
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


def add_nodes_recurs(to_add, parent2child, data_dict, tree):
    next_add = []
    for parent in to_add:
        pnode = tree.lookup_taxid(parent)
        for child in parent2child.get(parent, []):
            entry = data_dict.get(child, None)

            cnode = taxtree.TaxNode(child, name=entry['name'], rank=entry['rank'], parent=pnode)

            pnode.add_child(cnode)
            next_add.append(child)

    print("Processed {} nodes.".format(len(tree.taxnodes)))
    if next_add:
        tree = add_nodes_recurs(next_add, parent2child, data_dict, tree)

    return tree

def build_tree(data_dict):
    
    tree = taxtree.TaxTree()
    taxtree.TaxNode.set_default_tree(tree)
  
    # make a parent to id dict
    parent2child = {}
    for taxid, data in data_dict.items():
        parent2child[data['parent']] = parent2child.get(data['parent'], []) + [taxid]


    # add the root node
    taxtree.TaxNode(taxid="131567", name="cellular organisms", rank="root", parent=None)

    print("Adding nodes...", file=sys.stderr)
    tree = add_nodes_recurs(["131567"], parent2child, data_dict, tree)

    return tree



if __name__ == "__main__":
    
    data_dict = get_names_for_ids(sys.argv[1])

    get_parents_for_ids(sys.argv[2], data_dict)

    tree = build_tree(data_dict)
    
    tree.save_tree("tree.pickle")
