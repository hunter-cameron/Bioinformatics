

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

def get_parents_for_ids(nodes_f):
    parents_dict = {}
    with open(nodes_f, 'r') as IN:
        for line in IN:
            elements = line.split("\t|\t")
            elements[-1] = elements[-1].rstrip("\t|\n")

            id=elements[0]
            parent=elements[1]
            rank=elements[2]

            if not parent:
                print("Warning: taxid {} doesn't have parent".format(id))

            parents_dict[id] = {'parent': parent, 'rank': rank}
    return parents_dict


def build_tree(name_dict, parents_dict):
    
    tree = taxtree.TaxTree()
    taxtree.TaxNode.set_default_tree(tree)
    
    # add the root node
    taxtree.TaxNode(taxid="131567", name="cellular organisms", rank="root", parent=None)

    # holds already used taxids to prevent repeats (and also only get cellular organisms)
    previous = {"1": 1, "131567": 1}
    prev_tier = ["131567"]
    cont = True
    while cont:
        print("\n" + "-" * 10)
        print(prev_tier)
        curr_tier = []
        num_processed = 0
        for parent_id in prev_tier:
            p_node = tree.lookup_taxid(parent_id)
            for taxid, data in parents_dict.items():
                if data['parent'] == parent_id:
                    if taxid not in previous:
                        curr_tier.append(taxid)
                
                    name = name_dict.get(taxid, "")
                    parent = p_node
                    rank = data.get('rank', "")


                    print((taxid, name, rank, data['parent']))
                    c_node = taxtree.TaxNode(taxid, name, rank, parent)
                    p_node.add_child(c_node)

                    num_processed += 1
                    if num_processed / 1000 == int(num_processed / 1000):
                        print(("Num processed", num_processed))
        
        # check for loop end
        if curr_tier:
            for tid in curr_tier:
                previous[tid] = 1
            prev_tier = curr_tier
        else:
            break

    return tree



if __name__ == "__main__":
    
    names = get_names_for_ids(sys.argv[1])

    parents = get_parents_for_ids(sys.argv[2])

    tree = build_tree(names, parents)
    
    print(len(names))
    print(len(parents))
    print(len(tree.taxnodes))


    tree.save_tree("tree.pickle")
