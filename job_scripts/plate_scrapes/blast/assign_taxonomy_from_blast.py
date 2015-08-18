

from mypyli import taxtree
import sys
import os
import argparse
import pandas as pd

MIN_ID = 0
MIN_LEN = 0


def basename(path):
    path = os.path.splitext(path)[0]
    path = os.path.basename(path)
    return path

def parse(blast_f):
    """
    Simple blast file parser to yield a dict of values
    """
    taxid_counts = {}
    with open(blast_f, 'r') as IN:
        previous = ''
        for line in IN:
            blast_record = {key: value for key, value in zip(['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids'], line[:-1].split("\t"))}

            # convet all types to correct ones
            for key in ['pident', 'evalue', 'bitscore']:
                blast_record[key] = float(blast_record[key])
            for key in ['length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send']:
                blast_record[key] = int(blast_record[key])

            # check if it is a duplicate record, below min len or id
            if blast_record['qseqid'] == previous:
                continue
            elif blast_record['length'] < MIN_LEN:
                d_min_len += 1
                continue
            elif blast_record['pident'] < MIN_ID:
                d_min_id += 1
                continue

            previous = blast_record['qseqid']

            yield blast_record

def count_tax_ids(blast_f, method="length"):
    taxid_counts = {}
    for record in parse(blast_f):
        # this is getting the largest tid (which seems to be the more
        # specific one, usually)
        tid = record['staxids'].split(";")[-1]
        
        if method == "length":
            aln_length = int(record['length'])
            taxid_counts[tid] = taxid_counts.get(tid, 0) + aln_length
        elif method == "count":
            taxid_counts[tid] = taxid_counts.get(tid, 0) + 1
        else:
           raise ValueError("Method is invalid.")

    return taxid_counts

def assign_tax_ids(tree, iterable):
    assignments = {}
    for taxid in iterable:
        #assignments[taxid] = tree.lookup_taxid(taxid)
        try:
            assignments[taxid] = tree.lookup_taxid(taxid)
        except:
            print("Could not assign tax to id: {}".format(taxid), file=sys.stderr)
            continue

    return assignments
        
def get_best_taxonomy(taxstring_counts, class_perc=.50):
 

    # build tree of counts   
    tree_counts = {'count': 0, 'children': {}}
    for taxstring, count in taxstring_counts.items():
        levels = taxstring.split("; ")
        tree_counts['count'] += count
        current = tree_counts['children']
        for level in levels:
            # make the level if necessary
            current[level] = current.get(level, {'count': 0, 'children': {}})

            # add the count
            current[level]['count'] = current[level]['count'] + count

            # increment the current level
            current = current[level]['children']
      
    # traverse the tree to get the optimal branch
    total = tree_counts['count']
    current = tree_counts['children']
    taxonomy = []
    while current:
        #print(current, end="\n\n\n\n")
        tax = max(current, key=lambda key: current[key]['count'])
        
        if current[tax]['count'] / total >= class_perc:
            taxonomy.append(tax)
            current = current[tax]['children']
        else:
            return "; ".join(taxonomy)

    return "; ".join(taxonomy)


def make_tree(taxonomies):
    """ Makes a simple text tree of the taxonomies. """

    class SimpleTree(object):
        """ Simple inline tree class """

        def __init__(self, name, is_terminal=False):
            self.name = name
            self.is_terminal = is_terminal
            
            if self.is_terminal:
                self.count = 1
            else:
                self.count = 0
    
            self.children = {}

        def add_child(self, name):
            if name in self.children:
                return self.children[name]
            else:
                self.children[name] = SimpleTree(name)
                return self.children[name]

        def print_tree(self, fh=sys.stdout, level=0):
            """ Recursively print the tree collapsing non-terminal nodes. """
            if self.is_terminal:
                print("   "*level + self.name + "\t" + str(self.count), file=fh)

            else:
                print("   "*level + self.name, file=fh)

            for name in sorted(self.children):
                child = self.children[name]
                child.print_tree(fh, level+1)


        def print_collapsed_tree(self, fh=sys.stdout, level=0, collapsed=[]):
            """ Recursively print the tree collapsing non-terminal nodes. """
           
            next_collapsed = collapsed.copy()
        
            # we always want to print terminal nodes
            if self.is_terminal:
                col_text = "; ".join(collapsed + [self.name])
                print("   "*level + col_text + "\t" + str(self.count), file=fh)

            else:
                # for non terminal nodes, we only want to print nodes with > 1 child
                if len(self.children) > 1:
                    col_text = "; ".join(collapsed + [self.name])
                    print("   "*level + col_text, file=fh)

                else:
                    next_collapsed.append(self.name)

            for name in sorted(self.children):
                child = self.children[name]
                child.print_collapsed_tree(fh, level+1, next_collapsed)

    root = SimpleTree("root")

    current_tax = root
    for taxonomy in taxonomies:
        levels = taxonomy.split("; ")

        for level in taxonomy.split("; "):
            current_tax = current_tax.add_child(level)
            
        current_tax.is_terminal = True
        current_tax.count += 1

        current_tax = root

    root.print_tree(level=-1)

def main(blast_f, tree, checkm_f=None):
    classifications = {}
    for file in blast_f:

        taxid_counts = count_tax_ids(file)
    
        taxid_assign = assign_tax_ids(tree, taxid_counts)

        taxstr_counts = {taxid_assign[taxid].get_tax_string(): taxid_counts[taxid] for taxid in taxid_assign}

        taxonomy = get_best_taxonomy(taxstr_counts, class_perc=.80)

        # get the taxid of the assignment
        #node = tree.lookup_taxstring(taxonomy)
        #print(node.taxid)
        
        classifications[basename(file)] = taxonomy

    make_tree(classifications.values())

    for genome, tax in classifications.items():
        print("{}\t{}".format(genome, tax))

    if checkm_f:
        add_tax_assignments_to_checkm(classifications, checkm_f)

def add_tax_assignments_to_checkm(classifications, checkm_table_f):
    ser = pd.Series(classifications)
    df_blast = pd.DataFrame(ser, columns=['blast_tax'])
    
    df_checkm = pd.read_csv(checkm_table_f, header=0, index_col=0, sep="\t")

    merged = pd.merge(df_checkm, df_blast, left_index=True, right_index=True)

    merged.to_csv("checkm_with_blast_tax.txt", sep="\t", index_label="name")

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-blast", help="Blast files in output format 6", required=True, nargs="+")
    parser.add_argument("-tree", help="A TaxTree object", required=True)
    parser.add_argument("-checkm", help="a checkm table from my program")
    args = parser.parse_args()

    tree = taxtree.TaxTree.load_tree(args.tree)
    main(args.blast, tree, args.checkm)
