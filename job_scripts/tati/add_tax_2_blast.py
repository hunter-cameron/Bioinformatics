

from mypyli import taxtree
import sys

tree = taxtree.TaxTree.load_tree("/nas02/home/h/j/hjcamero/scripts/python/mypyli/taxtree_kingdom.pickle")

with open(sys.argv[1], 'r') as IN, open("blast_plus_tax.txt", 'w') as OUT:
    for line in IN:
        fields = line[:-1].split("\t")

        try:
            tax = tree.lookup_taxid(fields[-1].split(";")[-1])
            fields.append(tax.get_tax_string())
        except:
            fields.append("automatic lookup failed")

        OUT.write("\t".join(fields) + "\n")
