
# AUTHOR: Hunter Cameron
# DATE: 1/8/2015
# DESCRIPTION: Makes a table with kraken's assignment, confidence score, and contig length
# UPDATED:
#



from mypyli.kraken import KrakenIO
from mypyli import taxtree
import argparse 

def create_table(kraken_f, out_f, tree):
    KrakenIO.set_tree(tree)

    with open(kraken_f, 'r') as IN, open(out_f, 'w') as OUT:

        OUT.write("\t".join(["contig", "length", "taxonomy", "confidence"]) + "\n")
        for record in KrakenIO.parse(IN):
            if record.classified:
                tax = tree.lookup_taxid(record.taxid).get_tax_string()
                conf = record.get_kraken_confidence()
            else:
                tax = "unclassified"
                conf = "0/0"

            OUT.write("\t".join([record.name, str(record.length), tax, str(conf)]) + "\n")
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Makes a table with kraken's assignment, confidence score and contig length")
    parser.add_argument("-k", "--kraken", help="kraken output file", required=True)
    parser.add_argument("-o", "--out", help="name for output table", default="kraken_conf_table.txt")
    parser.add_argument("-t", "--taxtree", help="taxtree pickle file to load", required=True)

    args = parser.parse_args()

    tree = taxtree.TaxTree.load_tree(args.taxtree)
    create_table(args.kraken, args.out, tree)
