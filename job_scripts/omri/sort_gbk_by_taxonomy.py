
import argparse
import os
from mypyli import genbank

def get_taxonomy(gbk_f):
    return genbank.GenBank(gbk_f).taxonomy


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sorts a bunch of GenBank files by the taxonomy string.")
    parser.add_argument("-gbks", help="a bunch of gbk files", nargs="+")
    parser.add_argument("-level", help="numerical level 0=kingdom, 1=phylum, etc", type=int, required=True)
    args = parser.parse_args()
    

    tax_groups = {}
    for gbk in args.gbks:
        name = os.path.splitext(os.path.basename(gbk))[0]
        tax = get_taxonomy(gbk)

        try:
            tax_groups[tax[args.level]].append(name)
        except:
            tax_groups[tax[args.level]] = [name]


    for group in sorted(tax_groups):
        print(group)
        
        for genome in tax_groups[group]:
            print(genome)

        print()
