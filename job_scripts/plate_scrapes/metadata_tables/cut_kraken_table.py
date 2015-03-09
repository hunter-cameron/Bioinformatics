

from mypyli import taxstring

import argparse
import numpy as npy

def split_table(table, rank, bacteria_only, out):
    """
    2 potential methods:
        1. iterate over lines and dump taxonomy line by line (slow b/c requires calls to get taxonomy for each line)
        2. read it all in, then collapse with a single call to taxonomy, then print out (fast but requires memory)

    Currently, the second method is used.
    """


    # build a master dict of all the taxonomies and values at each position
    with open(table, 'r') as IN:
        header = IN.readline()
        head_elem = header[:-1].split("\t")
        master_dict = {k: [] for k in head_elem}
        for line in IN:
            ln_elem = line[:-1].split("\t")
            for head, ln in zip(head_elem, ln_elem):
                if head != "contig":
                    master_dict[head].append(int(ln))
                else:
                    master_dict[head].append(ln)
                   
    
    #for k, v in master_dict.items():
    #    print((k, len(v)))
    
    # collapse the dict at whatever level necessary
    new_dict = {}
    for taxon, arr in master_dict.items():
        if taxon == "contig":
            new_dict['contig'] = arr
            continue

        taxstr = taxstring.TaxString(tax=taxon)
        if bacteria_only:
            if taxstr.get_tax_at_rank(rank="kingdom") != "Bacteria":
                continue

        
        new_tax = taxstr.get_tax_string(trim_to=rank, truncate=False)
        if len(new_dict.get(new_tax, [])):
            #print(new_dict[new_tax])
            #print(arr)
            new_dict[new_tax] = npy.add(new_dict[new_tax], arr)
        else:
            new_dict[new_tax] = arr

    # print the new table
    with open(out, 'w') as OUT:
        contig_names = new_dict.pop('contig')
        
        sorted_taxa = sorted(new_dict)
        OUT.write("\t".join(["contig"] + sorted_taxa) + "\n")
        for indx in range(len(contig_names)):
            to_write = [contig_names[indx]]
            for taxon in sorted_taxa:
                to_write.append(str(new_dict[taxon][indx]))

            OUT.write("\t".join(to_write) + "\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remakes the master kraken taxonomy counts table at a given taxonomic level.")
    parser.add_argument("-mt", help="master table file", required=True)
    parser.add_argument("-rank", help="lowest rank to appear in the output table (anything lower will be merged)", required=True, choices=["kingdom", "phylum", "class", "order", "family", "genus", "species"])
    parser.add_argument("--bacteria_only", help="only keep entries from bacteria", action="store_true", default=False)
    parser.add_argument("-out", help="output path for new table, default = contig_taxonomy_table_at_[rank].txt")

    args = parser.parse_args()

    if args.out is None:
        out = "contig_taxonomy_table_at_{}.txt".format(args.rank)
    else:
        out = args.out

    split_table(table=args.mt, rank=args.rank, bacteria_only=args.bacteria_only, out=out)
