
import argparse
import os

def parse_ko_metadata(metadata_f):
    """ Parses the ko metadata file and returns a list of KOs """
    
    kos = []
    with open(metadata_f, 'r') as IN:
        for line in IN:
            # skip header
            if line.startswith("Trait\t"):
                continue

            kos.append(line.split("\t", 1)[0])

    return sorted(kos)


def count_KOs(ko_table_f, ko_list):
    """ Parses a JGI KO table and counts KOs on the list """ 
    ko_counts = {k: 0 for k in ko_list}
    ko_not_in_list = 0
    with open(ko_table_f, 'r') as IN:
        for line in IN:
            # skip header
            if line.startswith("KO ID"):
                continue

            elems = line.rstrip().split("\t")
            
            ko = elems[0]
            try:
                count = elems[-1]
            except IndexError:
                print((line,))
                sys.exit()
        
            # strip the "KO:" off the ko name
            ko = ko[3:]

            try:
                ko_counts[ko] = float(count)
            except KeyError:
                ko_not_in_list += 1

    print("{} KOs not in the ko_list.".format(ko_not_in_list))
    return ko_counts


def main(args):
    ko_list = parse_ko_metadata(args.ko_metadata)
   
    with open(args.out, 'w') as OUT:
        OUT.write("\t".join(["OTU_IDs"] + ko_list) + "\n")

        for ko_table in args.ko:
            name = os.path.splitext(os.path.basename(ko_table))[0]

            ko_counts = count_KOs(ko_table, ko_list)

            # count the number of 
            ko_not_counted = list(ko_counts.values()).count(0)

            print("{} KOs not found in genome.".format(ko_not_counted))
            print("{} KOs found in the genome.".format(len(list(ko_counts.values())) - ko_not_counted))
            print("")

            to_write = [name]
            [to_write.append(str(ko_counts[ko])) for ko in ko_list]
            OUT.write("\t".join(to_write) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts JGI's KO table into a PICRUSt trait table")
    parser.add_argument("-ko", help="one or more JGI KO tables", nargs="+")
    parser.add_argument("-ko_metadata", help="the KO metadata table from the PICRUSt deconstructed files", required=True)
    parser.add_argument("-out", help="output path for the trait table", default="trait_table.tab")

    args = parser.parse_args()

    main(args)
