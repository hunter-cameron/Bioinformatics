
import argparse


def main(args):
    significant_kos = {}
    for asaf_file in args.kos:
        lineage = asaf_file.split("_")[0]
        
        significant_kos[lineage] = []

        with open(asaf_file, 'r') as IN:
            # skip header
            IN.readline()

            for line in IN:
                fields = line.split("\t")

                # ko name is first field - '"KO:' and a closing quote
                ko_name = fields[0][4:-1]

                q_value = float(fields[-1])
                if q_value < .05:
                    significant_kos[lineage].append(ko_name)

    with open("significant_kos_by_lineage.tab", 'w') as OUT:
        OUT.write("{}\t{}\n".format("lineage", "ko_list"))
        
        for lineage in sorted(significant_kos):
            OUT.write("{}\t{}\n".format(lineage, ";".join(significant_kos[lineage])))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Picks KOs that are significant (q < .05) form Asaf's list of files")
    parser.add_argument("-kos", help="one or more of Asaf's siginificance files", nargs="+")

    args = parser.parse_args()
    main(args)


