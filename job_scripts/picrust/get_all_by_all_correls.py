
import argparse
import pandas
from mypyli import picrust


def all_by_all_correlation(table_f, out):

    # calculate using lower triangular algorithm
    ttm = picrust.TraitTableManager(table_f)
    correlations = {}
    names = []
    for ent1 in ttm:
        names.append(ent1.name)
        for ent2 in ttm:
            if ent1 == ent2:
                continue

            cor = ent1.correlation(ent2, traits=ttm.traits)
            correlations[ent1.name + ">" + ent2.name] = str(cor)

    # write as a square
    with open(out, 'w') as OUT:
        OUT.write("\t".join(["names"] + names) + "\n")

        for n1 in names:
            to_write = [n1]
            for n2 in names:
                try:
                    to_write.append(correlations[n1 + ">" + n2])
                except KeyError:
                    to_write.append(correlations[n2 + ">" + n1])


            OUT.write("\t".join(to_write)+ "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="calculates correlations between entries in PICRUSt trait tables")
    parser.add_argument("-table", help="PICRUSt trait table", required=True)
    parser.add_argument("-out", help="output path for the correlations")

    args = parser.parse_args()

    all_by_all_correlation(args.table, args.out)
