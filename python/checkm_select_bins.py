
import argparse
import pandas

parser = argparse.ArgumentParser(description="Subsets a checkm tab-separated outfile to include only entries that have the specified completeness/contamination level")
parser.add_argument("-checkm", help="the checkm out file", required=True)
parser.add_argument("-completeness", help="completeness value")
parser.add_argument("-comp_metric", help="the comparison to completeness to select [%(default)s]", choices=["<", "<=", "=", ">=", ">"], default=">=")
parser.add_argument("-contamination", help="contamination value")
parser.add_argument("-cont_metric", help="the comparison to contamination to select [%(default)s]", choices=["<", "<=", "=", ">=", ">"], default="<=")
parser.add_argument("-out", help="the output checkm tsv output [%(default)s]", default="checkm_subset.tsv")

args = parser.parse_args()

df = pandas.read_csv(args.checkm, sep="\t", header=0, index_col=0)

if args.completeness:
    df = df.query("Completeness {metric} {value}".format(metric=args.comp_metric, value=args.completeness))

if args.contamination:
    df = df.query("Contamination {metric} {value}".format(metric=args.cont_metric, value=args.contamination))

df = df.sort_index(0, 'Completeness', ascending=False)
df.to_csv(args.out, sep="\t", index_label="Bin Id")

