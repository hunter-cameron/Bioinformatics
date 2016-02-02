
import argparse
import pandas
from mypyli import plotmaster


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="plots the average coverage over all contigs for each sample from the depth file metaBAT uses")
    parser.add_argument("-depth", help="the depth file", required=True)
    parser.add_argument("-out", help="the output file name [%(default)s]", default="coverage_barplot.png")

    args = parser.parse_args()

    # read in the depth file
    df = pandas.read_csv(args.depth, sep="\t", index_col=0, header=0)

    # get just the data columns
    all_cols = df.columns
    data_cols = [col for col in df.columns if col.endswith(".bam")]
    data_cols = sorted(data_cols)

    df = df[data_cols]

    avg = df.mean(axis=0)

    ax = avg.plot(kind="bar")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Average Coverage")
    ax.set_title("Average Contig Coverage Per Sample")
    plotmaster.save_ax(ax, args.out)
