
import pandas
import matplotlib.pyplot as plt
import sys
import argparse
from mypyli import plotmaster



if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument("-counts", help="read counts file from my script in khmer section of job scripts", required=True)
    parser.add_argument("-out", help="name for the output graphic", default="read_counts.png")
    args = parser.parse_args()

    frame = pandas.read_csv(args.counts, header=0, index_col=0, sep="\t")

    frame.columns = ['raw reads (paired)', 'trimmo singles', 'trimmo paired', 'khmer singles', 'khmer paired']
    frame['trimmo'] = frame['trimmo singles'] + frame['trimmo paired']
    frame['total khmer'] = frame['khmer singles'] + frame['khmer paired']

    frame.sort_index(inplace=True)

    ax = plotmaster.stacked_bar_plot(frame, ['raw reads (paired)', 'trimmo', 'total khmer', 'khmer paired'], log=False)
    ax.set_title("Plate Scrapes Read Counts")
    ax.set_ylabel("number of reads")
    ax.set_xlabel("sample")

    lgd = ax.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0)
    fig = ax.get_figure()
    fig.savefig(args.out, bbox_extra_artists=(lgd,), bbox_inches='tight')
