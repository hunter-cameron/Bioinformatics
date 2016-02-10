
from __future__ import print_function
import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os

def read_data_table(data_f):
    """
    Converts raw CheckM data into a data frame.
    """
    # return a df with the first column as an index and the first row as headers
    return pd.read_csv(data_f, sep="\t", header=0, index_col=0)
    
def plot_comp_and_contam(frame, name):
    """ Makes a overlaid bar graph with completeness and contamination """
    x_vals = range(frame.shape[0])
    complete_vals = np.array(frame[["Completeness"]])
    contam_vals = np.array(frame[["Contamination"]])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("Completeness and Contamination Using CheckM")
    ax.bar(x_vals, complete_vals, 1, color='#deb0b0', align='center', label="completeness")

    #ax2 = ax.twinx()
    ax.bar(x_vals, contam_vals, 1, color='#b0c4de', align='center', label="contamination")

    names = list(frame.index)

    plt.xlim(x_vals[0] - 1, x_vals[-1] + 1)
    plt.xticks(x_vals, np.array(names), rotation="vertical", ha='center', fontsize='x-small')
    plt.xlabel("Genome")

    plt.ylabel("Percent")
    ax.yaxis.set_ticks_position("right")

    lgd = plt.legend(bbox_to_anchor=(1.1, 1), loc=2, borderaxespad=0)

    #plt.tight_layout()
    # save figure taking the legend into account to resize
    fig.savefig(name, bbox_extra_artists=(lgd,), bbox_inches='tight')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Make Completeness and Contamination bar plot from CheckM output file")
    parser.add_argument("-checkm", help="checkm tab output format", required=True)
    parser.add_argument("-base", help="base name for plots", required=True)
    parser.add_argument("-trim_comp", help="remove genomes below the specified completeness", type=float)
    parser.add_argument("-trim_cont", help="remove genomes above the specified contamination", type=float)
    parser.add_argument("-sort", help="sort genomes by completeness", action='store_true')
    args = parser.parse_args()

    checkm_table = read_data_table(args.checkm)

    if args.trim_comp:
        checkm_table = checkm_table.query('Completeness >= {}'.format(args.trim_comp))

    if args.trim_cont:
        checkm_table = checkm_table.query('Contamination <= {}'.format(args.trim_cont))

    if args.sort:
        checkm_table.sort_index(by='Completeness', inplace=True)

    plot_comp_and_contam(checkm_table, "{}.comp_and_cont.png".format(args.base))
