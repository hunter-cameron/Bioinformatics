
from __future__ import print_function
import argparse
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from mypyli import isolatedb
import os

def read_data_table(data_f):
    """
    Converts raw CheckM data into a data frame.
    Should have an option to accept a .csv file from CheckM
    """
    # return a df with the first column as an index and the first row as headers
    return pd.read_csv(data_f, sep="\t", header=0, index_col=0)
    

def make_fraction_matrix(df1, df2, columns):
    """
    Returns a matrix that contains the fraction the second data frame is of
    the first for each specified column
    """
    return df2[columns].divide(df1[columns], axis=0, fill_value=0)

def main2(ref_f, queries_f):
    ref = read_data_table(ref_f)

    queries = []
    for query_f in queries_f:
        queries.append(read_data_table(query_f))

    #print(ref[['completeness', 'contamination']])
    #print(queries[0][['completeness', 'contamination']])
    #sys.exit()
    for query, file in zip(queries, queries_f):

        frac_matr = make_fraction_matrix(ref, query, ['completeness', 'contamination', 'genome_size'])
        print(frac_matr)

        print("{} Summary:".format(file))
        print(frac_matr.sum(0, skipna=True))


def plot_size_by_num_markers(checkm_root, output_name):

    markers = get_markers_per_contig(checkm_root)

    
def get_markers_per_contig(checkm_root):
    """
    Returns an array of markers found on each contig (contig) is the key in a dict.
    
    Actual key will not be contig but 'bin_contig' in order to separate mutiple bins with
    contigs that share the same name.
    """
    marker_f = _get_marker_file(checkm_root)
    stdout = _execute_command("checkm qa --tab_table -o 5 {} {}".format(checkm_root + "/" + marker_f, checkm_root))


    markers_per_contig = {}
    for line in stdout.strip().split("\n"):
        if line.startswith("Bin Id"):
            continue
        
        bin, marker, gene = line.split("\t")

        # the gene will be the contig with the gene number appended with a _
        contig = "_".join(gene.split("_")[:-1])
    
        try:
            markers_per_contig["{}_{}".format(bin, contig)].append(marker)
        except LookupError:
            markers_per_contig["{}_{}".format(bin, contig)] = [marker]

    return markers_per_contig

def get_sizes_per_contig(checkm_root):
    """
    Returns an array of sizes of each contig (contig) is the key in a dict.
    
    Actual key will not be contig but 'bin_contig' in order to separate mutiple bins with
    contigs that share the same name.
    """

    sizes_per_contig = {}
    for line in stdout.strip().split("\n"):
        if line.startswith("Bin Id"):
            continue
        
        bin, marker, gene = line.split("\t")

        # the gene will be the contig with the gene number appended with a _
        contig = "_".join(gene.split("_")[:-1])
    
        try:
            markers_per_contig["{}_{}".format(bin, contig)].append(marker)
        except LookupError:
            markers_per_contig["{}_{}".format(bin, contig)] = [marker]

    return markers_per_contig

def _get_marker_file(checkm_root):

    marker_f = ''
    for file in os.listdir(checkm_root):
        if file.endswith(".ms"):
            if marker_f:
                raise ValueError("Expected only one marker set in directory. Got {} and {}".format(marker_f, file))
            else:
                marker_f = file

    return marker_f


def _execute_command(command):
    """
    Executes a command and returns the output if successful, else raises an error 
    """
    return subprocess.check_output(command.split(" "), stderr=open(os.devnull, 'w'))

def plot_comp_and_contam(frame, name):
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
    parser = argparse.ArgumentParser(description="This is essentially a plotting library wrapped into a single executable.")
    parser.add_argument("-root", help="output root directory used in the CheckM run")
    parser.add_argument("-checkm", help="the checkm table from the mypyli script", required=True)
    parser.add_argument("-type", help="type of plot to make", default="all", choices=['all'])
    parser.add_argument("-base", help="base name for plots", required=True)
    parser.add_argument("-conv", help="convert taxon oids to lab names", action='store_true')
    parser.add_argument("-trim_comp", help="remove genomes below the specified completeness", type=float)
    parser.add_argument("-trim_cont", help="remove genomes above the specified contamination", type=float)
    parser.add_argument("-sort", help="sort genomes by completeness", action='store_true')
    args = parser.parse_args()

    checkm_table = read_data_table(args.checkm)

    if args.trim_comp:
        checkm_table = checkm_table.query('Completeness >= args.trim_comp')

    if args.trim_cont:
        checkm_table = checkm_table.query('Contamination <= args.trim_cont')

    if args.sort:
        checkm_table.sort_index(by='Completeness', inplace=True)

    if args.conv:
        checkm_table = isolatedb.convert_dataframe(checkm_table)

    if args.type == 'all':

        plot_comp_and_contam(checkm_table, "{}.comp_and_cont.png".format(args.base))

    #markers = get_markers_per_contig(args.r)
    #print(markers)
