
from __future__ import print_function
import argparse
import sys
import pandas as pd
import subprocess
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



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="This is essentially a plotting library wrapped into a single executable.")
    parser.add_argument("-c", help="output directory used in the CheckM run")
    args = parser.parse_args()

    markers = get_markers_per_contig(args.c)
    print(markers)
