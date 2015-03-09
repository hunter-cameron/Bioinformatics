
import argparse
import sys
import pandas as pd


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

def main(ref_f, queries_f):
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compares CheckM outputs")
    parser.add_argument("-ref", help="reference CheckM dataset, this should include the same isolates as all other sets", required=True)
    parser.add_argument("-query", help="one or more CheckM datasets to compare to the ref", nargs="+", required=True)
    args = parser.parse_args()

    main(args.ref, args.query)
