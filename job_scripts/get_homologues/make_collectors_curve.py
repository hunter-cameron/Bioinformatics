
import argparse
from mypyli import homologues, plotmaster
from matplotlib import pyplot as plt


def subset_df(df, subset_f, inverse):
    
    # read the names
    subset_names = []
    with open(subset_f, 'r') as IN:
        for line in IN:
            subset_names.append(line.strip())

    if inverse:
        df.drop(subset_names, inplace=True)

    else:
        df = df.loc[subset_names]

    return df

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-matrix", help="GH pangenome matrix", required=True)
    parser.add_argument("-mode", help="mode to run in. either 'total' or 'core'", choices=["total", "core"], default="total")
    parser.add_argument("-subset", help="a list of genomes to include")
    parser.add_argument("-inverse", help="makes subset of list of genomes to remove", action="store_true")
    parser.add_argument("-out", help="name for the output figure", default="collectors_curve.png")

    args = parser.parse_args()

    # make df from pangenome matrix
    pg_matr_df = homologues.read_pangenome_matrix(args.matrix)

    # optionally filter the matrix
    if args.subset:
        pg_matr_df = subset_df(pg_matr_df, args.subset, args.inverse)

    if args.mode == "total":
        # generate a dataset to plot the collectors curve
        datapoints = homologues.get_bootstrapped_collectors_curve_datapoints(pg_matr_df)
    elif args.mode == "core":
        datapoints = homologues.get_bootstrapped_core_curve_datapoints(pg_matr_df)

    means = datapoints.mean()
    errors = datapoints.std()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    means.plot(yerr=errors, ax=ax, kind="line")
    ax.set_ylabel("# gene clusters")
    ax.set_xlabel("# genomes")
    fig.savefig(args.out)
