
import pandas
import sys
import argparse
import matplotlib.pyplot as plt

from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

def plot_trait_categories(df):
    """ Makes plots of each genome by NSTI for a category """

    for col in df.columns:
        if col == "NSTI":
            continue
        ax = df.plot(x="NSTI", y=col, title=col, kind="scatter", xlim=(-.01, max(df["NSTI"]) + .01), ylim=(-.1, 1.1)) 

        ax.set_xlabel("NSTI")
        ax.set_ylabel("Spearman Correlation")
        fig = ax.get_figure()
        fig.savefig(str(col).replace(" ", "_") + ".correlation.png")
        fig.close()

def plot_trait_panel_plot(df):
    """ Makes a pannelplot for each level 1 KO group with each level 2 as a pannel """

    # determine how many plots we need
    plots = {}
    for name in df.columns:
        try:
            lvl1, lvl2 = name.split(";")

        except ValueError:
            continue

        try:
            plots[lvl1].append(lvl2)
        except KeyError:
            plots[lvl1] = [lvl2]

    for plot_name in plots:
        sub_names = sorted(plots[plot_name])

        if len(sub_names) == 1:
            f, ax = plt.subplots()
            ax_matr = [ax]

        elif len(sub_names) == 2:
            f, ax_matr = plt.subplots(1, 2, sharey=True)

        elif len(sub_names) <= 4:
            f, ax_matr = plt.subplots(2, 2, sharex=True, sharey=True)

        elif len(sub_names) <= 9:
            f, ax_matr = plt.subplots(3, 3, sharex=True, sharey=True)

        elif len(sub_names) <= 12:
            f, ax_matr = plt.subplots(4, 3, sharex=True, sharey=True)

        elif len(sub_names) <= 16:
            f, ax_matr = plt.subplots(4, 4, sharex=True, sharey=True)

        else:
            raise ValueError("Too many subplots required.")

        ax_arr = ax_matr.flatten()
        for index, subname in enumerate(sub_names):
            ax = ax_arr[index]


            ax.scatter(x=df["NSTI"], y=df[";".join([plot_name, subname])])
            ax.set_title(subname)
        
            #ax.set_xlabel("Trait")
            #ax.set_xticks(xvals)
            #ax.set_xticklabels(df.index, fontsize="x-small", rotation="vertical", ha="center")
        
            ax.set_xlim(-.01, max(df["NSTI"]) + .01)
            ax.set_ylim(-.1, 1.1)
        
            #ax.set_ylabel("Spearman Correlation")
       
            if plot_name == "Metabolism":
                f.set_size_inches(15, 15)
            elif plot_name == "Human Diseases":
                f.set_size_inches(15, 15)

            f.savefig(str(plot_name).replace(" ", "_") + ".correlation.png")
            plt.close(f)


def plot_genomes(df):
    """ Makes plots of the correlation for each category for each genome """

    # drop the metadata and transpose
    df = df.drop("NSTI", 1).transpose()
    xvals = range(len(df.index))
    for col in df.columns:

        fig = plt.figure()

        ax = fig.add_subplot(111)

        ax.scatter(x=xvals, y=df[col])
        ax.set_title(col)
        
        ax.set_xlabel("Trait")
        ax.set_xticks(xvals)
        ax.set_xticklabels(df.index, fontsize="x-small", rotation="vertical", ha="center")
        
        ax.set_xlim(-1, xvals[-1]+1)
        
        ax.set_ylabel("Spearman Correlation")
        
        fig.savefig(str(col).replace(" ", "_") + ".correlation.png")
        plt.close(fig)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Makes plots from PICRUSt correlation tables.")
    parser.add_argument("-table", help="the correlation table; should have NSTI as a field", required=True)
    parser.add_argument("-rows", help="make plots from the rows", action="store_true")
    parser.add_argument("-cols", help="make plots from the cols", action="store_true")
    args = parser.parse_args()

    df = pandas.read_csv(args.table, sep="\t", index_col=0, header=0)

    if args.rows:
        plot_genomes(df)

    if args.cols:
        plot_trait_panel_plot(df)
