
import argparse
import pandas as pd
import matplotlib.pyplot as plt

def read_data(data_f):
    return pd.read_csv(data_f, sep="\t", header=0, index_col=0)

def reformat_table(table, field):
    #print(table)
    new_table = pd.DataFrame()
    for index in table.index.values:
        genome, cov = index.split("_cov")
        new_table.loc[genome, cov] = table.loc[index, field]

    #print(new_table)

    return(new_table)

def plot(table, title):
    table = table.sort_index(axis=1)
    
    ax = table.T.plot(sort_columns=True, title=title)
    ax.set_xlabel("assembly (number is coverage)")
    plt.xticks(range(len(table.columns.values)), table.columns.values)
    fig = ax.get_figure()
    fig.savefig("{}.png".format(title))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-checkm", help="the table output from mypyli's checkm script")
    parser.add_argument("-genes", help="a table of percentage of genes found")
    parser.add_argument("-tquast", help="transposed quast .tsv file")
    args = parser.parse_args()

    if args.checkm:

        table = read_data(args.checkm)
        reform_table = reformat_table(table, field='completeness')
        plot(reform_table, "checkm_completeness")

    if args.genes:
        table = read_data(args.genes)
        reform_table = reformat_table(table, field='present')
        plot(reform_table, "perc_genes_present_above_95id")

        reform_table = reformat_table(table, field='ordered')
        plot(reform_table, "perc_genes_ordered")

    if args.tquast:
        table = read_data(args.tquast)
    
        reform_table = reformat_table(table, field='N50')
        plot(reform_table, "quast_N50")

        reform_table = reformat_table(table, field='# contigs')
        plot(reform_table, "quast_num_contigs")

        reform_table = reformat_table(table, field='Total length')
        plot(reform_table, "quast_genome_length")

        reform_table = reformat_table(table, field='L50')
        plot(reform_table, "quast_L50")
