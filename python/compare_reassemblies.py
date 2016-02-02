
import argparse
import os
import sys
import pandas
from mypyli import plotmaster

def pair_genomes(quast_df):
    """ Returns a list of tuples with paired genomes. """
    
    # if the names are the same except with .orig and .new, to pair, I can sort
    pairs = []
    name_iter = iter(sorted(quast_df.index))
    while True:
        try:
            new = next(name_iter)
        except StopIteration:
            break

        try:
            orig = next(name_iter)
        except StopIteration:
            raise ValueError("There are an odd number of records in your quast table. Obviously they are not paired...")

        basename1 = orig.rsplit(".", 1)[0]
        basename2 = new.rsplit(".", 1)[0]

        if basename1 != basename2:
            raise ValueError("Names '{}' and '{}' were attempted to be paired but they don't appear to be the same. Make sure your pairs sort next to one another.".format(orig, new))

        pairs.append((orig, new))


    return pairs

def make_data_for_plot(df, pairs, fields=['Total length', 'N50', 'L50', '-# contigs', 'Largest Contig']):
    """ Makes a df with all the data, the fields list of a list of quast fields. Fields where better == smaller should be prepended with a '-' """

    data_dict = {}

    for orig, new in pairs:
        basename = orig.rsplit(".", 1)[0]
        data_dict[basename] = {}

        # get the fold difference between new and original
        fold_diff = (df.loc[new] - df.loc[orig]) / df.loc[orig]

        # collect the sort order, it will do this for each pair. But I don't care that much...
        sort_order = []
        for field in fields:
        # for fields where smaller is better, we want to inverse to keep all situations where the new asm is better > 0
            if field.startswith("-"):
                inverse = True
                field = field[1:]
            else:
                inverse = False

            value = fold_diff[field]

            if inverse:
                value = value * -1

            sort_order.append(field)
            data_dict[basename][field] = value

    asm_comp_df = pandas.DataFrame.from_dict(data_dict, orient='index')
   
    # sort the df by the fields given
    asm_comp_df = asm_comp_df[sort_order]


    # natural sort the df index
    def magic_sort(key):
        elems = key.split(".")

        sort_list = []
        for elem in elems:
            try:
                elem = int(elem)
            except ValueError:
                elem = elem.lower()

            sort_list.append(elem)

        return sort_list

    asm_comp_df = asm_comp_df.reindex(sorted(asm_comp_df.index, key=magic_sort))

    return asm_comp_df

def natural_sort(items):
        """ Returns an ordered list by a natural sort algorithm  """

        def convert(char):
            """ Attempts to convert a character into an integer """
            try:
                return int(char)
            except ValueError:
                return char

        def nat_sort(entry):
            """ Performs a natural sort that will sort text and numbers in a way that makes sense """
            return [convert(char) for char in entry]


        return sorted(items, key=nat_sort)

def plot_asm_comp(asm_comp_df, out):
    """ Makes a plot of the assembly comparison DataFrame """
    ax = plotmaster.scatter_plot(asm_comp_df)

    ax.axhline(y=0, linestyle="dashed")

    ax.set_title("Assembly Comparison")
    ax.set_xlabel("Genome")
    ax.set_ylabel("Fold Change")

    plotmaster.save_ax(ax, out)

def add_gene_counts_to_df(df, gene_files):
    """ Adds gene counts to the df by matching gene file basenames to index fields """
    gene_counts = {}
    for gene_file in gene_files:
        gene_base = os.path.splitext(os.path.basename(gene_file))[0]

        if gene_base in df.index:
            gene_count = 0
            with open(gene_file, 'r') as IN:
                for line in IN:
                    if line.startswith(">"):
                        gene_count += 1

            gene_counts[gene_base] = gene_count

        else:
            raise ValueError("No match was found for gene file '{}' in the quast table.".format(gene_file))

    # add the counts to the df
    count_list = []
    for name in df.index:
        try:
            count_list.append(gene_counts[name])
        except KeyError:
            raise ValueError("Genes file was not found for quast entry '{}'.".format(name))

    df["Gene count"] = count_list
    return df

def add_checkm_to_df(df, checkm_f):
    """ Adds checkm information to the data frame and returns """

    checkm_df = pandas.read_csv(checkm_f, sep="\t", index_col=0, header=0)

    if set(df.index) != set(checkm_df.index):
        raise ValueError("The checkm rows do not match the quast rows.")

    df['Completeness'] = checkm_df['Completeness']
    df['Contamination'] = checkm_df['Contamination']

    return df
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compares pairs of genomes on fields like N50, Genome length, longest contig, number of genes, etc.  Genomes will be paired by looking at the names in the quast file. It is expected to find a pair of entries with the same basename, one with '.orig' and one with '.new'. The plot is a fold-change plot where anything > 0 means the new assembly is better than the orig. (For fields like # contigs where fewer is better, the fold change is inversed to keep everything above the line meaning the new asm is better.)")
    parser.add_argument("-quast", help="transposed (genomes as rows) quast results file that includes all genomes", required=True)
    parser.add_argument("-genes", help=".faa files (one per assembly) with annotated genes", nargs="+")
    parser.add_argument("-checkm", help="checkm tab delimited file")
    parser.add_argument("-out", help="the path for the figure [%(default)s]", default="asm_comparison.png")
    args = parser.parse_args()

    # build the default df and set the default fields
    df = pandas.read_csv(args.quast, sep="\t", header=0, index_col=0)
    pairs = pair_genomes(df)
    fields = ['Total length', 'N50', '-L50', '-# contigs', 'Largest contig']

    # add genes if supplied
    if args.genes:
        df = add_gene_counts_to_df(df, args.genes)
        fields += ['Gene count']

    # add checkm data if supplied
    if args.checkm:
        df = add_checkm_to_df(df, args.checkm)
        fields += ['Completeness', '-Contamination']

    asm_comp_df = make_data_for_plot(df, pairs, fields=fields)
    plot_asm_comp(asm_comp_df, args.out)
