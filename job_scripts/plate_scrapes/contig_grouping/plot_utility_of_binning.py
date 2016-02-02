
import argparse
import pandas
from mypyli import plotmaster

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tries to answer the question: What impact did binning the plate scrapes have on what was ultimately assembled? by looking at where the coverage for the good contigs that were assembled comes from. The idea is that if there was sufficient coverage to assemble in a single sample, binning was no help to assemble the contig.")
    parser.add_argument("-depth", help="the depth file metaBAT uses", required=True)
    parser.add_argument("-out", help="output file path for the plot [%(default)s]", default="contig_coverage_barplot.png")
    parser.add_argument("-min_cov", help="the minimum fold coverage to be considered assemble-able [%(default)s]", default=20.0, type=float)
    parser.add_argument("-min_rel", help="minimum relative coverage for a single sample to require to claim other samples added a negligible amount [%(default)s]", default=.85, type=float)
    parser.add_argument("-min_len", help="minimum length of contigs to consider [%(default)s]", default=1500, type=int)


    args = parser.parse_args()


    # read in depth file
    df = pandas.read_csv(args.depth, sep="\t", header=0, index_col=0)
    
    # remove contigs that are too short
    #df["contigLen"] = df['contigLen'].astype(int)
    df = df[df["contigLen"] >= args.min_len]

    # get a list of columns that have coverage
    cov_columns = [col for col in df.columns if col.endswith(".bam")]

    categories = {
                'asm_in_multiple': 0,       # suf. coverage to asm in multiple samples
                'asm_in_single': 0,         # suf. coverage to asm in single sample
                'asm_by_group': 0,          # would not have assembled in any single sample but does in mult
                'no_asm_single': 0,         # below asm cov but most comes from a single sample
                'no_asm_multi': 0,          # below asm cov but cov comes from multiple samples
                }
    # loop over each row in the dataframe
    for index, row in df.iterrows():
        total_cov = row['totalAvgDepth']

        # contig had total cov required to asm
        if total_cov >= args.min_cov:

            # count the number of samples that had sufficient cov
            asm_count = 0
            for col in cov_columns:
                if row[col] >= args.min_cov:
                    asm_count += 1

            if asm_count > 1:
                categories['asm_in_multiple'] += 1
            elif asm_count == 1:
                categories['asm_in_single'] += 1
            else:
                categories['asm_by_group'] += 1

        else:
            # check if most of the cov is from a single sample
            for col in cov_columns:
                if row[col] / total_cov >= args.min_rel:
                    #print((total_cov, row[col]))
                    categories['no_asm_single'] += 1
                    break       # no need to check the rest

            else:
                categories['no_asm_multi'] += 1


    # convert to series and make a barplot
    category_s = pandas.Series(data=categories)
    
    # sort the series like I want 
    category_s = category_s[["asm_in_multiple", "asm_in_single", "asm_by_group", "no_asm_single", "no_asm_multi"]]

    ax = category_s.plot(kind="bar")
    ax.set_xlabel("Category")
    ax.set_ylabel("Number of Contigs (>= {}bp)".format(args.min_len))
    ax.set_title("Impact of Binning on Contig Coverage")
    plotmaster.save_ax(ax, args.out)
