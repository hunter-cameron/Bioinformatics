

import pandas as pd
from matplotlib_venn import venn2_unweighted, venn3_unweighted
import matplotlib.pyplot as plt
import argparse
import sys
import os

def read_mapping_file(map_f):
    """ Reads a mapping.csv file and returns a pandas dataframe. """
    frame = pd.read_csv(map_f, sep="\t", header=0, na_values=['None'])

    return frame

def read_sample_sheet(samp_f):
    """ Reads the sample sheet and returns a dict of sample and name. """
    sample_dict = {}
    with open(samp_f, 'r') as IN:
        for line in IN:
            sample, name = line[:-1].split("\t")
            sample_dict[sample] = name
    
    return sample_dict

def bin_samples(frames, how="both"):
    """ 
    Returns pooled data frames as a new dict (labeled with pool names)
    Bins based on how which is one of ['both', 'chamber', 'site']
    Chamber is either gnoto or open.
    Site is either root or soil.
    Both is gnot_soil, gnot_root, open_soil, open_root.
    """
    bins = {}
    for sample in frames:
        if "_" in sample:   # checks for libraries, unknown
            code, site, chamber = sample.split("_")
        else:
            continue

        # find the bin type
        if how == "both":
            s_type = "_".join([site, chamber]).lower()
        elif how == "site":
            s_type = site.lower()
        elif how == "chamber":
            s_type = chamber.lower()
        else:
            raise ValueError("How must be one of 'both', 'site', 'chamber'")

        print("{} classified as {}".format(sample, s_type))

        cur_frame = bins.get(s_type, None)
        if cur_frame is not None:
            bins[s_type] = cur_frame.append(frames[sample], ignore_index=True)
        else:
            bins[s_type] = frames[sample]

    return bins

def rewrite_sample_sheet(samp_f):
    """ Rewrites the sample sheet to have a consistent code_type_chamber format"""
    with open(samp_f, 'r') as IN, open("new_samp_sheet.txt", 'w') as OUT:
        for line in IN:
            sample, name = line[:-1].split("\t")
            if "_" in name:
                code, pt1, pt2 = name.split("_")
                if pt1[0] in ["S", "R", "s", "r"]:
                    name = "_".join([code, pt1, pt2])
                else:
                    name = "_".join([code, pt2, pt1])
            OUT.write("\t".join([sample, name]) + "\n")

def make_gene_table(frames, names):
    """ Makes a table of gene counts for all the frames in frames. """

    fr_dict = {}
    for frame, name in zip(frames, names):
        counts = frame['gene'].value_counts(sort=False)
        fr_dict[name] = counts

    comb_frame = pd.concat(fr_dict, axis=1)
    comb_frame.fillna(0, inplace=True)

    comb_frame.to_csv("gene_table.csv", sep="\t", index_label="gene")



    


def venn_diagram(frame1, frame2, field, names, frame3=None):
    """ 
    Makes a Venn diagram for two frames (optionally 3) based on the
    field specified. Return as axis object. 
    """
    set1 = set(frame1[field].dropna())
    #print(sorted(set1))
    #print()
    set2 = set(frame2[field].dropna())
    #print(sorted(set2))
    #print()

    if frame3 is None:
        set3 = None
    else:
        set3 = set(frame3[field].dropna())
        #print(sorted(set3))

    if set3:
        if len(names) != 3:
            raise ValueError("Length of names does not match the number of frames.")

        venn3_unweighted([set1, set2, set3], tuple(names))

    else:
        if len(names) != 2:
            raise ValueError("Length of names does not match the number of frames.")

        venn2_unweighted([set1, set2], tuple(names))
    
    return plt.gca()

def make_venn_diagrams(frames):
    ### this is just to make a bunch of venn diagrams
    
    # small lib to big lib 
    names = ["PfTn9", "PfTn10", "PfTnYel"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=frames[names[1]], frame3=frames[names[2]], names=names, field="gene")
    ax.set_title("Partial Libraries to Full Library")
    ax.get_figure().savefig("lib9_10_Yel.venn_diagram.genes.png")
    plt.clf()

    # STPCR to standard
    names = ["PfTnYel", "PfTnYelSTPCR"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=frames[names[1]], names=names, field="gene")
    ax.set_title("Standard PCR to STPCR")
    ax.get_figure().savefig("Yel_YelSTPCR.venn_diagram.genes.png")
    plt.clf()

    # STPCR to standard to SAMPLE lib
    names = ["PfTnYel", "PfTnYelSTPCR", "PfTnInput"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=frames[names[1]], frame3=frames[names[2]], names=names, field="gene")
    ax.set_title("Full Input Libraries")
    ax.get_figure().savefig("Yel_YelSTPCR_Input.venn_diagram.genes.png")
    plt.clf()

    
    ## Pool the samples
    bins = bin_samples(frames)

    # library to GNOT samples
    names = ["PfTnYel", "soil_gnot", "root_gnot"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=bins[names[1]], frame3=bins[names[2]], names=names, field="gene")
    ax.set_title("Library to Gnotobiotic Samples")
    ax.get_figure().savefig("Yel_gnoto.genes.png")
    plt.clf()

    # library to open samples
    names = ["PfTnYel", "soil_open", "root_open"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=bins[names[1]], frame3=bins[names[2]], names=names, field="gene")
    ax.set_title("Library to Open Samples")
    ax.get_figure().savefig("Yel_open.genes.png")
    plt.clf()

    # library to soil samples
    names = ["PfTnYel", "soil_gnot", "soil_open"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=bins[names[1]], frame3=bins[names[2]], names=names, field="gene")
    ax.set_title("Library to Soil Samples")
    ax.get_figure().savefig("Yel_soil.genes.png")
    plt.clf()

    # library to root samples
    names = ["PfTnYel", "root_gnot", "root_open"]
    ax = venn_diagram(frame1=frames[names[0]], frame2=bins[names[1]], frame3=bins[names[2]], names=names, field="gene")
    ax.set_title("Library to Root Samples")
    ax.get_figure().savefig("Yel_root.genes.png")
    plt.clf()



    #random = ['C1_Root_GNOT', 'A5_Root_GNOT', 'C2_root_GNOT']

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", nargs="+", help="one or more mapping.csv files", required=True)
    parser.add_argument("-s", help="sample sheet", required=True)
    args = parser.parse_args()


    #rewrite_sample_sheet(args.s)
    #sys.exit()
    sample_dict = read_sample_sheet(args.s)

    frames = {}
    for map_file in args.m:
        # get the basename
        sample = os.path.basename(map_file).split(".")[0]
        # remove the well code
        sample = sample.split("_")[0]
        
        # lookup the real name
        sample_name = sample_dict[sample]
        
        # index the dataframe with the name
        frames[sample_name] = read_mapping_file(map_file)

    ### AT THIS POINT, THE ENVIRONMENT IS SET UP FOR INTERACTIVE MODE

    bins = bin_samples(frames)
    for_table = [frames['PfTnYel'], bins['soil_gnot'], bins['root_gnot'], bins['soil_open'], bins['root_open']]
    counts = make_gene_table(for_table, names=['PfTnYel', 'soil_gnot', 'root_gnot', 'soil_open', 'root_open'])
    sys.exit()
    
    
    [print(frame) for frame in frames]

    make_venn_diagrams(frames)


