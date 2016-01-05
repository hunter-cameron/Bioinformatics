
import argparse
import pandas
import numpy as np
import sys
import random

def read_pangenome_matrix(matrix_f):
    """ Reads a pangenome matrix into a df """

    df = pandas.DataFrame.from_csv(matrix_f, header=0, sep="\t", index_col=0)
    df.fillna(0, inplace=True)
    return df


def get_bootstrapped_collectors_curve_datapoints(pg_matr_df, bootstrap=100):
    """ Gets datapoints for a bootstrapped collectors curve """

    # make boolean version of df
    boolean_masks = {}
    for indx in pg_matr_df.index:
        boolean_masks[indx] = pg_matr_df.loc[indx, :].astype(bool)
    
    #
    ## Build an array of collectors curve datapoints
    #
    bootstrapped_data = {}
    for i in range(1, bootstrap+1):
        # make an all False collector's mask
        collectors_mask = np.array([False for count in range(len(pg_matr_df.columns))], dtype=bool)

        # make a copy of boolean masks
        b_masks = boolean_masks.copy()

        bootstrapped_data[i] = {}

        added_genomes = 0
        while b_masks:
            #print("Genomes left: {}".format(len(b_masks)))

            selected_mask = random.choice(list(b_masks.keys()))
            added_genomes += 1
            
            collectors_mask = np.logical_or(collectors_mask, b_masks[selected_mask])
  
            bootstrapped_data[i][added_genomes] = np.sum(collectors_mask)

            del b_masks[selected_mask]


    # convert the nested dict data structure to a df
    df = pandas.DataFrame.from_dict(bootstrapped_data, orient="index")
    return df

def get_collectors_curve_datapoints(pg_matr_df):
    """ Makes an optimal collectors curve using numpy boolean arrays. """

    # make boolean version of df
    boolean_masks = {}
    for indx in pg_matr_df.index:
        boolean_masks[indx] = pg_matr_df.loc[indx, :].astype(bool)
    
    #
    ## Build an array of collectors curve datapoints
    #

    # make an all False collector's mask
    collectors_mask = np.array([False for count in range(len(pg_matr_df.columns))], dtype=bool)
    
    # keep counts on how many clusters each genome adds
    collectors_data = []
    while boolean_masks:
        print("Genomes left: {}".format(len(boolean_masks)))

        # tuple to store (# additions, key)
        best = (0, None)

        # pick the best genome
        for key in boolean_masks:
            # compare to collectors mask by getting the sum of a bool array of test > mask
            # In numpy True = 1 and False = 0
            additions = np.sum(boolean_masks[key] > collectors_mask)
            
            # update best
            if additions >= best[0]:
                best = (additions, key)

        # update the mask with the best and remove the best
        #print(("collectors", list(collectors_mask)))
        #print(("addition..", list(boolean_masks[best[1]])))
        collectors_mask = np.logical_or(collectors_mask, boolean_masks[best[1]])
        #print(("new collec", list(collectors_mask)))
        #print(np.sum(collectors_mask))
        #print()
        collectors_data.append(np.sum(collectors_mask))
        
        del boolean_masks[best[1]]

    return collectors_data


def get_core_curve_datapoints(pg_matr_df):
    """ Makes an optimal curve of core clusters using numpy boolean arrays 
    
    NOTE: this curve may not be completely optimal because the genomes with the most clusters (the starting genome) may not have the most core clusters.
    
    """

    # make boolean version of df
    boolean_masks = {}
    for indx in pg_matr_df.index:
        boolean_masks[indx] = pg_matr_df.loc[indx, :].astype(bool)
    
    #
    ## Build an array of core curve datapoints
    #

    # make an all True collector's mask
    core_mask = np.array([True for count in range(len(pg_matr_df.columns))], dtype=bool)
    
    # keep counts on how many core clusters each genome removes
    core_data = []
    while boolean_masks:
        print("Genomes left: {}".format(len(boolean_masks)))

        # tuple to store (# core, key)
        best = (0, None)

        # pick the best genome
        for key in boolean_masks:
            # compare to core mask by getting the sum of a bool array of test > mask
            # In numpy True = 1 and False = 0
            num_core = np.sum(np.logical_and(boolean_masks[key], core_mask))
            print(num_core)
            
            # update best
            if num_core >= best[0]:
                best = (num_core, key)

        # update the mask with the best and remove the best
        #print(("core", list(core_mask)))
        #print(("addition..", list(boolean_masks[best[1]])))
        core_mask = np.logical_and(core_mask, boolean_masks[best[1]])
        #print(("new core", list(core_mask)))
        #print(np.sum(core_mask))
        #print()
        core_data.append(np.sum(core_mask))
        
        del boolean_masks[best[1]]

    return core_data

def get_bootstrapped_core_curve_datapoints(pg_matr_df, bootstrap=100):
    """ Gets datapoints for a bootstrapped collectors curve """

    # make boolean version of df
    boolean_masks = {}
    for indx in pg_matr_df.index:
        boolean_masks[indx] = pg_matr_df.loc[indx, :].astype(bool)
    
    #
    ## Build an array of collectors curve datapoints
    #
    bootstrapped_data = {}
    for i in range(1, bootstrap+1):
        # make an all False collector's mask
        core_mask = np.array([True for count in range(len(pg_matr_df.columns))], dtype=bool)

        # make a copy of boolean masks
        b_masks = boolean_masks.copy()

        bootstrapped_data[i] = {}

        added_genomes = 0
        while b_masks:
            #print("Genomes left: {}".format(len(b_masks)))

            selected_mask = random.choice(list(b_masks.keys()))
            added_genomes += 1
            
            core_mask = np.logical_and(core_mask, b_masks[selected_mask])
  
            bootstrapped_data[i][added_genomes] = np.sum(core_mask)

            del b_masks[selected_mask]


    # convert the nested dict data structure to a df
    df = pandas.DataFrame.from_dict(bootstrapped_data, orient="index")
    return df

def find_duplicates(pg_matr_df, min_diff=1000):
    """ Finds duplicate genomes based on a pangenome matrix """

    # calculate distance in terms of gene clusters using a lower right triangle algorithm
    genomes = pg_matr_df.index
    differences = {genome: {} for genome in genomes}

    # tests = n**2 + 2 / 2 + n (+ n because we do not need to test self)
    tests_to_run = ((len(genomes) ** 2 + len(genomes)) / 2) - len(genomes)
    tests_ran = 0
    print("Beginning to look for duplicates.", file=sys.stderr)
    for g1 in genomes:
        for g2 in genomes:
            if g1 == g2:
                break

            difference = 0
            for c1, c2 in zip(pg_matr_df.loc[g1, :], pg_matr_df.loc[g2, :]):
                if c1 != c2:
                    difference += 1

            differences[g1][g2] = difference

            tests_ran += 1
            # small progress indicator if running over 1000 tests
            #print(tests_ran, file=sys.stderr)
            if tests_ran % 1000 == 0:
                print("Progress...{}%".format((tests_ran / tests_to_run) * 100), end="\r", file=sys.stderr)

            #print((g1, g2, difference))

    duplicates = []
    for g1 in differences:
        for g2 in differences[g1]:
            if g1 == g2:
                continue
            if differences[g1][g2] < min_diff:
                duplicates.append((g1, g2))

    return duplicates




if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-matrix", help="a GH pangenome matrix")

    args = parser.parse_args()

    pg_matr_df = read_pangenome_matrix(args.matrix)

    #data = get_collectors_curve_datapoints(pg_matr_df)
    data = find_duplicates(pg_matr_df)
    print(data)
