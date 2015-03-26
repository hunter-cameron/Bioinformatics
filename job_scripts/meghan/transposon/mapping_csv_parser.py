

import pandas as pd
import argparse
import sys
import os

def read_mapping_file(map_f):
    frame = pd.read_csv(map_f, sep="\t", header=0)

    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", nargs="+", help="one or more mapping.csv files", required=True)
    args = parser.parse_args()

    

    frames = {}
    for map_file in args.m:
        # this current method of parsing the name will not work 
        # for relative paths
        frame_name = os.path.basename(map_file.split(".")[0])
        frames[frame_name] = read_mapping_file(map_file)
