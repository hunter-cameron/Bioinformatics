#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import pandas
import numpy as np

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--otu_table", "-in", help="Matrix with samples as rows and OTU's as columns", type=argparse.FileType('r'), required=True)
    parser.add_argument("-r", "--rarefy", help="Number to rarefy to. Samples with fewer reads than the rarefaction will be reduced to 0.", type=int)
    parser.add_argument("-t", "--transpose", help="Use this flag to denote that the OTU table has OTU's as rows and samples as columns", action="store_true")
    parser.add_argument("-f", "--format_matrix", help="Matrix with samples as rows and experimental groups as columns", type=argparse.FileType('r'))
    parser.add_argument("-d", "--dendrogram", help="Cluster by dendogram on samples, OTUs, or both", choices=["samples", "OTUs", "both"])

    args = parser.parse_args()
    print(args)


    #####
    #Error Checking for the input files, etc.
    #####

    if (args.format_matrix):
        parse_format(args.format_matrix)

    data = parse_OTU_table(args.otu_table)

    make_heatmap(data)


def parse_OTU_table(file):
    data = pandas.read_csv(file, sep="\t", header=0, index_col=0)
    print(data)
    return data

def parse_format(format_matrix):
    print("Format matrix option hasn't been implemented yet.")

def rarefy(data, r=0, seed=0):
    data = np.matrix(data)

    

def make_heatmap(data):
    #make a figure prepared with subplot
    print(data)
    figure, sub = plt.subplots()

    #make the actual heatmap
    heatmap = sub.pcolor(data, cmap=plt.cm.Blues, alpha=0.8)

    figure = plt.gcf()     #.gcf() returns a reference to the current figure
    figure.set_size_inches(8, 11)

    #turn off frame
    sub.set_frame_on(False)

    #put major ticks in middle. why?
    sub.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
    sub.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)

    #more natural table-like display?
    sub.invert_yaxis()
    sub.xaxis.tick_top()

    #set labels
    #list of labels

    sub.set_xticklabels(data.columns, minor=False)
    sub.set_yticklabels(data.index, minor=False)

    #rotate x labels
    plt.xticks(rotation=90)

    sub.grid(False)

    #Turn off all ticks
    sub = plt.gca()

    for t in sub.xaxis.get_major_ticks():
        t.tick10n = False
        t.tick20n = False
    for t in sub.yaxis.get_major_ticks():
        t.tick10n = False
        t.tick20n = False

    plt.show()


if __name__ == "__main__":
    main()




"""
import matplotlib.pyplot as plt
import pandas as pd
from urllib2 import urlopen
import numpy as np
%pylab inline

page = urlopen("http://datasets.flowingdata.com/ppg2008.csv")
nba = pd.read_csv(page, index_col=0)

# Normalize data columns
nba_norm = (nba - nba.mean()) / (nba.max() - nba.min())

# Sort data according to Points, lowest to highest
# This was just a design choice made by Yau
# inplace=False (default) ->thanks SO user d1337
nba_sort = nba_norm.sort('PTS', ascending=True)

nba_sort['PTS'].head(10)

# Plot it out
fig, ax = plt.subplots()
heatmap = ax.pcolor(nba_sort, cmap=plt.cm.Blues, alpha=0.8)

# Format
fig = plt.gcf()
fig.set_size_inches(8, 11)

# turn off the frame
ax.set_frame_on(False)

# put the major ticks at the middle of each cell
ax.set_yticks(np.arange(nba_sort.shape[0]) + 0.5, minor=False)
ax.set_xticks(np.arange(nba_sort.shape[1]) + 0.5, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

# Set the labels

# label source:https://en.wikipedia.org/wiki/Basketball_statistics
labels = [
    'Games', 'Minutes', 'Points', 'Field goals made', 'Field goal attempts', 'Field goal percentage', 'Free throws made', 'Free throws attempts', 'Free throws percentage',
    'Three-pointers made', 'Three-point attempt', 'Three-point percentage', 'Offensive rebounds', 'Defensive rebounds', 'Total rebounds', 'Assists', 'Steals', 'Blocks', 'Turnover', 'Personal foul']

# note I could have used nba_sort.columns but made "labels" instead
ax.set_xticklabels(labels, minor=False)
ax.set_yticklabels(nba_sort.index, minor=False)

# rotate the
plt.xticks(rotation=90)

ax.grid(False)

# Turn off all the ticks
ax = plt.gca()

for t in ax.xaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False
for t in ax.yaxis.get_major_ticks():
    t.tick1On = False
    t.tick2On = False

"""

