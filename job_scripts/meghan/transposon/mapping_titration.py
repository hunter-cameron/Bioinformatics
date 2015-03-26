
import pandas as pd
import sys
import matplotlib.pyplot as plt

headers = ['gene_cov', 'gene_hits', 'num_mapped_Q>20', 'uniq_genes', 'uniq_pos']
stats_f = sys.argv[1]
title = sys.argv[2]
name = "titration_figure.png"

frame = pd.read_csv(stats_f, sep="\t", header=None, names=headers, index_col=0)

series = frame['uniq_genes']

series.index = [int(indx.split("genomic")[1]) for indx in series.index]

series = series.sort_index()
    
print(series)

ax = series.T.plot(title=title)
ax.set_xlabel("% of filtered reads")
ax.set_ylabel("# of unique genes")
#plt.xticks(range(len(series.name)), series.name)
ax.set_ylim(ymin=0)
fig = ax.get_figure()
fig.savefig(name)
