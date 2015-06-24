
import numpy as np
import matplotlib.pyplot as plt

def histogram(data, bins):
    """ 
    Returns a matplotlib ax object or a histogram, specializes in 
    making histograms with uneven bins
    """
    
    print(data[0:99])

    hist, h_bins = np.histogram(data[0:99], bins=bins)
    fig = plt.figure()
    ax = fig.add_subplot(111)

    print(hist)
    print(h_bins)

    ax.bar(range(len(hist)), height=hist, align='edge', width=1)
    ax.set_xlim([0, len(hist)])       # x vals with a 1 pad on each side
    ax.set_xticks(range(len(h_bins)))
    #ax.set_xticks([0] + [str(bin) for bin in bins])
    ax.set_xticklabels([0] + [str(bin) for bin in bins])
    
    return ax


def stacked_bar_plot(data, stack_cols, colors=["blue", "green", "red", "cyan", "magenta", "yellow", "black", "white"], log=False):
    """ 
    Returns an ax with a stacked bar plot using columns in stack cols
    (in the order they appear in stack cols - should be tallest first)
    """
    x_vals = range(data.shape[0])
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for indx, col in enumerate(stack_cols):
        if log:
            ax.bar(x_vals, data[col], color=colors[indx], align='center', label=col, log=True)
        else:
            ax.bar(x_vals, data[col], color=colors[indx], align='center', label=col)

    names = data.index.values

    plt.xlim(x_vals[0] - 1, x_vals[-1] + 1)
    plt.xticks(x_vals, names, rotation="vertical", ha='center', fontsize='x-small')

    return ax

def savefig(ax, filename, height=8.5, width=11):
    fig = ax.get_figure()
    fig.tight_layout()
    fig.savefig(filename)
    #fig.savefig(filename, bbox_extra_artists=fig.artists, bbox_inches='tight')
