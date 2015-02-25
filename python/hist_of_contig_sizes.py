

import sys
import numpy
import matplotlib.pyplot as plt

def get_size_per_contig(fasta):
    header = ""
    seq = ""
    all_lens = []
    with open(fasta, 'r') as IN:
        for line in IN:
            if line.startswith(">"):
                # make the tuple
                if not header == "":
                    all_lens.append((header, len(seq)))
                
                header = line[:-1]
                seq = ""
            else:
                seq += line[:-1]

    return all_lens


def make_histogram(array):
    plt.hist(array, bins=(0, 1000, 10000, 100000, 500000, 1000000, 5000000))
    ax = plt.gca()
    ax.set_autoscale_on(False)
    plt.title("Contig Lengths")
    plt.xlabel("Length")
    plt.ylabel("Count")
    plt.show()

def make_text_hist(array):
    stdout = sys.stdout
    sys.stdout = sys.stderr
    (hist, bin_edges) = numpy.histogram(array, bins=(0, 1000, 5000, 10000, 100000, 500000, 1000000), density=False)
    
    sys.stdout = stdout
    print("\t".join(["count", "bin"]))
    for count, bin in zip(hist, bin_edges[1:]):
        print("\t".join([str(count), str(bin)]))
    return(hist, bin_edges)
    
    




if __name__ == "__main__":
    len_info = get_size_per_contig(sys.argv[1])

    # generate a list of just sizes from tuple
    lengths = [seq[1] for seq in len_info]

    print("Total base pairs: {}".format(str(sum(lengths))))

    #print("\n".join([str(tup) for tup in len_info]))
    make_text_hist(lengths)
