
import argparse
import subprocess
import sys

def count_kmers(fasta_f, k, prefix, threads=1, hash_size="10M"):
    """ Uses Jellyfish to count k-mers and returns a path to the generated index"""

    # first run count
    cmd = "jellyfish count --output={output} --threads={threads} --mer-len={k} --size={hash_size} {fasta_f}".format(output=prefix + ".jf", threads=threads, k=k, hash_size=hash_size, fasta_f=fasta_f)

    print(cmd)
    code = subprocess.call(cmd.split(" "))

    if code:
        raise ValueError("Running Jellyfish failed with error code: {}".format(code))

    # check for multiple files and merge if found
    # so far there hasn't been any so I'm not going to do this...

    return prefix + ".jf"
    
def dump_kmer_counts(jellyfish_f, output_f, min_count=0, max_count=999999999999):
    """ Dumps a tab separated file of counts for each kmer """

    cmd = "jellyfish dump -L {min_count} -U {max_count} -o {output_f} -c -t {jellyfish_f}".format(min_count=min_count, max_count=max_count, output_f=output_f, jellyfish_f=jellyfish_f)

    code = subprocess.call(cmd.split(" "))

    if code:
        raise ValueError("Dumping Jellyfish results failed with error code: {}".format(code))

    return output_f


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-fasta", help="fasta file to count k-mers for", required=True)
    parser.add_argument("-k", help="length of k", required=True, type=int)
    parser.add_argument("-prefix", help="prefix for output files", default="jellyfish_counts")
    args = parser.parse_args()

    jf_file = count_kmers(args.fasta, args.k, args.prefix)
    counts_f = dump_kmer_counts(jf_file, args.prefix + ".kmer_counts.txt")

    print("Wrote counts to {}".format(counts_f))

    kmer_counts = {}
    with open(counts_f, 'r') as IN:
        for line in IN:
            kmer, count = line[:-1].split("\t")

            kmer_counts[kmer] = count

    print(len(kmer_counts))
    print(sys.getsizeof(kmer_counts))
