
import sys

suffix_list = ['raw.fastq.gz', 'raw.trimmo.singles.fastq', 'raw.trimmo.interleaved.fastq', 'trimmo.khmer_singles.fastq', 'trimmo.khmer_interleaved.fastq']
out_f = 'sample_read_counts.txt'
count_dict = {}
for ct_file in sys.argv[1:]:
    with open(ct_file, 'r') as IN:
        for line in IN:
            try:
                fasta, count = line[:-1].split("\t")
                prefix, suffix = fasta.split(".", 1)
            except:
                continue

            if suffix in suffix_list:
                try:
                    count_dict[prefix][suffix] = count
                except KeyError:
                    count_dict[prefix] = {}
                    count_dict[prefix][suffix] = count

with open(out_f, 'w') as OUT:
    OUT.write("\t".join(['sample'] + suffix_list) + "\n")
    for fasta in count_dict:
        OUT.write("\t".join([fasta] + [count_dict[fasta][suffix] for suffix in suffix_list]) + "\n")


