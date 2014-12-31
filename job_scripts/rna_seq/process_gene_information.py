
import argparse


"""A script that does all sorts of processing for the genes file with various filtering options 
       
        MODES:

        make_matrix - makes a matrix of gene counts per genome, per sample
        get_gene_counts - makes a list of all the genes sorted by the number of reads for that gene
"""

class FileReader(object):
    
    def __init__(self, file_path, qc_params):
        self.file_path = file_path
        self.qc_params = qc_params

    def get_next_line(self):
        with open(self.file_path, 'r') as IN:
            # skip header and then cycle through
            IN.readline()
            for line in IN:
                elements = line[:-1].split("\t")
                params = {
                        'sample': elements[0],
                        'genome': elements[1],
                        'contig': elements[2],
                        'read': elements[3],
                        'mapq':elements[4],
                        'locus_tag': elements[5],
                        'start': elements[6],
                        'end': elements[7],
                        'product':elements[8]
                        }
                line = Line(params)

                if self.line_qc(line):
                    yield line
                else:
                    continue

    def line_qc(self, line):
        """ Returns a True or a False depending on whether or not the read passes QC standards """
        
        # ignore things that didn't fall on a gene
        if "Nongene" in line.locus_tag:
            return False
           
        # check if ignore tRNA is enabled
        if self.qc_params['ignore_tRNA']:
            if "tRNA" in line.product:
                return False


        # after all checks
        return True


class Line(object):
    
    
    def __init__(self, params):
        self.sample = params['sample']
        self.genome = params['genome']
        self.contig = params['contig']
        self.read = params['read']
        self.mapq = params['mapq']
        self.locus_tag = params['locus_tag']
        self.start = params['start']
        self.end = params['end']
        self.product = params['product']



#
##
### Mode functions
##
#

def make_matrix(file_reader):
    """ Makes a matrix with the number of times each gene is found in each genome for each sample. """
    data = {}
    samples = {}
    for line in file_reader.get_next_line():
        if line.genome not in data:
            data[line.genome] = {}

        if line.product not in data[line.genome]:
            data[line.genome][line.product] = {}

        data[line.genome][line.product][line.sample] = data[line.genome][line.product].get(line.sample, 0) + 1
        samples[line.sample] = 1

    # print data
    samples = sorted(samples.keys())
    with open("gene_matrix.txt", 'w') as OUT:
        OUT.write("\t".join(["Genome", "Gene"] + samples) + "\n")
        for genome in sorted(data.keys()):
            for product in sorted(data[genome].keys()):
                OUT.write(genome + "\t" + product)
                for sample in samples:
                    count =  data[genome][product].get(sample, 0)

                    OUT.write("\t" + str(count))
                OUT.write("\n")

def get_genes_per_all(file_reader):
    """ Makes a matrix with the number of times each gene is found in each contig """
    data = {}
    for line in file_reader.get_next_line():
        # add product if needed
        if line.product not in data:
            data[line.product] = {"count": 0, "g_locs": [], "c_locs": []}
            
        # add locationif needed
        if line.genome not in data[line.product]["g_locs"]:
            data[line.product]["g_locs"].append(line.genome)
        if line.contig not in data[line.product]["c_locs"]:
            data[line.product]["c_locs"].append(line.contig)

        data[line.product]["count"] += 1

    # print data
    with open("genome_counts.txt", 'w') as GEN, open("contig_counts.txt", 'w') as CON:
        GEN.write("\t".join(["Gene", "Total Count", "Num Genomes", "Per Genome"]) + "\n")
        CON.write("\t".join(["Gene", "Total Count", "Num Contigs", "Per Contig"]) + "\n")

        for product, info in data.items():
            GEN.write("\t".join([product, str(info["count"]), str(len(info["g_locs"])), str(info["count"] / len(info["g_locs"]))]) + "\n")
            CON.write("\t".join([product, str(info["count"]), str(len(info["c_locs"])), str(info["count"] / len(info["c_locs"]))]) + "\n")


def make_subSeq_matrix(file_reader):
    """ Makes a matrix for the R package subSeq, This matrix needs to be transposed before it can be used. """
    data = {}
    samples = {}
    for line in file_reader.get_next_line():
        if line.genome not in data:
            data[line.genome] = {}

        if line.locus_tag not in data[line.genome]:
            data[line.genome][line.locus_tag] = {}

        data[line.genome][line.locus_tag][line.sample] = data[line.genome][line.locus_tag].get(line.sample, 0) + 1
        samples[line.sample] = 1

    # print data
    samples = sorted(samples.keys())
    with open("subSeq_gene_matrix.txt", 'w') as OUT:
        OUT.write("\t".join(["Gene"] + samples) + "\n")
        for genome in sorted(data.keys()):
            for tag in sorted(data[genome].keys()):
                OUT.write(genome + tag)
                for sample in samples:
                    count =  data[genome][tag].get(sample, 0)

                    OUT.write("\t" + str(count))
                OUT.write("\n")

def get_gene_counts(file_reader):
    """ Counts the number of times each gene occurs in the file """

    counts= {}
    for line in file_reader.get_next_line():
        counts[line.product] = counts.get(line.product, 0) + 1

    with open("gene_counts.txt", 'w') as OUT:
        for key in sorted(counts, key=counts.get, reverse=True):
            OUT.write(key + "\t" + str(counts[key]) + "\n")





if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-file", help="large gene file output from earlier script", required=True)
    parser.add_argument("-mode", choices=["make_matrix", "make_subSeq_matrix", "get_genes_per_all", "get_gene_counts"], required=True, help= "mode in which the program is to be ran")
    parser.add_argument("--ignore_tRNA", help="ignore genes that have a tRNA product", action='store_true')
    parser.add_argument("-genomes_file", help="file of genomes to include (default includes all)")
    parser.add_argument("--exclude", help="exclude the genomes in the -genomes_file", action="store_true")
    args = parser.parse_args()


    # set the Quality Control parameters
    qc_params = {}


    if args.genomes_file:
        genomes = []
        with open(args.genomes_file) as IN:
            for line in IN:
                genomes.append(line[:-1])

        if args.exclude:
            qc_params['exclude'] = genomes
        else:
            qc_params['include'] = genomes



    qc_params['ignore_tRNA'] = args.ignore_tRNA
    

    reader = FileReader(args.file, qc_params)

    if args.mode == "make_matrix":
        make_matrix(reader)
    
    elif args.mode == "get_gene_counts":
        get_gene_counts(reader)

    elif args.mode == "make_subSeq_matrix":
        make_subSeq_matrix(reader)
    
    elif args.mode == "get_genes_per_all":
        get_genes_per_all(reader)


