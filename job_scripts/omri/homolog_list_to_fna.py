
import argparse
from Bio import SeqIO

import sys
import os

class Cluster(object):
    def __init__(self, name):
        self.name = name
        self.out_name = self.name + ".fna"

        self.contigs = {}

        self.overwrite = True

    def parse_cluster_faa(self, cluster_faa_f, min_len=100):
        """ Parses the cluster .faa file to get headers to get from the .fna file """

        with open(cluster_faa_f, 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                if len(record.seq) < min_len:
                    continue

                header = record.description
                
                # parse genome_name from the header
                # basically gets text between the first set of []
                genome_name = header.split("[", 1)[1].split("]", 1)[0]

                self.add_contig(genome_name, header)

    def add_contig(self, genome, contig):
        self.contigs[contig] = {"genome": genome, "written": 0}

    def check_cluster(self, genome_list):
        """ 
        Makes sure there is one and only one entry for each genomes on the list

        Returns True if so, False otherwise
        """
        
        genomes = [self.contigs[contig]["genome"] for contig in self.contigs]

        # check for duplicates
        if len(set(genomes)) < len(genomes):
            #print(("Duplicate!: ", genomes))
            return False

        # check that all genomes in the genome list are present
        for genome in genome_list:
            if genome not in genomes:
                #print(("Missing!: {} ".format(genome), genomes))
                return False

        return True


        self.fh = open(self.name + ".fna", 'a')

    def write_seq(self, seq):
        """ I open and close the FH on each write to avoid having too many file handles open """
        
        # overwrite on the first open
        if self.overwrite:
            mode = 'w'
            self.overwrite = False
        else:
            mode = 'a'

        # there were some seqs that were blank -- don't count those
        if seq.seq:
            self.contigs[seq.description]["written"] = 1

            with open(self.out_name, mode) as IN:
                SeqIO.write(seq, IN, "fasta")


    def check_all_written(self):
        """ Checks if all contigs have been written """

        values = [self.contigs[contig]["written"] for contig in self.contigs]
        if 0 not in values:
            return True
        else:
            return False



class Genome(object):
    def __init__(self, name, filename=None):
        self.name = name

        self.filename = filename

        self.contig2cluster = {}

    def add_contig(self, contig, cluster):
        # make sure the contig is not in any other clusters
        assert contig not in self.contig2cluster

        self.contig2cluster[contig] = cluster

    def extract_fna(self, homologues_dir):
        with open(os.path.join(homologues_dir, self.filename + ".fna"), 'r') as IN:
            for record in SeqIO.parse(IN, "fasta"):
                try:
                    cluster = self.contig2cluster[record.description]
                except KeyError:
                    continue

                cluster.write_seq(record)



# Functions to get the genomes we are looking for
def make_genome_objects(input_dir):
    """ Returns a dict that links genome name (to be parsed from GH files) to a Genome object """
    filename_to_genomes = parse_genome_name_from_gbks(input_dir)

    name_to_genome = {}
    for filename, genome_name in filename_to_genomes.items():
        genome = Genome(genome_name, filename)
        
        name_to_genome[genome_name] = genome

    return name_to_genome

def parse_genome_name_from_gbks(gbk_dir):
    """ Returns a dict of the file name: genome name """

    filename_to_genome_name = {}
    for path in os.listdir(gbk_dir):
        basename, ext = os.path.splitext(path)

        if ext == ".gbk":
            genome_name = _get_organism_from_gbk(os.path.join(gbk_dir, path))
            filename_to_genome_name[path] = genome_name

    return filename_to_genome_name
                           
def _get_organism_from_gbk(gbk_f):
            
    # open the .gbk file and parse the first source tag for organism
    with open(gbk_f, 'r') as IN:
        for record in SeqIO.parse(IN, "genbank"):
            for feat in record.features:
                if feat.type == "source":
                    return feat.qualifiers["organism"][0]
 

def parse_cluster_list(clusters_f):
    clusters = []
    with open(clusters_f, 'r') as IN:
        for line in IN:
            clusters.append(line[:-1])

    return clusters

def process_cluster_faa(intersect_dir, cluster_list, names_to_genomes):
    """ Processes all the cluster files and makes Cluster objects with info """
    clusters = []
    for clus_faa in cluster_list:
        name = os.path.splitext(clus_faa)[0]

        cluster = Cluster(name)
        cluster.parse_cluster_faa(os.path.join(intersect_dir, clus_faa))

        if cluster.check_cluster(names_to_genomes):
            # add the contig to the Genome object now that we know it is a good cluster
            for contig in cluster.contigs:
                genome = names_to_genomes[cluster.contigs[contig]["genome"]]
                genome.add_contig(contig, cluster)

            clusters.append(cluster)
            
    return clusters


def extract_fna_sequences(genomes, homologues_dir):
    for genome in genomes:
        print("Extracting sequences from {}...".format(genome.name))
        genome.extract_fna(homologues_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Loops up original nucleotide sequences for a list of clusters.")
    parser.add_argument("-clusters", help="file of clusters from get_homologues. The core list, for example")
    parser.add_argument("-intersect", help="the intersections directory that contains all the files in the list")
    parser.add_argument("-homologues_dir", help="the output directory from get_homologues")
    parser.add_argument("-input", help="the input directory of gbk file for get_homologues")
    
    args = parser.parse_args()


    """ The idea here is to make a map of contigs in a genome to the cluster they belong to so we can process each genome once and write each contig to a different filehandle corresponding to its cluster """

    # get a map from GH genome names to Genome objects
    names_to_genomes = make_genome_objects(args.input)


    # add cluster information to Genomes
    cluster_list = parse_cluster_list(args.clusters)
    clusters = process_cluster_faa(args.intersect, cluster_list, names_to_genomes)
    
    print("{} good clusters detected.".format(len(clusters)))

    genomes = names_to_genomes.values()

    # sanity check
    for genome in genomes:
        print("{} contigs (1 per cluster) in genome {}".format(len(genome.contig2cluster), genome.name))

    extract_fna_sequences(genomes, args.homologues_dir)


    # make sure all sequences were written to each cluster
    for cluster in clusters:
        if not cluster.check_all_written():
            print("Warning cluster {} is missing sequences! Deleting file.".format(cluster.name))
            try:
                os.remove(cluster.out_name)
            except:
                continue
