
import argparse

from mypyli import utilities, picrust

"""
def parse_predicted_trait_table(old_trait_table, new_trait_table):
    "" Returns a dict of genomes and a dict of trait_descriptions ""
    
    genome_dict = {}
    trait_names = []
    trait_metadata = {}
    
    with open(old_trait_table, 'r') as IN, open(new_trait_table, 'w') as OUT:
        for indx, line in enumerate(IN):
            fields = line[:-1].split("\t")

            # check for comment lines
            if fields[0].startswith("#"):
                # check for the header line
                if fields[-1] == "metadata_NSTI":
                    trait_names = fields[1:-1]

                    # remove the comment from the first field
                    fields[0] = fields[0][1:]
                    OUT.write("\t".join(fields[:-1]) + "\n")
                continue


            elif fields[0].startswith("metadata"):
                trait_metadata[fields[0]] = fields[1:] 
                continue


            # process each genome
            try:
                if float(fields[-1]) == 0:
                    OUT.write("\t".join(fields[:-1]) + "\n")
                    genome_dict[fields[0]] = None
            except:
                print(("Exception", indx+1, line))

    # build a dict of trait metadata
    trait_dict = {}
    for index, name in enumerate(trait_names):
        trait_dict[name] = {field: trait_metadata[field][index] for field in trait_metadata}


    return genome_dict, trait_dict
"""

def parse_predicted_trait_table(old_trait_table, new_trait_table):
    """ Returns a dict of genomes and a dict of trait_descriptions """
    
    genome_dict = {}
    trait_names = []
    trait_metadata = {}

    ttm = picrust.TraitTableManager(old_trait_table)

    # get an ordered list of traits excluding any metadata fields
    ordered_traits = ttm.get_ordered_traits()
    ordered_traits = [trait for trait in ordered_traits if not trait.startswith("metadata")]

    with open(new_trait_table, 'w') as OUT:
        # write header line
        OUT.write("\t".join(["OTU_IDs"] + ordered_traits) + "\n")

        # process each entry in the precalculated trait table
        for entry in ttm:

            # check for the metadata fields and extract them separately
            if entry.name.startswith("metadata"):
                trait_metadata[entry.name] = entry.traits
            else:
                if float(entry.traits["metadata_NSTI"]) == 0:
                    # set entry in genome dict for later use
                    genome_dict[entry.name] = None

                    # write the entry to the new trait table
                    picrust.TraitTableManager.write_entry(entry, OUT, ordered_traits)

    # build a dict of trait metadata
    trait_dict = {}
    for index, name in enumerate(trait_names):
        trait_dict[name] = {field: trait_metadata[field][index] for field in trait_metadata}


    return genome_dict, trait_dict


def map_genomes_to_copynumbers(copy_number_f, genomes):
    """ Looks up copynumbers for genomes using a supplied copy_number_f """
    with open(copy_number_f, 'r') as IN:
        for line in IN:
            name, copynum = line[:-1].split("\t")

            if name in genomes:
                genomes[name] = copynum

def write_copynumber_file(genomes, out):
    with open(out, 'w') as OUT:
        OUT.write("OTU_IDS\t16S_rRNA_Count\n")
        for genome, copynumber in genomes.items():
            OUT.write("{}\t{}\n".format(genome, copynumber))

def write_trait_metadata(traits, out):
    headers = sorted(traits[list(traits.keys())[0]])
    with open(out, 'w') as OUT:
        OUT.write("\t".join(["Trait"] + headers) + "\n")
        for trait in traits:
            to_write = []
            for name in headers:
                to_write.append(traits[trait][name])

            OUT.write("\t".join([trait] + to_write) + "\n")

def write_new_marker_fasta(old_marker, genomes, out):
    """ Writes a fasta of only the selected genomes to make tree building easier """
    utilities.write_fasta_by_header(old_marker, headers=genomes.keys(), out=out)

def main(args):

    genomes, traits = parse_predicted_trait_table(args.input, args.prefix + ".traits.tab")

    map_genomes_to_copynumbers(args.copynumber, genomes)

    # write all of the outputs
    write_copynumber_file(genomes, args.prefix + ".marker_copynumbers.tab")
    write_trait_metadata(traits, args.prefix + ".trait_metadata.tab")
    write_new_marker_fasta(args.fasta, genomes, args.prefix + ".markers.fasta")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extracts a trait table, copynumber file, sequence file, and trait metadata file from a precalcalculated file. Useful for going from a predicted trait table to genomes that have data to use as a base for adding new OTUs. It works by looking for NSTI (nearest sequenced taxon index) of 0 to indicate that it was sequenced.")

    parser.add_argument("-input", help="input predicted trait table")
    parser.add_argument("-copynumber", help="tab delimited copynumber file")
    parser.add_argument("-fasta", help="fasta format file of the marker gene")
    parser.add_argument("-prefix", help="prefix for output files")

    args = parser.parse_args()

    main(args)
