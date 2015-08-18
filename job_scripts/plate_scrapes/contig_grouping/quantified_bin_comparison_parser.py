
import argparse
import subprocess
import pandas
import sys

from mypyli import heatmap

class Comparison(object):
    FASTA_DIR = ""
    MUMMER_DIR = ""
    
    def __init__(self, qry, ref, gANI, AF):
        self.qry = qry
        self.ref = ref
        self.gANI = gANI
        self.AF = AF


    def make_dotplot(self, prefix):
        ref = self.FASTA_DIR + "/" + self.ref + ".fasta"
        qry = self.FASTA_DIR + "/" + self.qry + ".fasta"
        
        # run nucmer 
        command = [   
            self.MUMMER_DIR + "/" + "nucmer",
            "--prefix", prefix,
            ref, qry
        ]

        print("\nRunning:\n  {}".format(" ".join(command)))
        
        if subprocess.call(command):
            raise Exception("Running nucmer failed!")

        # filter the delta file
        delta = prefix + ".delta"
        filt_delta = delta.rsplit(".", 1)[0] + ".filt" + ".delta"
        
        command = [
            self.MUMMER_DIR + "/" + "delta-filter",
            delta,
            "-q",
            "-i", "70",     # identity requirement
            "-l", "1000",   # length requirement
        ]
    
        print("\nRunning:\n  {}".format(" ".join(command)))
        with open(filt_delta, 'w') as OUT:
        
            if subprocess.call(command, stdout=OUT):
                raise Exception("Running delta-filter failed!")

        # make the plot
        command = [
            self.MUMMER_DIR + "/" + "mummerplot",
            "--png",
            filt_delta,
            "-R", ref,
            "-Q", qry,
            "--layout",
            "--prefix", prefix
        ]

        print("\nRunning:\n  {}".format(" ".join(command)))
        
        if subprocess.call(command):
            raise Exception("Running mummerplot failed!")


def parse_matrix(input_f, skip_self=True):
    """ Generator for Comparison objects resulting from the matrix """
    with open(input_f, "r") as IN:
        headers = IN.readline()
        refs = headers[:-1].split("\t")[1:]

        for line in IN:
            genome, values = line[:-1].split("\t", 1)
            for value, ref in zip(values.split("\t"), refs):

                if skip_self:
                    if ref == genome:
                        continue

                # basically parses the text representation of a tuple
                # this removes the parens and the splits on ", "
                gANI, AF = value[1:-1].split(", ")
                yield Comparison(genome, ref, float(gANI), float(AF))

def make_dotplots(input_f, min_af=.60):
    for comp in parse_matrix(input_f):
        if comp.AF >= min_af:
            comp.make_dotplot(comp.qry + "__--__" + comp.ref)


def write_sorted_hits(input_f, output_f="sorted_comparisons.txt", sort_first="AF"):
    comparisons = [comparison for comparison in parse_matrix(input_f)]

    def sort_function(comparison):
        if sort_first == "AF":
            return (comparison.AF, comparison.gANI)
        else:
            return (comparison.gANI, comparison.AF)

    with open(output_f, 'w') as OUT:
        for comparison in sorted(comparisons, key=sort_function, reverse=True):
            OUT.write("{c.qry}__->__{c.ref}\t{c.gANI}\t{c.AF}\n".format(c=comparison))

def make_heatmap(input_f):
    tmp_dict = {}
    for comp in parse_matrix(input_f, skip_self=False):
        try:
            tmp_dict[comp.qry][comp.ref] = comp.AF
        except KeyError:
            tmp_dict[comp.qry] = {comp.ref: comp.AF}


    df = pandas.DataFrame.from_dict(data=tmp_dict, orient='index')
    #print(df)
    #sys.exit()


    hm = heatmap.Heatmap(df)
    hm.test()
    hm.plot_heatmap((.05, .15, .8, .7))
    hm.plot_legend((.05, .90, .2, .05))

"""
    hm = heatmap.HiearchicalHeatmap
    hm.frame = df
    hm.row_method = None
    hm.column_method = None
    fig, axm, axcb, cb = hm.plot(hm)
    fig.savefig("figure.png")
"""


def dereplicate(input_f, min_AF=.80):
    """ 
    Prints a list of dereplicated genomes; replicated is >= min_AF 
    
    """
    class Organism(object):
        def __init__(self):
            self.assemblies = set()
            self.comparisons = []

        def add_assembly(self, assembly):
            """ Good to add a single assembly for genomes without replications """
            self.assemblies.update({assembly})

        def add_comparison(self, comparison):
            """ Use to add a comparison; good for replicated organisms. """
            # add the comparison
            self.comparisons.append(comparison)

            # update the assemblies set
            self.assemblies.update({comparison.qry, comparison.ref})
   
        def pick_representative(self, min_AF=.8):
            """ Picks an assembly to be the representative for the organism. """
            # is not replicated, just return the seq
            if len(self.assemblies) == 1:
                return self.assemblies[0]

            # if replicated, pick the one with the lowest AF to the others (above some minimum) because that means
            # that thet genome has more content than the others. However, if we go too low, it is likely
            # the bin would just be a junk catch-all bin
            
            # tuple with the asm, AF of the lowest AF found so far
            lowest_AF = ("", 1)

            for comp in self.comparisons:
                if comp.AF < min_AF:
                    continue

                if comp.AF < lowest_AF[1]:
                    lowest_AF = (comp.qry, comp.AF)

            return lowest_AF[0] 
    
    organisms = []

    # this is a dict for all assemblies that will later be used to determine if that asm has been considered
    assemblies = {}    
    for comp in parse_matrix(input_f):
        assemblies[comp.qry] = 0
        if comp.AF >= min_AF:

            for organism in organisms:
                if comp.qry in organism.assemblies or comp.ref in organism.assemblies:
                    organism.add_comparison(comp)
                    break

            else:
                org = Organism()
                org.add_comparison(comp)
                organisms.append(org)

    dereplicated_asms = []
    # get best representative for replicated organisms
    for organism in organisms:
        # dissolve single-way hits
        if len(organism.comparisons) < 2:
            continue

        rep = organism.pick_representative()
        print("{} selected as representative for {}".format(rep, organism.assemblies), file=sys.stderr)
        
        # tell the dict that we have used this asm
        for asm in organism.assemblies:
            assemblies[asm] = 1

        dereplicated_asms.append(rep)

    
    # add all the non-replicated assemblies
    for asm, status in assemblies.items():
        if status == 0:
            dereplicated_asms.append(asm)

    
    print("\n".join(sorted(dereplicated_asms)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parses the matrix from calculate_ANI_and_AF.py and performs whatever functions are not commented out. This could be used by others but it is not set up for that currently.")
    parser.add_argument("-i", help="table of gANI and AF", required=True)
    parser.add_argument("-f", help="directory with all the fastas")
    parser.add_argument("-m", help="mummer directory")
    parser.add_argument("-min_af", help="minimum AF for mummerplotting", type=float)

    args = parser.parse_args()

    if args.f:
        Comparison.FASTA_DIR = args.f
    if args.m:
        Comparison.MUMMER_DIR = args.m

    #make_dotplots(args.i, args.min_af)
    write_sorted_hits(args.i)
    #make_heatmap(args.i)
    #dereplicate(args.i)
