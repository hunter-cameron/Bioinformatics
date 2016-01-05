
import sys
import os
import subprocess
import logging
import time
import re

import pandas

logging.basicConfig()
LOG = logging.getLogger(__name__)

class TraitTableEntry(object):
    """ A single entry in a trait table """

    def __init__(self, name):
        self.name = name
        self.traits = {}

    def __str__(self):
        return "TraitTableEntry {}".format(self.name)

    def add_trait(self, trait, value):
        """ Checks traits to make sure it doesn't already exist in the dict and adds it """
        if trait in self.traits:
            raise ValueError("{} already has a trait called '{}'.".format(str(self), trait))
        else:
            # see if we can convert the trait into a number
            try:
                value = float(value)
            except ValueError:
                pass

            self.traits[trait] = value

    def correlation(self, other, traits=None):
        """ 
        Finds the correlation between self and other for the listed traits 
        
        If traits is not suppiled, uses all the traits from self.

        Only uses traits that both have. I'm not sure if this is the intuitive default behavior.
        It might make sense to throw an error if a trait isn't found.
        Or to return a list of the ultimate traits used?

        It may be unintuitive because if I pass the same set of traits to multiple comparisons,
        each comparison may actually use different traits and not tell me.
        """
        
        if traits is None:
            traits = list(self.traits.keys())

        pandas_dict = {}
        for trait in traits:
            try:
                st = self.traits[trait]
            except KeyError:
                continue

            try:
                ot = other.traits[trait]
            except KeyError:
                continue
     

            pandas_dict[trait] = {"self": st, "other": ot}

        if not pandas_dict:
            raise ValueError("No traits were shared between the entries.")

        
        df = pandas.DataFrame.from_dict(pandas_dict, orient="index")

        corr = df.corr(method="spearman")

        return corr.loc["self", "other"]
    

class TraitTableManager(object):
    """ A class for parsing and manipulating trait tables """

    def __init__(self, trait_table_f):
        self.trait_table_f = trait_table_f

        # get headers
        with open(self.trait_table_f, 'r') as IN:
            headers = IN.readline().rstrip().split("\t")

            # set the entry header replacing a comment line if present
            self.entry_header = headers[0].replace("#", "")
            self.traits = headers[1:]

    def __iter__(self):
        """ Yields a TraitTableEntry for each line in a trait table """

        with open(self.trait_table_f, 'r') as IN:
            # skip header line
            IN.readline()

            for line in IN:
                # skip blank lines
                if not line:
                    continue

                try:
                    name, trait_values = line.rstrip().split("\t", 1)
                except ValueError:
                    print((line,))

                tte = TraitTableEntry(name)

                for index, val in enumerate(trait_values.split("\t")):
                    tte.add_trait(self.traits[index], val)

                yield tte

    def get_ordered_traits(self, metadata_last=True):
        """ Returns an ordered list of traits by a natural sort algorithm that optionally sends metadata to the back. """

        def convert(char):
            """ Attempts to convert a character into an integer """
            try:
                return int(char)
            except ValueError:
                return char

        def nat_sort(entry):
            """ Performs a natural sort that will sort text and numbers in a way that makes sense """
            if metadata_last:
                if entry.startswith("metadata"):
                    # I append "~" to the beginning because of its high Unicode value
                    entry = "~~~" + "metadata"
            
            return [convert(char) for char in entry]
            

        return sorted(self.traits, key=nat_sort)

    def get_subset(self, subset_names, remove=False):
        """ 
        A filter around the iter method that only gets entries in the subset_names list (or removes them) 


        Something is wrong with this method and it sometimes returns incorrect results!!!
        """
        
        to_find = len(subset_names)
        found = 0
        for entry in self:

            # check if we have exhausted the list to speed up the search
            if found == to_find:
                if remove:
                    yield entry
                else:
                    return

            if entry.name in subset_names:

                if not remove:
                    found += 1
                    yield entry
            else:
                if remove:
                    found += 1
                    yield entry



    @staticmethod
    def write_entry(entry, fh, traits):
        to_write = [entry.name]

        for trait in traits:
            try:
                to_write.append(str(entry.traits[trait]))
            except KeyError:
                LOG.warning("Entry {} doesn't have trait {}. Writting 'NA'".format(str(entry), trait))
                to_write.append("NA")

        fh.write("\t".join(to_write) + "\n")


class PicrustExecuter(object):
    """ Runs PICRUSt """

    job_id = 0

    @classmethod
    def predict_traits_wf(cls, tree, trait_table, type="trait", base_dir=None):
        """ Runs the predict_traits_wf. Returns a name and an output path """
        # make a directory to hold the analysis
        if base_dir is None:
            base_dir = os.getcwd() + "/" + "picrust_project"
        else:
            base_dir = os.path.abspath(base_dir)

        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)

        if type == "trait":
            format_dir = base_dir + "/" + "format_trait"
            predict_out = base_dir + "/" + "predicted_traits.tab"
        elif type == "marker":
            format_dir = base_dir + "/" + "format_marker"
            predict_out = base_dir + "/" + "predicted_markers.tab"
        else:
            raise ValueError("type must be one of 'trait', 'marker'")
        
        format_cmd = cls._get_format_command(trait_table, tree, format_dir)

        # formatted paths
        fmt_table = format_dir + "/" + "trait_table.tab"
        fmt_tree = format_dir + "/" + "reference_tree.newick"
        prun_tree = format_dir + "/" + "pruned_tree.newick"
        
        
        asr_out = format_dir + "/" + "asr.tab"
        reconstruct_cmd = cls._get_asr_command(fmt_table, prun_tree, asr_out)

        predict_cmd = cls._get_predict_traits_command(fmt_table, asr_out, fmt_tree, predict_out)


        # link all the necessary commands into a single command
        super_command = "; ".join([format_cmd, reconstruct_cmd, predict_cmd])

        job_name = "picrust_cmd{}".format(cls.job_id)
        subprocess.call([   "bsub",
                            "-o", "{}/auto_picrust.out".format(base_dir),
                            "-e", "{}/auto_picrust.err".format(base_dir),
                            "-J", job_name,
                            super_command
                            ])
        cls.job_id += 1

        return job_name, predict_out

    def predict_metagenome(cls, otu_table, copy_numbers, trait_table, base_dir=None):
        # make a directory to hold the analysis
        if base_dir is None:
            base_dir = os.getcwd() + "/" + "picrust_project"
        else:
            base_dir = os.path.abspath(base_dir)

        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)

        norm_out = base_dir + "/" + "normalized_OTU_table.biom" 
        norm_cmd = cls._get_normalize_command(otu_table, copy_numbers, norm_out)

        predict_out = base_dir + "/" + "predicted_metagenome.tab"
        predict_cmd = cls._get_predict_metagenome_command(norm_out, trait_table, out=predict_out)

        # link all the necessary commands into a single command
        super_command = "; ".join([norm_cmd, predict_cmd])

        job_name = "picrust_cmd{}".format(cls.job_id)
        subprocess.call([   "bsub",
                            "-o", "{}/auto_picrust.out".format(base_dir),
                            "-e", "{}/auto_picrust.err".format(base_dir),
                            "-J", job_name,
                            super_command
                            ])
        cls.job_id += 1

        return job_name, predict_out

    @classmethod
    def wait_for_job(cls, job_name="picrust_cmd*"):
        """ waits for job to complete, checks every 10 seconds """
        while cls._job_running(job_name):
            time.sleep(10)

    @staticmethod
    def _job_running(job_name="picrust_cmd*"):

        output = subprocess.check_output([
                    "bjobs",
                    "-J", "picrust_cmd*"
                ])
        #print(output)
        if output:
            return True
        else:
            return False

    @staticmethod
    def _get_format_command(trait_tab, tree, out):
        exe = subprocess.check_output(["which", "format_tree_and_trait_table.py"]).strip()
        format_files = "python {exe} -t {tree} -i {trait_tab} -o {out}".format(exe=exe, tree=tree, trait_tab=trait_tab, out=out)

        return format_files

    @staticmethod
    def _get_asr_command(trait_table, tree, out):
        exe = subprocess.check_output(["which", "ancestral_state_reconstruction.py"]).strip()
        asr = "python {exe} -i {trait_table} -t {tree} -o {out}".format(exe=exe, trait_table=trait_table, tree=tree, out=out)

        return asr

    @staticmethod
    def _get_predict_traits_command(trait_table, asr_table, tree, out):
        exe = subprocess.check_output(["which", "predict_traits.py"]).strip()
        predict = "python {exe} -i {trait_table} -t {tree} -r {asr_table} -o {out} -a".format(exe=exe, trait_table=trait_table, asr_table=asr_table, tree=tree, out=out)

        return predict

    @staticmethod
    def _get_normalize_command(otu_table, copy_numbers, out):
        exe = subprocess.check_output(["which", "normalize_by_copy_number.py"]).strip()
        normalize = "python {exe} -i {otu_table} -c {copy_numbers} -o {out}".format(exe=exe, otu_table=otu_table, copy_numbers=copy_numbers, out=out)

        return normalize

    @staticmethod
    def _get_predict_metagenome_command(otu_table, trait_table, out):
        exe = subprocess.check_output(["which", "predict_metagenomes.py"]).strip()
        predict = "python {exe} -i {otu_table} -c {trait_table} -o {out} -f".format(exe=exe, otu_table=otu_table, trait_table=trait_table, out=out)
       
        return predict


def get_ko_by_function(ko_metadata_f, level=2):
    """ 
    Level 1 is the top level. 
    Level 2 is an intermediate level (Corresponding approx to COGs)
    Level 3 is the pathway level
    """
    if level not in [1, 2, 3]:
        raise ValueError("Level must be 1, 2, or 3.")

    data = {}
    with open(ko_metadata_f, 'r') as IN:
        # skip header line
        IN.readline()

        for line in IN:
            fields = line.rstrip().split("\t")

            ko_name = fields[0]
            ko_pathways = fields[-1]
            
            # multiple pathways sep by "|"
            for pathway in ko_pathways.split("|"):
                levels = pathway.split(";")
                
                try:
                    data[";".join(levels[:level])].append(ko_name)
                except KeyError:
                    data[";".join(levels[:level])] = [ko_name]
                except IndexError:
                    LOG.warning("{} did not have a pathway at the requested level.".format(ko_name))

    return data

def get_plant_associated_kos(plant_associated_f):
    """ Reads in a database of plant associated kos; returns a dict of lineage: KOs """

    data = {}
    with open(plant_associated_f, 'r') as IN:
        # skip header
        IN.readline()
        
        for line in IN:
            lineage, kos = line[:-1].split("\t")
            
            if kos:
                data[lineage] = kos.split(";")
            else:
                data[lineage] = []


    return data

if __name__ == "__main__":

    args = sys.argv[1:]

    test_executer(args) 
