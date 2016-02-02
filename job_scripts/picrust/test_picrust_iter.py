
import subprocess
import sys
import argparse
import os
import time

import pandas 
import numpy as np

class GenomeProfile(object):
    def __init__(self, name):
        self.name = name
        self.kos = {}

    def add_ko(self, name, count):
        self.kos[name] = float(count)

    def process(self, trait_tab, tree, base_dir=None):

        # make a directory to hold the analysis
        if base_dir is None:
            base_dir = self.name

        if not os.path.isdir(base_dir):
            os.mkdir(base_dir)

        # store the expected path for the output for later
        self.predicted_path = base_dir + "/" + "ko_predicted.tab"
 
        if os.path.isfile(self.predicted_path):
            return

        else:

            # write the trait table excluding self
            self._write_all_but_self(trait_tab, base_dir + "/" + "edited_trait_table.tab")

            tree = os.path.abspath(tree)

            self.run_picrust(base_dir, tree, base_dir + "/" + "edited_trait_table.tab")


    def run_picrust(cls, base_dir, tree, trait_tab):
        """ Wraps all these commands into one bsub call to process """



        # find the path
        exe = subprocess.check_output(["which", "format_tree_and_trait_table.py"]).strip()
        format_files = "python {exe} -t {tree} -i {trait_tab} -o {base_dir}/format/".format(exe=exe, tree=tree, trait_tab=trait_tab, base_dir=base_dir)

        
        exe = subprocess.check_output(["which", "ancestral_state_reconstruction.py"]).strip()
        asr = "python {exe} -i {base_dir}/format/trait_table.tab -t {base_dir}/format/pruned_tree.newick -o {base_dir}/asr/ko_asr.tab".format(exe=exe, base_dir=base_dir)


        exe = subprocess.check_output(["which", "predict_traits.py"]).strip()
        predict = "python {exe} -i {base_dir}/format/trait_table.tab -t {base_dir}/format/reference_tree.newick -r {base_dir}/asr/ko_asr.tab -o {base_dir}/ko_predicted.tab -a".format(exe=exe, base_dir=base_dir)

    
        # link all the necessary commands into a single command
        super_command = "; ".join([format_files, asr, predict])

        subprocess.call([   "bsub", 
                            "-o", "{}/out".format(base_dir),
                            "-e", "{}/err".format(base_dir),
                            "-J", "picrust_test",
                            super_command
                            ])

    def compare_results(self, traits=None):
        """ Returns the correlation between observed and predicted and the NSTI value """

        predicted_kos = self._read_predicted_table()
   
        dict_for_pandas = {}      
        for ko in sorted(self.kos):
            if traits:
                if ko not in traits:
                    continue
            dict_for_pandas[ko] = {'observed': self.kos[ko], 'predicted': predicted_kos[ko]}


        df = pandas.DataFrame.from_dict(dict_for_pandas, orient="index")

        corr = df.corr()

        return corr.loc["observed", "predicted"], predicted_kos["metadata_NSTI"]
        
    def _write_all_but_self(self, trait_tab_f, out_f):
        with open(trait_tab_f, 'r') as IN, open(out_f, 'w') as OUT:
            for line in IN:
                if line.startswith(self.name):
                    continue
                else:
                    OUT.write(line)

    def _read_predicted_table(self):
        """ Reads the predicted table and returns a dict of KOs corresponding to this genome """

        with open(self.predicted_path, 'r') as IN: 
            headers = IN.readline()[:-1].split("\t")[1:]

            for line in IN:
                if line.startswith(self.name):
                    ko_dict = {}
                    for index, count in enumerate(line[:-1].split("\t")[1:]):
                        ko_dict[headers[index]] = float(count)

                    return ko_dict

            else:
                raise ValueError("Genome '{}' was not found in the predicted file.".format(self.name))
    

def parse_trait_table(trait_tab_f):
    """ Parses the trait table and returns a GenomeProfiles """

    profiles = []
    with open(trait_tab_f, 'r') as IN:
        headers = IN.readline()[:-1].split("\t")[1:]

        for line in IN:
            name, kos = line[:-1].split("\t", 1)
            
            gp = GenomeProfile(name)
           
            for indx, ko_count in enumerate(kos.split("\t")):
                gp.add_ko(headers[indx], ko_count)
    
            profiles.append(gp)

    return profiles



def _jobs_still_running():

    output = subprocess.check_output([   
                    "bjobs", 
                    "-J", "picrust_test"
                ])
    print(output)
    if output:
        return True
    else:
        return False

def trait_list_by_variance(trait_table):

    df = pandas.read_csv(trait_table, sep="\t", header=0, index_col=0)

    var_dict = {}

    for trait in df.columns:
        var = df[trait].var()

        var_dict[trait] = var

    return var_dict 

def main(args):

    output_dir = os.getcwd()

    profiles = parse_trait_table(args.i)

    for profile in profiles:
        profile.process(args.i, args.t, output_dir + "/" + profile.name)

    # wait for all jobs to finish
    while _jobs_still_running():
        time.sleep(10)

    """
    # get the 50% most variable traits
    var_traits = trait_list_by_variance(args.i)
    high_var = []
    for trait in sorted(var_traits, key=lambda trait: var_traits[trait], reverse=True):
        #print("{}\t{}".format(trait, str(var_traits[trait])))
        if var_traits[trait] > 1:
            high_var.append(trait)

    """
    with open("PICRUSt_test_results.tab", 'w') as OUT:

        OUT.write("\t".join(["name", "correlation", "NSTI"]) + "\n")
        for profile in profiles:
            corr, NSTI = profile.compare_results()
        
            OUT.write("\t".join([profile.name, str(corr), str(NSTI)]) + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="tests PICRUSt by dropping one genome from the table and then predicting it using the others")
    parser.add_argument("-t", help="newick tree of all the isolates")
    parser.add_argument("-i", help="trait table containing all the isolates")
    args = parser.parse_args()

    main(args)


