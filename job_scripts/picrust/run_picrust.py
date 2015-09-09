
import argparse
import os

from mypyli import picrust


def predict_traits(tree, ref_traits, type, base_dir):
    pe = picrust.PicrustExecuter()

    return pe.predict_traits_wf(tree, ref_traits, type=type, base_dir=base_dir)


def predict_metagenome(otu_table, copy_numbers, trait_table, base_dir):
    pe = picrust.PicrustExecuter()
    
    return pe.predict_metagenome(otu_table, copy_numbers, trait_table, base_dir)

def main(args):

    # set up the output dir if needed
    if not os.path.isdir(args.out):
        os.mkdir(args.out)

    if args.wf == "predict_traits" or args.wf == "both":
        
        jobs = []
        for index, path in enumerate([args.traits, args.marker_counts]):
            if path:
                if index == 0:
                    job_name, new_traits = predict_traits(args.tree, args.traits, "trait", args.out)
                    print("Storing trait predictions to: {}".format(new_traits))
                    jobs.append(job_name)
                    args.traits = new_traits
                else:
                    job_name, new_markers = predict_traits(args.tree, args.marker_counts, "marker", args.out)
                    print("Storing marker predictions to: {}".format(new_markers))
                    jobs.append(job_name)
                    args.marker_counts = new_markers

        for job in jobs:
            picrust.PicrustExecuter.wait_for_job(job)

        print("Trait prediction complete.")

    if args.wf == "predict_metagenome" or args.wf == "both":
        job_name, predicted_metagenome = predict_metagenome(args.otu_table, args.marker_counts, args.traits, args.out)

        picrust.PicrustExecuter.wait_for_job(job_name)

        print("Predicted metagenome stored as {}".format(predicted_metagenome))

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="runs PICRUSt")
    parser.add_argument("-wf", help="choice of workflow to run", choices=["predict_traits", "predict_metagenome", "both"], required=True)

    parser.add_argument("-tree", help="Newick format tree of all OTUs")
    parser.add_argument("-traits", help="trait table")
    parser.add_argument("-marker_counts", help="counts of marker genes")
    parser.add_argument("-otu_table", help="otu table to use for metagenome prediction")
    parser.add_argument("-out", help="directory for output", default=os.getcwd())
    args = parser.parse_args()

    main(args)
