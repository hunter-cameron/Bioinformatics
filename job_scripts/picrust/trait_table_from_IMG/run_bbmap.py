
import argparse
import subprocess
import tempfile
import os

def main(args):

    if args.i:
        if len(args.i) == 1:
            status = run_bbmap(args.i[0], args.ref, args.prefix)
        else:
            if args.i[0].endswith(".gz"):
                suffix = ".fastq.gz"
            else:
                suffix = ".fastq"


            handle, filename = tempfile.mkstemp(suffix=suffix, dir=os.getcwd())
            cmd = "cat " + " ".join(args.i) + " > " + filename

            # dangerous command -- runs directly as shell, could be injected
            subprocess.call(cmd, shell=True)

            status = run_bbmap(filename, args.ref, args.prefix)
            #os.unlink(filename)


        if status:
            raise Exception("BBMap Failed to Finish Successfully")

        parse_coverage(args.prefix + ".covstats.txt", args.prefix)


def run_bbmap(fq, ref, prefix):
    cmd = ["bbmap.sh", 
            "in={}".format(fq),
            "ref={}".format(ref), 
            "nodisk=t", 
            "ambiguous=random",
            "covstats={}.covstats.txt".format(prefix),
            "out={}.sam".format(prefix)
            ]
    status = subprocess.call(cmd)

    return status

def parse_coverage(covstats_f, prefix):
    """ Writes a coverage file """
    annot_cov = {}

    with open(covstats_f, 'r') as IN:
        # get rid of header
        IN.readline()

        for line in IN:
            scaffold, cov = line.split("\t")[:2]

            annot = scaffold.split("annot=")[-1]

            annot_cov[annot] = annot_cov.get(annot, 0) + float(cov)

    with open(prefix + ".annot_cov", 'w') as OUT:
        for annot in sorted(annot_cov):
            OUT.write("{}\t{}\n".format(annot, str(int(annot_cov[annot]))))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="runs BBMap on one or more FASTQ files. --- I no longer remember what this is. I'm adding it to GitHub as a legacy script. I think it may be used to get counts of features through coverage.")
    parser.add_argument("-i", help="one or more interleaved read sets", nargs="*")
    parser.add_argument("-ref", help="reference FASTA file", required=True)
    parser.add_argument("-prefix", help="out prefix for the SAM file and the Covstats", required=True)

    args = parser.parse_args()

    main(args)
