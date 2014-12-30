

import sys
import os
import argparse

from Bio import SeqIO
import re
import subprocess
import tempfile
import time

def parse_query_from_command(command):
    match = re.search("-query ([^ ]+)", command)
    if match:
        query = match.groups(1)[0]
    else:
        raise AssertionError("Could not find -query [query]. Make sure -query is included in the BLAST command.")

    return query

def parse_out_from_command(command):
    match = re.search("-out ([^ ]+)", command)
    if match:
        return match.groups(1)[0]
    else:
        print("-out not found")
        return None

def split_fasta_file(fasta, bins, out):
    if not os.path.isdir(out):
        os.makedirs(out)

    # count seqs in fasta -- probably a faster way to do this (without calling to grep)
    count = 0
    with open(fasta, 'rU') as IN:
        for line in IN:
            if line.startswith(">"):
                count += 1

    bin_size = (count // bins) + 1
    
    # split the seqs
    print("Splitting FASTA file...", file=sys.stderr)
    files = []
    with open(fasta, 'rU') as IN:
        record_iter = SeqIO.parse(IN, "fasta")

        for indx, bin in enumerate(_bin_generator(record_iter, bin_size)):
            out_file = out + fasta.split("/")[-1] + "_" + str(indx)
            with open(out_file, 'w') as OUT:
                print("   Splitting to {}".format(out_file), file=sys.stderr)
                SeqIO.write(bin, OUT, 'fasta')
                files.append(out_file)

    print("", file=sys.stderr)
    return files

def _bin_generator(iterator, bin_size):
    bin = []
    for item in iterator:
        bin.append(item)
        if len(bin) >= bin_size:
            yield bin
            bin = []
    # return any left overs
    if bin:
        yield bin

def lsf_blast(fa_file, command, queue, out, rand_id):
    
    bsub_command = "bsub -o parallel_blast.out -e parallel_blast.err -J {} -q {}".format(rand_id, queue)

    blast_command = sub_blast_command(command=command, query=fa_file, out=fa_file + ".blastout")
    
    # replace and double quotes that may have been used with singles
    blast_command = blast_command.replace("'", "\"")

    subprocess.call(bsub_command.split(" ") + _str2subprocess(blast_command))

    print("   Submitted Job: {}".format(bsub_command + " " + blast_command), file=sys.stderr)
    return fa_file, fa_file + ".blastout"

def sub_blast_command(command, query, out):
    
    old_query = parse_query_from_command(command)
    old_out = parse_out_from_command(command)

    command = command.replace("-query {}".format(old_query), "-query {}".format(query))
    command = command.replace("-out {}".format(old_out), "-out {}".format(out))

    return command

def wait_until_finished(job_name):
    """ Checks if jobs are still running every 5 minutes """
    print("Waiting for jobs to run...", file=sys.stderr, end="\n\n")
    output = 1
    while output:
        time.sleep(3)
        output = subprocess.check_output(["bjobs", "-J", job_name], stderr=open("/dev/null", 'w'))
    return True

def concatenate_split_results(in_outs, output):
    print("Concatenating split results...", file=sys.stderr)
    
    with open(output, 'w') as OUT:
        for fa, result in in_outs:
            with open(result, 'r') as IN:
                for line in IN:
                    OUT.write(line)
    print("   Concatenated to: {}".format(output), file=sys.stderr, end="\n\n")

def check_blast_command(command, fasta, out):
    """ Checks if the BLAST command is capable of running using a single sequence"""

    print("Checking BLAST command...", file=sys.stderr, end="\n\n")
    # write a file with a single fasta sequence
    with open(fasta, 'r') as IN, open(out + "/" + "test.fa", 'w') as OUT:
        record_iter = SeqIO.parse(IN, "fasta")
        record = next(record_iter)
        SeqIO.write(record, OUT, "fasta")

    command = sub_blast_command(command, out + "/" + "test.fa", "/dev/null")
    try:
        subprocess.check_call(_str2subprocess(command), stderr=sys.stdout)
    except subprocess.CalledProcessError as e:
        print("\n\nChecking BLAST command returned an error. Make sure you can run the BLAST command from the command line.")
        raise e

def _str2subprocess(command):
    """ converts a string command to subprocess' argument array format """
        
    command = command.replace("\"", "'")
    elements = command.split(" ")
    final_e = []
    temp = ""
    cat = False
    for element in elements:

        if element.startswith("'"):
            cat = True
            element = element.replace("'", "")

        if cat:
            temp += (element + " ")
        else:
            final_e.append(element)

        if element.endswith("'"):
            temp = temp[:-2]
            final_e.append(temp)
            temp = ""
            cat = False

    return final_e

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parallelize the running of a BLAST search and return the results as if it was a serial search! All for the low price of $19.99!")
    parser.add_argument("-n", help="number of nodes to use", type=int, default=10)
    parser.add_argument("-c", help="blast command; hint: wrap this in quotes")
    parser.add_argument("-out", help="directory for the output", default="./")
    parser.add_argument("-q", help="the queue", default="week")
    args = parser.parse_args()

    query = parse_query_from_command(args.c)
   
    # make a temp dir and use it as a random id for the job names
    tmp_out = tempfile.mkdtemp(prefix="fasplit_", dir=os.path.abspath(args.out))
    rand_id = tmp_out.split("/")[-1]
    
    files = split_fasta_file(os.path.abspath(query), args.n, tmp_out + "/")

    check_blast_command(args.c, os.path.abspath(query), tmp_out)

    print("Submitting BLAST jobs...", file=sys.stderr)
    in_outs = []
    for fa_file in files:
        in_out_tup = lsf_blast(fa_file, args.c, args.q, tmp_out + "/", rand_id)
        in_outs.append(in_out_tup)

    print("", file=sys.stderr)

    wait_until_finished(rand_id)

    # get the path the user supplied in the command for the output
    output = parse_out_from_command(args.c)
    if output:
        output = os.path.abspath(output)
    else:
        output = os.path.abspath(args.out) + "full_blast_results.txt"


    concatenate_split_results(in_outs, output)


    print("Job Successfully Completed!")
