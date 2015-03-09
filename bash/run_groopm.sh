#!/usr/bin/env bash

# Author: Hunter Cameron
# Date: 2/18/2015
# Updated:
#
# Description: A script to run the pipeline of GroopM to bin metagenomic reads.


set -eu -o pipefail      # exit at any non-0 return, unset variables, pipe fails

function HELP {
    echo "$(basename $0) -f <contigs.fasta> -b <index_bam_dir/> -d <groopm_db.gm> -o <output_dir_for_fasta_bins/>"
    exit 0
}

fasta=''
bam_dir=''
gm_database='groopm_db.gm'
output_dir=$PWD

while getopts "f:b:d:o:h" arg; do

    case $arg in
        f) # set fasta
            fasta=$OPTARG
            ;;

        b) # set bam dir
            bam_dir=$OPTARG
            ;;
        
        d) # set database name
            gm_database=$OPTARG
            ;;

        o) # set output dir
            output_dir=$OPTARG
            ;;

        h)  # show help
            HELP
            ;;
        *) # unrecognized
            echo -e "\nOption \"$OPTARG\" not allowed."
            HELP
            ;;
    esac
done



# check variables and make sure read input files exist
for var in fasta bam_dir gm_database output_dir; do
    if [ "${!var}" == '' ]; then
        echo -e "\n$var was not specified."
        HELP
    else
        echo "$var = ${!var}"
        case $var in 
            fasta)
                if [ ! -e ${!var} ]; then
                    echo -e "\nError - Path does not exist: ${!var}"
                fi
                ;;

            bam_dir)
                if [ ! -d ${!var} ]; then
                    echo -d "\nError - Not a directory: ${!var}"
                fi
                ;;
        esac

    fi
done

echo -e "fasta = $fasta"
echo -e "bam dir = $bam_dir"
echo -e "groopm database = $gm_database"
echo -e "output directory = $output_dir"



# make output directory if needed
if [ ! -d $output_dir ]; then
    echo "Making output directory structure\n   mkdir -p $output_dir"
    mkdir -p $output_dir
fi


## Begin calling groopm commands
parse="groopm parse $gm_database $fasta $bam_dir/*.bam"
core="groopm core $gm_database"
recruit="groopm recruit $gm_database"
extract="groopm extract -o $output_dir $gm_database $fasta"


for cmd in parse core recruit extract; do
    command="${!cmd}"

    echo "Running $command"
    $command

    if [ $? == 0 ]; then
        echo -e "!!!    Successfully completed!!!\n"
    fi
done

