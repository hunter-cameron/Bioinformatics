#!/usr/bin/env bash

# Author: Hunter Cameron
# Date: 12/16/2014
# Updated:
#
# Future Updates:
#   - Allow multiple read files to be specified
#
#
#
#
# Description: A script to run bowtie2 in a reproducible way. Allows for single reads and paired end reads to both be supplied as arguments.


set -eu -o pipefail      # exit at any non-0 return, unset variables, pipe fails

function HELP {
echo -e "Count the number of sequences in all the fastq files in a directory or file.\n    Gzipped fastq also welcome.\n    Assumes 1 sequence per 4 lines.\n\nUSAGE\n$(basename $0) -i <input>"
    exit 0
}

input=''
while getopts ":i:h" arg; do

    case $arg in
        i) # set input
            input=$OPTARG
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
for var in input; do
    if [ "${!var}" == '' ]; then
        echo -e "\n$var was not specified."
        HELP
    else
        #echo "$var = ${!var}"
        case $var in 
            input)
                if [ ! -e ${!var} ]; then
                    echo -e "\nError - Path does not exist: ${!var}"
                fi
                ;;
        esac

    fi
done

echo -e "file\t# sequences"

if [ -d $input ]; then
    # loop through directory gleaning every file with .fastq .fastq.gz .fq .fq.gz
    files=$(find $input -maxdepth 1 -regextype posix-egrep -regex ".*\.(fastq|fq)")
    for file in $files; do
        echo -e "$(basename $file)\t$((`sed -n '$=' $file` / 4))"
    done

    gz_files=$(find $input -maxdepth 1 -regextype posix-egrep -regex ".*\.(fastq.gz|fq.gz)")

    for file in $gz_files; do
        echo -e "$(basename $file)\t$((`zcat $file | sed -n '$='` / 4))"
    done
else
    case $(basename $input) in
        *"gz")      # it is a .gz file
            echo -e "$(basename $input)\t$((`zcat $input | sed -n '$='` / 4))"
            ;;
        *)  # not gz
            echo -e "$(basename $input)\t$((`sed -n '$=' $input` / 4))"
            ;;
    esac
fi
