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
    echo "$(basename $0) -i <index_prefix> -f <fwd_reads> -r <rev_reads> [-s <single_reads>] -o <output_sam -t <threads>"
    exit 0
}

indx=''
fwd=''
rev=''
single=''
out=''
threads=1
while getopts ":i:f:r:s:o:t:h" arg; do

    case $arg in
        i) # set index prefix
            indx=$OPTARG
            ;;
        f) # set forward reads
            fwd=$OPTARG
            ;;
        r) # set reverse reads
            rev=$OPTARG
            ;;
        s) # set single reads
            single=$OPTARG
            ;;

        o) # set output sam
            out=$OPTARG
            ;;

        t) # set number threads
            threads=$OPTARG
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
for var in indx fwd rev out threads; do
    if [ "${!var}" == '' ]; then
        echo -e "\n$var was not specified."
        HELP
    else
        echo "$var = ${!var}"
        case $var in 
            fwd|rev|single)
                if [ ! -e ${!var} ]; then
                    echo -e "\nError - Path does not exist: ${!var}"
                fi
                ;;
        esac

    fi
done

echo "singles = $single"



# get output directory
name="$(basename $out)"
dir=${out/$name/}

# make the directory structure if needed
if [ ! -d $dir ]; then
    echo "Making output directory structure\n   mkdir -p $dir"
    mkdir -p $dir
fi

# map the reads, option for paried only or paired and single
if [ single == '' ]; then

    bowtie2 --end-to-end --sensitive --phred33 --threads $threads --seed 20140418 -x $indx -1 $fwd -2 $rev -S $out

else
    
    bowtie2 --end-to-end --sensitive --phred33 --threads $threads --seed 20140418 -x $indx -1 $fwd -2 $rev -S $dir/tmp_sam_pe.sam


    bowtie2 --end-to-end --sensitive --phred33 --threads $threads --seed 20140418 -x $indx -U $single -S $dir/tmp_sam_se.sam

    # merge the sam files
    cat $dir/tmp_sam_pe.sam <(grep -v '^@' $dir/tmp_sam_se.sam) > $out

    # clean up temp files
    # rm $dir/tmp_sam_pe.sam $dir/tmp_sam_se.sam

fi
