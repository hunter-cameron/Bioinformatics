

# Concoct recommends that you split contigs longer than 10kb (using their split_fasta script) to help mitigate assembly errors.
# Also suggests running bowtie2 for mapping with picard to remove duplicates
# This script does neither :)

# Requires bedtools 

set -eu -o pipefail      # exit at any non-0 return, unset variables, pipe fails

# add the GSL library
export LD_LIBRARY_PATH=/proj/dangl_lab/apps/CONCOCT-0.4.0/GSL/lib:$LD_LIBRARY_PATH


function HELP {
    echo "
    ===============================================================================
    |    run_concoct.sh   -   A pipeline to run CONCOCT to bin reads using BAM cov |
    ===============================================================================
    

    This script runs the CONCOCT pipeline on all the bam files (files that match *.bam) in the supplied directory. Also makes a sample names file based on the sample name prefix (for file CL_Col_mCL_r1.norm.bam;; sample name = CL_Col_mCL_r1). 

    USAGE:

    $(basename $0) -c CONCOCT <scripts dir> -r <reference_fasta> -b <BAM directory> -d <output_directory>
            

    "
    exit 1
}




scripts=/proj/dangl_lab/apps/CONCOCT-0.4.0/scripts      # path to concoct scripts
bam_dir=''      # directory with bam files
out=$(pwd)
fasta_f=''

while getopts ":c:b:d:r:h" arg; do

    case $arg in
        c) # set scripts dir
            scripts=$OPTARG
            ;;
        
        b) # set bam_dir
            bam_dir=$OPTARG
            ;;

        d) # set working directory
            working_dir=$OPTARG
            ;;

        r) # set ref_fasta
            ref_fasta=$OPTARG
            ;;

        ?) # unknowns
            echo -e "\nOption \"$OPTARG\" not allowed."
            HELP
            ;;
    esac
done


# check variables and make sure read input files exist
for var in scripts bam_dir ref_fasta out; do
    if [ "${!var}" == '' ]; then
        echo -e "\n$var was not specified."
        HELP
    else
        # check paths exist
        case $var in 
            scripts|bam_dir|ref_fasta)
                if [ ! -e ${!var} ]; then
                    echo -e "\nError - Path does not exist: ${!var}"
                fi
                ;;

            # check output dir
            out)
                if [ ! -d ${!var} ]; then
                    echo -e "\nMaking directory ${!var}"
                    mkdir -p ${!var}
                fi
                ;;
        esac

    fi
done





echo "Using reference fasta: $ref_fasta"
echo "Using BAM directory: $bam_dir"
echo "Using output directory; $out"
echo "Using CONCOCT scripts dir: $scripts"

sample_names_file=$out/sample_names.txt

cd $bam_dir

# write sample names file based on the file prefixes
> $sample_names_file
for file in *.bam; do
    prefix=${file/.*/}
    echo $prefix >> $sample_names_file
done

if [ ! -e $out/concoct_coverage_mean.csv ]; then
    # generate the input table for coverage
    python $scripts/gen_input_table.py --samplenames $sample_names_file $ref_fasta *.bam > $out/concoct_coverage_table.csv

    # trim coverage file to just mean cov
    cut -f1,3-26 $out/concoct_coverage_table.csv > $out/concoct_coverage_mean.csv

fi


# run concoct
if [ ! -e $out/clustering_gt1000.csv ]; then 
    concoct --coverage_file $out/concoct_coverage_mean.csv --composition_file $ref_fasta -c 1000 -b $out
fi


# I haven't found that anything below this actually does anything
exit 0 

# validate initial clustering
#if [ ! -e $out/clustering_eval_gt1000.csv ]; then
#    echo "Validating original clustering..."
#    perl $scripts/Validate.pl --cfile $out/clustering_gt1000.csv --ffile $ref_fasta --ofile $out/clustering_eval_gt1000.csv
#fi

# generate the input table for linkage
if [ ! -e $out/concoct_linkage_table.csv ]; then
    python $scripts/bam_to_linkage.py --regionlength 500 --fullsearch --samplenames $sample_names_file $ref_fasta *.bam > $out/concoct_linkage_table.csv
fi

# add linkage info to clustering
if [ ! -e $out/clustering_gt1000_link.csv ]; then
    perl $scripts/ClusterLinkNOverlap.pl --cfile $out/clustering_gt1000.csv --lfile $out/concoct_linkage_table.csv --covfile $out/concoct_coverage_mean.csv --ofile $out/clustering_gt1000_link.csv
fi

# removed because it requires species classification of the contigs
# validate new clustering
#if [ ! -e $out/clustering_eval_link.csv ]; then
#    echo "Validating linkage clustering..."
#    perl $scripts/Validate.pl --cfile $out/clustering_gt1000_link.csv --ffile $ref_fasta --ofile $out/clustering_eval_link.csv
#fi


echo "Successfully Completed."

exit 0
echo "

For Cluster Validation:
    N = number of contigs clustered
    M = number with labels
    S = number of unique labels
    K = number of clusters
    Rec = recall
    Prec = precision
    NMI RAND = normalized mutual information RAND score
    AdjRAND = adjusted RAND score
"
