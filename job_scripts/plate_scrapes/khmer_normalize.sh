


# k = 20 recommended for filtering; k = 32 recommended for partitioning
# though k = 32 is the default in the filtering program (program is newer than recommedation)

# normalize C=20 -> filter C=100;; normalize C=10 -> filter C=50  for metagenomes

# keep n=4 and vary x as needed "x=15e9 = 60G more than ever needed"
#
#  --- FROM THE KHMER DOCS--- >> These numbers are in bits... Khmer now expects the number as bytes
#For digital normalization, we recommend:
#
#        -x 2e9 for any amount of sequencing for a single microbial genome, MDA-amplified or single colony.
#        -x 4e9 for up to a billion mRNAseq reads from any organism. Past that, increase it.
#        -x 8e9 for most eukaryotic genome samples.
#        -x 8e9 will also handle most “simple” metagenomic samples (HMP on down)
#        For metagenomic samples that are more complex, such as soil or marine, start as high as possible. For example, we are using -x 64e9 for ~300 Gbp of soil reads.



set -eu -o pipefail      # exit at any non-0 return, unset variables, pipe fails

function HELP {
    echo "
===============================================================================
|    khmer_normalize.sh - A pipeline to normalize filtered metagenomic reads. |
===============================================================================
 

    $(basename $0) -p <paired_reads> -s <single_reads> -d <working_directory>
    

"
    exit 1
}


paired_read_f=''
single_read_f=''
working_dir=''

while getopts ":p:s:d:h" arg; do

    case $arg in
        p) # set paired file
            paired_read_f=$OPTARG
            ;;
        
        s) # set single file
            single_read_f=$OPTARG
            ;;

        d) # set working directory
            working_dir=$OPTARG
            ;;

        ?) # handle unknowns
            echo -e "\nOption \"$OPTARG\" not allowed."
            HELP
    esac
done


# check variables and make sure read input files exist
for var in paired_read_f single_read_f working_dir; do
    if [ "${!var}" == '' ]; then
        echo -e "\n$var was not specified."
        HELP
    else
        echo "$var = ${!var}"
        case $var in 
            paired_read_f|single_read_f)
                if [ ! -e ${!var} ]; then
                    echo -e "\nError - Path does not exist: ${!var}"
                fi
                ;;

            working_dir)
                if [ ! -d ${!var} ]; then
                    echo -e "\nMaking directory ${!var}"
                    mkdir -p ${!var}
                fi
                ;;
        esac

    fi
done


# want to do this with links in working directories because some of the output files are just adding extensions on filenames

# make links
base_paired_f=$(basename $paired_read_f)
base_single_f=$(basename $single_read_f)

# line rel2abs path -- will overwite existing links
ln -sf $(cd $(dirname $paired_read_f); pwd)/$(basename $paired_read_f) $working_dir/$base_paired_f
ln -sf $(cd $(dirname $single_read_f); pwd)/$(basename $single_read_f) $working_dir/$base_single_f

# get the prefix for the new reads
prefix=$(basename $paired_read_f)
prefix=${prefix/.*/}
echo $prefix

cd $working_dir

# check if already completed
if [ -e $prefix.khmer_interleaved.fastq.gz ] && [ -e $prefix.khmer_singles.fastq.gz ]; then
    echo -e "\n\nFound final files in outdir $working_dir. Delete the final files of the entire directory to re-run analysis."
    exit 0
fi



# normalize to median coverage 20 at k=20
if [ ! -e k20_C20.kh ]; then
    echo -e "\n\nNormalizing by median..."
    normalize-by-median.py -p --ksize 20 -C 20 --n_tables 4 -x 8e9 --savetable k20_C20.kh $base_paired_f 1>&2 
    normalize-by-median.py --loadtable k20_C20.kh --savetable k20_C20.kh $base_single_f 1>&2

else
    echo -e "\n\nFound .kh file...skipping"

fi

if [ ! -e $base_paired_f.keep.abundfilt ]; then
    echo "Filtering abundance..."
    filter-abund.py --variable-coverage --normalize-to 20 k20_C20.kh  $base_paired_f.keep $base_single_f.keep 1>&2
else 
    echo "Found abundfilt files...skipping"
fi


# extract reads that are still paired after filtering
if [ ! -e $base_paired_f.keep.abundfilt.se ]; then
    echo "Extracting paired reads..."
    extract-paired-reads.py $base_paired_f.keep.abundfilt 1>&2
else
    echo "Found extracted paired reads...skipping"
fi


# concatenate and zip the single reads      -- removed the gzip from each 
if [ ! -e $prefix.khmer_singles.fastq.gz ]; then
    echo "Concatenating single reads..."
    #cat *single.fastq.keep.abundfilt *keep.abundfilt.se | gzip -c > $prefix.khmer_singles.fastq.gz
    cat *single.fastq.keep.abundfilt *keep.abundfilt.se > $prefix.khmer_singles.fastq.gz
else
    echo "Found zipped single reads...skipping"
fi

# rename and zip paired reads
if [ ! -e $prefix.khmer_interleaved.fastq.gz ]; then
    echo "Zipping paired reads..."
    mv *interleaved*keep.abundfilt.pe $prefix.khmer_interleaved.fastq
    #gzip $prefix.khmer_interleaved.fastq
else
    echo "Found zipped paired reads...skipping"
fi



# remove intermediate files
rm *.abundfilt
rm *.keep
rm *.abundfilt.se
rm *.kh
rm $base_paired_f
rm $base_single_f

echo "Successfully Completed"
