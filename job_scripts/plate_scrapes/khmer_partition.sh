


set -eu -o pipefail      # exit at any non-0 return, unset variables, pipe fails

function HELP {
    echo "
===============================================================================
|    khmer_partition.sh - A pipeline to partition filtered metagenomic reads. |
===============================================================================
 

$(basename $0) -k <.kh_file> -s <single_reads> -d <working_directory>
    


This script will dump several *.khmer_partition.* files. 

One of them may be large due to connections in Illumina reads. There is a way to deal with this in the doc/user of khmer's github.

At any rate, the partitions can be assembled individually though Velvet doesn't like the partition headers. 

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
                    echo "\nMaking directory ${!var}"
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

### Workflow from khmer's github page in user
if [ ! -e $prefix.ct ]; then
    load-graph.py -k 32 -N 4 -x 16e9 $prefix $fq.gz
fi

# partition the graph into subsets
if [ ! -e $prefix.subset.*.pmap ]; then
    partition-graph.py --threads 4 -s 1e5 $prefix
fi

# merge partitions
if [ ! -e $prefix.pmap.merged ]; then
    merge-partitions.py $prefix
fi

# annotate read with partition
if [ ! -e $base_paired_f.part ]; then
    annotate-partitions.py $prefix $base_paired_f $base_single_f
fi


# extract the partitions
extract-partitions.py $prefix.khmer_partition $base_paired_f.part $base_single.part
