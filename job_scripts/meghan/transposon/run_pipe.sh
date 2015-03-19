
set -e

ref_fasta='/proj/dangl_lab/hunter/tasks/meghan_transposon_files/WCS417.fna'
ref_gbk='/proj/dangl_lab/hunter/tasks/meghan_transposon_files/WCS417r.gbk'

function USAGE {
    echo -e "
    $(basename $0) genomic_frags prefix

    Prefix is the prefix for all of the output files. If no prefix is given,
    the prefix used will be the name of the input fasta file without the ext.
    "
    exit 1
}


# check for input fasta
if [ -z $1 ]; then
    USAGE
else
    frags=$1
fi


# check for/set prefix
if [ -z $2 ]; then
    prefix=$(echo $(basename $in_fasta) | sed 's/\..*$//')
    echo $prefix
else
    prefix=$2
fi

#python /nas02/home/h/j/hjcamero/scripts/job_scripts/meghan/transposon/sequence_miner.py $in_fasta $prefix.genomic.fasta

#bbmap.sh ref=$ref_fasta in=$frags threads=4 ambiguous=best sam=1.4 out=$prefix.bbmap.sam local=f k=8 nodisk
bbmap.sh ref=$ref_fasta in=$frags threads=4 ambiguous=best sam=1.4 out=$prefix.bbmap.sam local=f nodisk

python /nas02/home/h/j/hjcamero/scripts/job_scripts/meghan/transposon/sam2mappos.py $prefix.bbmap.sam $ref_gbk > $prefix.mapping.csv

mv counts_per_gene.png $prefix.counts_per_gene.png
mv plot_by_position.png $prefix.counts_per_position.png
mv hist_of_mapping_quality.png $prefix.hist_of_mapping_quality.png
