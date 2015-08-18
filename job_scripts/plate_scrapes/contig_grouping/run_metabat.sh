
fasta=''
bams=''



MBIN="/netscr/hjcamero/meta_meta_assembly/meta_binning/metabat/bin"

function HELP {
	echo "
USAGE:

	run_metabat.sh -i assembly.fasta -b 'mybamdir/*.bam'

You may have to set the MBIN variable in the code to the bin/ of the metabat installation.

Run this with a complete node because metabat will use every core available (8).
Num threads can be specifed as a metabat option but it is not currently.
"
	
	exit 1
}

while getopts "i:b:" arg; do

    case $arg in
        i) # set scripts dir
            fasta=$OPTARG
            ;;
        
        b) # set bam_dir
            bams=$OPTARG
            ;;

        ?) # unknowns
            echo -e "\nOption \"$OPTARG\" not allowed."
            HELP
            ;;
    esac
done


for var in fasta bams; do
    if [ "${!var}" == '' ]; then
        echo -e "\n$var was not specified."
        HELP	
    fi
done


# make a coverage file
$MBIN/jgi_summarize_bam_contig_depths --outputDepth depth.txt $bams

# first binning with sensitive parameters (stringent)
$MBIN/metabat -i $fasta -a depth.txt -o bin_sens --sensitive -l -v --saveTNF saved.tnf --saveDistance saved.gprob

# second binning with specific parameters (build the bins as much as possible)
$MBIN/metabat -i $fasta -a depth.txt -o bin_spec --specific -l -v --saveTNF saved.tnf --saveDistance saved.gprob

