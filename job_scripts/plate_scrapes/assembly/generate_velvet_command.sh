

# generates velvet command that will run a whole assembly
# this is for my convenience, don't use this script unless you know exactly what it is doing

echo "data dir = $1"
echo "asm dir = $2"


if [ -z $1 ]; then
    echo "Positional argument 1, data dir, is required."
    exit 1
fi

if [ -z $2 ]; then
    echo "Positional argument 2, asm dir, is required."
fi

base_cmd="bsub -o stdout -e stderr -q bigmem -n 8 -R 'span[hosts=1]' -M 500 -J meta_Assembly \"velveth $2 101"

indx=1
cd $1
for file in *khmer_interleaved*.fastq.gz; do
    base_cmd+=" -shortPaired$indx -fastq.gz -interleaved $1/$file"
    ((indx+=1))
done

indx=1
for file in *khmer_singles*.fastq.gz; do
    base_cmd+=" -short$indx -fastq.gz $1/$file"
    ((indx+=1))
done


base_cmd+="; velvetg $2 -exp_cov auto -cov_cutoff auto\""

echo $base_cmd
