#!/usr/bin/env bash

set -e      # exit at any non-0 return

echo Ref: $1
echo FQ1: $2
echo FQ2: $3
echo outdir: $4


if [ ! -d $4 ]; then
    mkdir $4
fi

# make the DB
if [ ! -e $4/_database.bwt ]; then
    echo Running... bwa index -p $4/_database $fasta
    bwa index -p $4/_database $1
fi


# align the reads
if [ ! -e $4/alignment.sam ]; then
    echo Running... bwa mem -M  $4/_database $2 $3 > $4/alignment.sam
    bwa mem -M  $4/_database $2 $3 > $4/alignment.sam
fi

# convert sam to bam
if [ ! -e $4/alignment.bam ]; then
    echo Running... samtools view -S -b -o $4/alignment.bam $4/alignment.sam
    samtools view -S -b -o $4/alignment.bam $4/alignment.sam
fi

# sort bam
if [ ! -e $4/sorted_alignment.bam ]; then
    echo Running... samtools sort $4/alignment.bam $4/sorted_alignment
    samtools sort $4/alignment.bam $4/sorted_alignment
fi

# copy fasta so the faidx won't put it in the directory with the fasta
if [ ! -e $4/copy_scaffold.fasta ]; then
    echo Running... cp $1 $4/copy_scaffold.fasta
    cp $1 $4/copy_scaffold.fasta
fi

# begin making the genome file
if [ ! -e $4/copy_scaffold.fasta.fai ]; then
    echo Running... samtools faidx $4/copy_scaffold.fasta
    samtools faidx $4/copy_scaffold.fasta
fi

# cut out the columns I want
if [ ! -e $4/scaffold.genome ]; then
    echo Running... cut -f 1-2 $4/copy_scaffold.fasta.fai > $4/scaffold.genome
    cut -f 1-2 $4/copy_scaffold.fasta.fai > $4/scaffold.genome
fi

# get the coverage
if [ ! -e $4/scaffold_depth.txt ]; then
    echo Running... bedtools genomecov -d -ibam $4/sorted_alignment.bam -g $4/genome.genome > $4/scaffold_depth.txt
    bedtools genomecov -d -ibam $4/sorted_alignment.bam -g $4/genome.genome > $4/scaffold_depth.txt
fi

### stripped perl script to parse contig coverage for mean and median
### code still looks comprehensible... I need to work on my perl-ing! (:
echo Running perl script...
perl -e '

open(my $IN, "<", $ARGV[0]);

my %contig_stats;
my @bases;
my $previous = "";
my $contig == "";
while( readline $IN ) {
    chomp;
    ($contig, undef, my $coverage) = split "\t";
   
    if ($contig ne $previous) {
        if (@bases) {
            $contig_stats->{$contig}{"median"} = median(\@bases);
            $contig_stats->{$contig}{"mean"} = $contig_stats->{$contig}{"total"} / $contig_stats->{$contig}{"count"};
        }
        @bases = [];
        $previous = $contig;

    }
    $counts{$contig}->{"count"}++;
    $counts{$contig}->{"total"} += $coverage;
    push(@bases, $coverage);
}

# calculate for the last contig
$contig_stats->{$contig}{"median"} = median(\@bases);
$contig_stats->{$contig}{"mean"} = $contig_stats->{$contig}{"total"} / $contig_stats->{$contig}{"count"};


print "contig\tmean_cov\tmed_cov\n";
foreach ( sort keys %contig_stats ) {
    print "$_\t", $contig_stats->{$_}{"mean"}, "\t", $contig_stats->{$_}{"median"}, "\n";
}

sub median {
    my @data = @{$_};
    @data = sort {$a <=> $b} @data;
    my $median;
    my $mid = int @data / 2;
    if ( @data % 2 ) {
        #if data set is odd, use middle
        return $data[$mid];
    }
    else {
        #if even, average 2 middle; use -1 because scalar @ is base 1 and @ index is base 0 
        return ($data[$mid - 1] + $data[$mid]) / 2;
    }
}



' $4/scaffold_depth.txt > $4/depth_per_contig.txt


