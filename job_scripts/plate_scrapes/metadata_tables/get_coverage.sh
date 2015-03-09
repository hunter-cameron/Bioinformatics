#!/usr/bin/env bash

set -e      # exit at any non-0 return

echo FASTA: $1
echo SAM: $2
echo outdir: $3


if [ ! -d $3 ]; then
    mkdir $3
fi

# convert sam to bam
if [ ! -e $3/alignment.bam ]; then
    echo Running... samtools view -S -b -o $3/alignment.bam $2
    samtools view -S -b -o $3/alignment.bam $2
fi

# sort bam
if [ ! -e $3/sorted_alignment.bam ]; then
    echo Running... samtools sort $3/alignment.bam $3/sorted_alignment
    samtools sort $3/alignment.bam $3/sorted_alignment
fi

# copy fasta so the faidx won't put it in the directory with the fasta
if [ ! -e $3/copy_scaffold.fasta ]; then
    echo Running... cp $1 $3/copy_scaffold.fasta
    cp $1 $3/copy_scaffold.fasta
fi

# begin making the genome file
if [ ! -e $3/copy_scaffold.fasta.fai ]; then
    echo Running... samtools faidx $3/copy_scaffold.fasta
    samtools faidx $3/copy_scaffold.fasta
fi

# cut out the columns I want
if [ ! -e $3/scaffold.genome ]; then
    echo Running... cut -f 1-2 $3/copy_scaffold.fasta.fai > $3/scaffold.genome
    cut -f 1-2 $3/copy_scaffold.fasta.fai > $3/scaffold.genome
fi

# get the number of reads that map back to each contig
if [ ! -e $3/scaffold_read_counts.txt ]; then
    echo Running... samtools index $3/sorted_alignment.bam
    samtools index $3/sorted_alignment.bam 
    echo "Running... samtools idxstats $3/sorted_alignment.bam | cut -f 1,3"
    samtools idxstats $3/sorted_alignment.bam | cut -f 1,3 > $3/scaffold_read_counts.txt
fi

# get the coverage
if [ ! -e $3/scaffold_depth.txt ]; then
    echo Running... bedtools genomecov -d -ibam $3/sorted_alignment.bam -g $3/genome.genome > $3/scaffold_depth.txt
    bedtools genomecov -d -ibam $3/sorted_alignment.bam -g $3/genome.genome > $3/scaffold_depth.txt
fi

### stripped perl script to parse contig coverage for mean and median
echo Running perl script...

perl -e '

use strict;
use warnings;
open(my $IN, "<", $ARGV[0]);

my %contig_stats;
my @bases;
my $previous = "";
my $contig = "";
while( readline $IN ) {
    chomp;
    ($contig, undef, my $coverage) = split "\t";
  
    if ($contig ne $previous) {
        if (scalar @bases) {
            $contig_stats{$previous}->{"median"} = median(\@bases);
            $contig_stats{$previous}->{"mean"} = $contig_stats{$previous}->{"total"} / $contig_stats{$previous}->{"count"};
        }
        else{
            if ($previous) {
                $contig_stats{$previous}->{"median"} = 0;
                $contig_stats{$previous}->{"mean"} = 0;
            }
        }
        @bases = [];
        $previous = $contig;

    }
    $contig_stats{$contig}->{"count"}++;
    $contig_stats{$contig}->{"total"} += $coverage;
    push(@bases, $coverage);
}

# calculate for the last contig
$contig_stats{$contig}->{"median"} = median(\@bases);
$contig_stats{$contig}->{"mean"} = $contig_stats{$contig}->{"total"} / $contig_stats{$contig}->{"count"};


close $IN;

# add the info for the read mapping
open($IN, "<", $ARGV[1]);

while ( readline $IN ) {
    chomp;
    my ( $contig, $count ) = split "\t";
    next if $contig eq "*";

    $contig_stats{$contig}->{"reads"} = $count;
}



print "contig\tmean_cov\tmed_cov\tnum_reads\n";
foreach ( sort keys %contig_stats ) {
    print "$_\t", $contig_stats{$_}->{"mean"}, "\t", $contig_stats{$_}->{"median"}, "\t", $contig_stats{$_}->{"reads"}, "\n";
}

sub median {
    my ($data) = @_;
    my @data = sort {$a <=> $b} @{$data};
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

' $3/scaffold_depth.txt $3/scaffold_read_counts.txt > $3/final_coverage_stats.txt


