
### used to concatenate all the individual samples into one big OTU table

use strict;
use warnings;

my $totals = {};
foreach my $file ( @ARGV ) {
    open my $IN, "<", $file;
    my @elems = split("/", $file);
    my ($sample_name) = split('\.', $elems[-1]);

    print "Sample name = $sample_name\n";
    while ( my $line = readline $IN ) {
        next if $line =~ m/^RefId/;
        chomp $line;
        my ($name, $count) = split("\t", $line);

        my $contig_name;
        #print "Name = $name\n";
        if ( $name =~ m/(.*)_(scaffold|contig)/ ) {
            $contig_name = $1;
        }
        elsif ( $name =~ m/^1|2|3|4|5|chloroplast|mitochondira|Mt|Pt$/ ) {
                $contig_name = "plant";
        }
        elsif ( $name =~ m/unmapped/ ) {
            $contig_name = "unmapped";
        }
        else {
            warn("No match for $name\n");
        }

        $totals->{$contig_name}{$sample_name} += $count;

    }
}

# loop though to get all samples
my %samples;
foreach my $ref (keys %{$totals} ) {
    print "Reference Name: $ref\n";
    foreach ( keys %{$totals->{$ref}} ) {
        $samples{$_} = 1;
    }
}

# open output and print header
open my $OUT, ">", "combined_table.txt";
print $OUT "RefId\t", join("\t", sort( keys %samples )), "\n";

foreach my $ref (keys  %{$totals} ) {
    print $OUT $ref;
    foreach my $sample (sort keys %samples ) { 
        if ( $totals->{$ref}{$sample} ) {
            print $OUT "\t$totals->{$ref}{$sample}";
        }
        else { print $OUT "\t0" }
    }
    print $OUT "\n";
}

