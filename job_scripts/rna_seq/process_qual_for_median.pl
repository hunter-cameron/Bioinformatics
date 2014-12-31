
### used to concatenate contigs of quality results

use strict;
use warnings;

my $quality = {};

open my $OUT, ">", "mean_qualities.txt";

foreach my $file ( @ARGV ) {
    open my $IN, "<", $file;
    my @elems = split("/", $file);
    my ($sample_name) = split('\.', $elems[-1]);

    print "Sample name = $sample_name\n";
    
    while ( my $line = readline $IN ) {
        chomp $line;
        my ($read, $name_pos, $qual) = split("\t", $line);

        my ($name, $pos) = split("\_p", $name_pos);

        my $contig_name;
        my $contig_number;
        #print "Name = $name\n";
        if ( $name =~ m/(.*)_(scaffold|contig)_?(\d+)/ ) {
            $contig_name = $1;
            $contig_number = $3;
        }
        elsif ( $name =~ m/^(1|2|3|4|5|chloroplast|mitochondira|Mt|Pt)$/ ) {
                next;       # skip plant reads
                $contig_name = "plant";
                $contig_number = $1;
        }
        else {
            warn("No match for $name\n");
        }
        
        $quality->{$sample_name}{$contig_name}{total} += $qual;
        $quality->{$sample_name}{$contig_name}{count}++;
         
    }

}


# print quality
    
my %temp;
foreach my $sample ( keys %{$quality} ) {
    foreach ( keys %{$quality->{$sample}} ) {
        $temp{$_} = 1;
    }
}
    
my @isolates = sort { $a <=> $b } keys %temp;
print $OUT "Sample\t", join("\t", @isolates), "\n";

foreach my $sample( sort keys %{$quality} ) {
    print $OUT "$sample";
    foreach my $isolate (@isolates) {
        if ( defined $quality->{$sample}{$isolate}{count} ) {
            my $mean = $quality->{$sample}{$isolate}{total} / $quality->{$sample}{$isolate}{count};
            print $OUT "\t$mean";
        }
        else {
            print $OUT "\t0";
        }
    }
    print $OUT "\n";
}
