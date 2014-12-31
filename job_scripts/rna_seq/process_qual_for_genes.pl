
### used to concatenate contigs of quality results

use strict;
use warnings;

foreach my $file ( @ARGV ) {
    open my $IN, "<", $file;
    my @elems = split("/", $file);
    my ($sample_name) = split('\.', $elems[-1]);

    print "Sample name = $sample_name\n";
    
    my $quality = {};
    open my $OUT, ">", "$sample_name\_map_quality.txt";
    #print $OUT "Genome\tPosition\tquality\n";
    while ( my $line = readline $IN ) {
        chomp $line;
        my ($query, $name_pos, $qual) = split("\t", $line);

        my ($name, $pos) = split("\_p", $name_pos);

        my $contig_name;
        my $contig_number;
        #print "Name = $name\n";
        if ( $name =~ m/(.*_)(scaffold|contig)(_?\d+\.?\d+?)/ ) {       # this regex is experimental now (adding trailing digits
            $contig_name = $1 . $2 . $3;
        }
        elsif ( $name =~ m/^(1|2|3|4|5|chloroplast|mitochondria|Mt|Pt)$/ ) {
                next;       # skip plant reads
                $contig_name = "plant";
                $contig_number = $1;
        }
        else {
            warn("No match for $name\n");
        }
        
        #$quality->{$contig_name}{$qual}++;
        print $OUT "$sample_name\t$query\t$contig_name\t$pos\t$qual\n";
         
    }

}



__END__
# print quality
    
    my %temp;
    foreach my $contig ( keys %{$quality} ) {
        foreach ( keys %{$quality->{$contig}} ) {
            $temp{$_} = 1;
        }
    }
    
    my @quality = sort { $a <=> $b } keys %temp;
    print $OUT "Contig\t", join("\t", @quality), "\n";

    foreach my $contig ( keys %{$quality} ) {
        print $OUT "$contig";
        foreach (@quality) {
            if ( defined $quality->{$contig}{$_} ) {
                print $OUT "\t$quality->{$contig}{$_}";
            }
            else {
                print $OUT "\t0";
            }
        }
        print $OUT "\n";
    }
}
