#checks the format of the stats files

my ($folder) = @ARGV;

open my $IN, "<", "$folder/genome_stats.txt";

my @headers = qw(name size med_cov avg_cov code_density);
while ( my $line = readline $IN ) {
    chomp $line;

    my @fields = split "\t", $line;

    for ( my $i = 0; $i < @fields; $i++ ) {
        if ( ! defined $fields[$i] ) {
            print "$fields[0]:: $headers[$i]\n";
        }
    }
}

close $IN;



#check contigs
open $IN, "<", "$folder/contig_stats.txt";

@headers = qw(name size med_cov avg_cov code_density gc_content ref_Genome);
while ( my $line = readline $IN ) {
    chomp $line;

    my @fields = split "\t", $line;

    for ( my $i = 0; $i < @fields; $i++ ) {
        if ( ! defined $fields[$i] ) {
            print "$fields[6]-$fields[0]:: $headers[$i]\n";
        }
    }
}

close $IN;






#check genes
open $IN, "<", "$folder/gene_stats.txt";

@headers = qw(name size start_pos end_pos med_cov avg_cov product ref_Genome ref_Contig);
while ( my $line = readline $IN ) {
    chomp $line;

    my @fields = split "\t", $line;

    for ( my $i = 0; $i < @fields; $i++ ) {
        if ( ! defined $fields[$i] ) {
            print "$fields[7]-$fields[0]:: $headers[$i]\n";
        }
    }
}

close $IN;



print "DONE\n";
