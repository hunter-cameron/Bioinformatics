#!/usr/bin/perl

###Eliminates contigs as plasmids by blasting to a database. High coverage, high identity matches indicate that the contig is chromosomal. Whatever does not match well is more likely to be plasmid. 


#reccommended database /nas02/data/blast/nt
#perhaps refseq_genomic if it wasn't so big

use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);
use Cwd qw(getcwd);

use Directory_obj;
use parallel_blast 'parallel_blast';

use Data::Dumper;

my $USAGE = "Usage: chromosome_reversal_forked.pl --db <database> --query <query.fasta> --out my_output.txt\n";

my $prefix = getcwd();

#Set up options
my $database;
my $query;
my $num_forks = 10;

my $MIN_LENGTH = 500;
my $MIN_IDENT = 90;
my $MIN_COV = 95;
my $out;
GetOptions(     'db=s' => \$database,
                'query=s' => \$query,
                'out=s' => \$out,
          );

if ( any {!defined} $database, $query ) {
    die "$USAGE";
}


my $directory_obj = Directory_obj->new({'directory' => $query,});

while ( my $file = $directory_obj->next_file() ) {
    my $filename = $directory_obj->get_filename();
    main($database, $file, $filename);
    last; 
}

#print_results($alignments);


    
sub main {
    my ($db, $query, $filename) = @_;

    if ( ! defined $out ) {

        $out = "$prefix/$filename\_reversal.txt";
        print "OUT not defined\n";
    }


    #make a hash of all sequences
    my %sequences;
    open my $IN, "<", $query;

    while( my $line =  readline $IN ) {
        chomp $line;
        if ( $line =~ m/\A>(.+)\z/ ) {  #get the header
            $sequences{$1} = 1;
        }
    }


    my $results = blast($db, $query, $out);



    my $plasmids = 0;
    my $total = 1;
    foreach my $key ( keys %sequences ) {
       
        #print "Key = $key\n";
        if (! defined $results->{$key}) {
            $results->{$key}{prediction} = "Plasmid";
            $results->{$key}{pident} = 0;
            $results->{$key}{qcovs} = 0;
            $plasmids++;
            $total++;
            next;

        }

        if ( $results->{$key}{pident} >= $MIN_IDENT and $results->{$key}{qcovs} >= $MIN_COV ) {
            $results->{$key}{prediction} = "Chromosome";
        }

        else {
            $results->{$key}{prediction} = "Plasmid";
            $plasmids++;
        }

        $total++;
        #print results

    }

    my $validity = ($total - $plasmids) / $total;

    open my $OUT, ">", "$out";
    print $OUT "Contig\tCoverage\tIdentity\tPrediction\t\%Chromosomal = $validity \n";

    foreach my $key ( keys %{$results} ) {


        print $OUT join "\t", ( $key,
                                $results->{$key}{qcovs},
                                $results->{$key}{pident},
                                $results->{$key}{prediction},
                                ), "\n";
    }



        #print Dumper $results;
        close $OUT;


}


sub blast {
    
    my ($database, $query, $out) = @_;


    #-culling_limit 1
    #prevents repeat aligns (removes a hit that is completely enveloped by a higher scoring hit
    my $command = join " ", (   "blastn",
                                "-task megablast",
                                "-db $database",
                                "-max_target_seqs 1",
                                #"-dust",   #remove low complexity from input
                                #"-num_threads 3",
                                "-query $query",
                                "-out $out",
                                "-outfmt '6 std qcovs'",
                                );

    #print "$command \n";

    parallel_blast::parallel_blast({num_forks => 10, command => $command});
     


    #if ( system($command) != 0 ) {
    #if ( system( "blastn -task megablast -max_target_seqs 10 -dust -db $database -query $query -out 'temp_results.out' -outfmt '6 std qcovs' -num_threads 3" ) != 0 )
         #{
     #   print "Blast could not be done successfully:\n blastn -db $database -query $query -out 'temp_results.out' -outfmt '6 std qcovs'\n";
     #   return 0;
    #}

    #process each line of the blast results
    open my $IN, "<", "$out";
    my %results;
    while ( my $line = readline $IN ) {
        chomp $line;
                                                                                         
        #print "$.\n";
        my @ln_elements = split "\t", $line;

        my $qseqid = $ln_elements[0];
        my $pident = $ln_elements[2];
        my $qcovs = $ln_elements[12];

        #my %ln_elements = ( 'qseqid'    =>  $ln_elements[0],
                            #'sseqid'    =>  $ln_elements[1],
         #                   'pident'    =>  $ln_elements[2],
                            #'length'    =>  $ln_elements[3],
                            #'mismatch'  =>  $ln_elements[4],
                            #'gapopen'   =>  $ln_elements[5],
                            #'qstart'    =>  $ln_elements[6],
                            #'qend'      =>  $ln_elements[7],
                            #'sstart'    =>  $ln_elements[8],
                            #'send'      =>  $ln_elements[9],
                            #'evalue'    =>  $ln_elements[10],
                            #'bitscore'  =>  $ln_elements[11],
          #                  'qcovs'    =>  $ln_elements[12],
        
                            #);

        if ( !defined $results{$qseqid}{pident} ) {
            $results{$qseqid}{pident} = $pident;
            $results{$qseqid}{qcovs} = $qcovs;
            next;
        }

        #check if there is a better match
        if ( $results{$qseqid}{pident} <= $pident and $results{$qseqid}{qcovs} <= $qcovs ) {
            $results{$qseqid}{pident} = $pident;
            $results{$qseqid}{qcovs} = $qcovs;
        }
    }
    close $IN;
    #unlink "temp_results.out";

    return \%results;
}


           
