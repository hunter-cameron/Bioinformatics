#!/usr/bin/perl

use strict;
use warnings;

use BioUtils::FastaIO;
use BioUtils::FastaSeq;

@ARGV == 2 or die "Usage: cbar_test.pl <original.fasta> <new.fasta>\n";

my %lt1k =          ( 'name' => '<1k' );
my %btw1k_5k =      ( 'name' => '1k::5k' );
my %btw5k_10k =     ( 'name' => '5k::10k' );
my %btw10k_20k =    ( 'name' => '10k::20l' );
my %btw20k_50k =    ( 'name' => '20k::50k' );
my %btw50k_100k =   ( 'name' => '50k::100k' );
my %gt100k =        ( 'name' => '>100k' );
my %total =         ( 'name' => 'total' );


for ( my $i = 0; $i < 100; $i++) {
    my $in = BioUtils::FastaIO->new({stream_type => '<', file => $ARGV[0]});
    my $out = BioUtils::FastaIO->new({stream_type => '>', file => $ARGV[1]});

    while ( my $seq_obj = $in->get_next_seq() ) {
        #limit it to ~100 output contigs (total plasmids in file = 4378)
        if ( int(rand 45) != 1) {
            next;
        }
        else {
            my $header = $seq_obj->get_header();
            my $seq = $seq_obj->get_seq();
    
            my $length = length $seq;

            my $start = int(rand $length);
            my $num_bases = int(rand($length - $start));
        
            $seq = substr $seq, $start, $num_bases;
            
            $header = $header . " | _start_$start";

            $seq_obj->set_header($header);
            $seq_obj->set_seq($seq);
    
            $out->write_seq($seq_obj);
        }

    }

    `/nas02/home/h/j/hjcamero/scripts/plasgraph/cBar.1.2/cBar.pl $ARGV[1] cbout.txt`;
    open my $IN, "<", "cbout.txt";
    while ( my $line = readline $IN ) {
        chomp($line);
        my ($id, $length, $type) = split "\t", $line;
        next if ( $length eq "Length" );
        my $key;
        print "$type\n";
        if ( $type eq "Plasmid" ) {
            $key = "correct";
        }
        else {
            $key = "wrong"
        }

        print "$length\n";

        $lt1k{$key}++ if ( $length < 1000 );
        $btw1k_5k{$key}++ if ( $length >= 1000 and $length < 5000 );
        $btw5k_10k{$key}++ if ( $length >= 5000 and $length < 10000 );
        $btw10k_20k{$key}++ if ( $length >= 10000 and $length < 20000 );
        $btw20k_50k{$key}++ if ( $length >= 20000 and $length < 50000 );
        $btw50k_100k{$key}++ if ( $length >= 50000 and $length < 100000 );
        $gt100k{$key}++ if ( $length >= 100000 );

        $total{$key}++


        }
    close $IN;
}

open my $OUT, ">", "results.txt";
my @hrefs = (\%lt1k, \%btw1k_5k, \%btw5k_10k, \%btw10k_20k, \%btw20k_50k, \%btw50k_100k, \%gt100k, \%total);

print $OUT "Range\tCorrect\tWrong\tPercent\n";

foreach my $href (@hrefs) {
    my $total = $href->{'correct'} + $href->{'wrong'};
    my $perc_correct = ($href->{'correct'} / $total) * 100;
    print $OUT $href->{'name'}, "\t", $href->{'correct'}, "\t", $href->{'wrong'}, "\t", $perc_correct, "\n";

    delete $href->{'correct'};
    delete $href->{'wrong'};
    
}

close $OUT;

