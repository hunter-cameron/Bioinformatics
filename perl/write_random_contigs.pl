#!/usr/bin/perl

use strict;
use warnings;

use BioUtils::FastaIO;
use BioUtils::FastaSeq;

@ARGV == 2 or die "Usage: write_random_contigs.pl <original.fasta> <new.fasta>\n";

my $in = BioUtils::FastaIO->new({stream_type => '<', file => $ARGV[0]});
my $out = BioUtils::FastaIO->new({stream_type => '>', file => $ARGV[1]});

while ( my $seq_obj = $in->get_next_seq() ) {
    #limit it to ~100 output contigs (total plasmids in file = 4378)
    #if ( int(rand 45) != 1) {
        #next;
    #}
    #else {
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
    #}
}
