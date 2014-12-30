#!/usr/bin/perl

use strict;
use warnings;

use BioUtils::FastaIO;
use BioUtils::FastaSeq;

@ARGV == 2 or die "Usage: add_length_to_header.pl <original.fasta> <new.fasta>\n";

my $in = BioUtils::FastaIO->new({stream_type => '<', file => $ARGV[0]});
my $out = BioUtils::FastaIO->new({stream_type => '>', file => $ARGV[1]});

while ( my $seq_obj = $in->get_next_seq() ) {
    my $header = $seq_obj->get_header();
    my $seq = $seq_obj->get_seq();
    
    my $length = length $seq;

    $header = $header . " | len_$length";

    $seq_obj->set_header($header);
    
    $out->write_seq($seq_obj);
}

