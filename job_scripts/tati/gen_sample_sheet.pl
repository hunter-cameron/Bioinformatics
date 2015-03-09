
use strict;
use warnings;

die "program_name.pl <metadata>" unless @ARGV == 1;

open(my $IN, "<", $ARGV[0]);

while( my $line = readline($IN) ) {
    chomp $line;

    my ($id, $barcode, $fwd_shift, $rev_shift, $indx_seq) = split("\t", $line);

    chop $indx_seq;
    my $sample_barcode = join("-", ($fwd_shift, $barcode, $rev_shift, $indx_seq));

    print $sample_barcode . "\t" . $id . "\n";
}

