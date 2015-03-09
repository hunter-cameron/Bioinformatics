
use strict;
use warnings;

use BioUtils::FastqIO;

# output file
my $out = BioUtils::FastqIO->new({stream_type => ">", file => "multiplexed_seqs.fastq"});

# loop through all fastqs
foreach my $file (@ARGV) {
    my $in = BioUtils::FastqIO->new({stream_type => '<', file => $file});
    
    # loop through all seqs -- could be done with a regex instead of splits
    while( my $seq_obj = $in->get_next_seq() ) {
        my $header = $seq_obj->get_header();
        #print($header);
        my ($name, $read) = split(" ", $header);
        #print($read);
        (undef, undef, undef, my $well_code) = split(":", $read);

        # append well_code to the beinning of the header
        $seq_obj->set_header($well_code . "_" . $header);

        $out->write_seq($seq_obj)
    }
}

