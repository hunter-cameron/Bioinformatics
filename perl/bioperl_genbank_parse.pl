#!/usr/bin/perl

### Gets the sequence features listed in @get_tags and writes to a file


use strict;
use warnings;

use Bio::SeqIO;


my @get_tags = qw(locus_tag product translation);

open(my $OUT, ">", "output.txt");

@ARGV == 1 or die "Usage: bioperl_genbank_parse.pl <myfile.gbk>\n";

#create bioperl writer SeqIO
my $seqio_object = Bio::SeqIO->new(-file => $ARGV[0]);
my @parse_data = ();

#read through each sequence
while(my $seq_object = $seqio_object->next_seq) {
    
    foreach my $feat ($seq_object->get_SeqFeatures) {
        print("feat = $feat \n");
        #if ($feat->primary_tag eq "CDS") {
            
            #get sub-features
            foreach my $sub_feat ($feat->get_SeqFeatures) {
                print("\t subfeat = $sub_feat\n");
                foreach my $tag (@get_tags) {
                    #if ($sub_feat->has_tag($tag)) {
                        print($OUT "$tag\t=" . $sub_feat->get_tag_values($tag) . "\n");
                    #}
                }
                #separate each entry
                #print($OUT "=" x 20 . "\n\n");
            }
        #}
    }
    #separate each contig
    print($OUT "\n" . "=" x20 . "\n\n");

    last;
}
