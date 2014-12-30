#!/usr/bin/perl

use strict;
use warnings;

use Directory_obj;


my $usage = "Usage: genbank_write_by_header.pl <folder of genbanks>\n";

@ARGV == 1 or die "$usage \n";

my $directory_gbk = Directory_obj->new({'directory' => $ARGV[0]});

open(my $OUT, ">", "hypothetical_proteins.fasta");
while ( my $gbk = $directory_gbk->next_file() ) {
    my $filename = $directory_gbk->get_filename();
    print "gbk = $filename\n";
      
    use Bio::SeqIO;


    my $seqio_object = Bio::SeqIO->new( -file => $gbk ); 
    while ( my $seq_object = $seqio_object->next_seq() ) {

        for my $feat_object ($seq_object->get_SeqFeatures) {

            if ($feat_object->primary_tag() eq "CDS" or $feat_object->primary_tag() eq "rRNA" ) {
        
                #get gene tag to use as header
                my ($gene_tag) = $feat_object->get_tag_values('locus_tag') if ($feat_object->has_tag('locus_tag' ));
                 

                #$data->{$contig}{genes}{$gene_tag}{start_pos} = $feat_object->location()->start();
                #$data->{$contig}{genes}{$gene_tag}{end_pos} = $feat_object->location()->end();

                #get the product
                my ($product) = $feat_object->get_tag_values("product") if ($feat_object->has_tag("product"));
                if ($product =~ m/hypothetical protein/) {
                    my ($sequence) = $feat_object->get_tag_values("translation") if ($feat_object->has_tag("translation"));
                    print($OUT  ">$gene_tag\n$sequence\n");
                }   
            }
        }
    }
}
close $OUT;
