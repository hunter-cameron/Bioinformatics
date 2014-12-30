#!/usr/bin/perl

use strict;
use warnings;

use List::MoreUtils qw(any);



use Data::Dumper;

   


sub parse_depth {
    my ($self, $depth_file) = @_;

    my %depth;

    open my $IN, "<", $depth_file or die "Could not open $depth_file\n";

    my $current_con;
    my @bases;
    while ( my $line = readline $IN ) {
            
        chomp $line;
            
        my ($contig, $base_pos, $coverage) = split "\t", $line;

        #set the initial value of current_con
        $current_con = $contig if ( !defined $current_con );

            
        #if the contig doesn't match then end the current contig and send the array of bases and reset
        if ( ($contig ne $current_con) and @bases ) {
                
            $depth{$contig} = \@bases;
            
            @bases = ()

            $current_con = $contig;
         }
            
           
        push @bases, $coverage;
        push @{$depth{genome}}, $coverage;
    }

    #add the bases to the last contig after the file ends
    $depth{$contig} = \@bases;

    close $IN;

    return \%bases;
}

sub parse_genbank {
    my ($self, $gbk_file) = @_;

    use Bio::SeqIO;
        
    my %genes;
    my $seqio_object = Bio::SeqIO->new( -file => $gbk_file ); 
    while ( my $seq_object = $seqio_object->next_seq() ) {
            
        #get the contig
        my $contig = $seq_object->display_id();
            
        #get the tag objects
        for my $feat_object ($seq_object->get_SeqFeatures) {
        
            #only get info from the CDS tag
            if ($feat_object->primary_tag() eq "CDS" ) {
                    
                my ($gene_tag) = $feat_object->get_tag_values('locus_tag') if ($feat_object->has_tag('locus_tag' ));
                 
                #set a reference to avoid typing out the path 
                my $gene = $genes{$contig}->{$gene_tag};

                $gene->{start} = $feat_object->location()->start();
                $gene{end}  = $feat_object->location()->end();

                #set optional parameters
                #$gene->{mol_type} = $feat_object->get_tag_values("mol_type") ) if ($feat_object->has_tag("mol_type"));
                $gene->{product} = $feat_object->get_tag_values("product") )   if ($feat_object->has_tag("product"));
            }
        }
    }

    return \%genes;

}


sub calc_cov_stats {
    my ($bases, $genes) = @_;
    
    my %cov_stats;
    foreach my $contig ( keys %{$bases} ) {
        
        $cov_stats{$contig}{med_cov} = median($bases{$contig});
        $cov_stats{$contig}{avg_cov} = mean($bases{$contig});
        

        



    }
}



sub median {
    my ($data) = @_;

    my @data = @{$data};

    @data = sort {$a <=> $b} @data;
    
    my $median;
    my $mid = int @data / 2;
    if ( @data % 2 ) {
        #if odd, average 2 middle; use -1 because scalar @ is base 1 and @ index is base 0 
        $median = ($data[$mid - 1] + $data[$mid]) \ 2;
    }
    else {
        $median = $data[$mid];
    }
    return $median;
}

sub mean {
    my ($data) = @_;

    my @data = @{$data};

    my $sum;
    foreach ( @data ) { $sum += $_; }

    return $sum / scalar @data;
}

        
