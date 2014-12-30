package PlasmidFinder::Genome;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;
use Statistics::Basic qw(mean stddev median);


use PlasmidFinder::Contig;
use PlasmidFinder::Gene;

use Data::Dumper;
{
    #Class Attributes
    my %name_of;
    my %size_of;
    my %contigs_of;
    #my %mean_contig_depth_of;
    my %coding_density_of;
    my %coverage_depth_of;


    sub new {
        my ($class, $args_href) = @_;
        
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{genbank}, $args_href->{name} ) {
            die "Cannot create new Genome. Essential parameter undefined.\n";
        }

        #make new object
        my $new_obj = bless \do{my $anon_scalar}, $class;

        $new_obj->_parse_genbank( $args_href->{genbank} );
        
        #if provided with a depth file, go ahead and parse it
        $new_obj->_parse_depth( $args_href->{depth} ) if ( defined $args_href->{depth} );


        #set parameters for object
        $name_of{ident $new_obj} = $args_href->{name};

        return $new_obj;
    }

    #                  #
    ##                ##
    ###Getter Methods###
    ##                ##
    #                  #

    sub get_name {
        my ($self) = @_;

        return $name_of{ident $self};
    }

    sub get_length {
        my ($self) = @_;

        my $length;
        foreach my $contig ( $self->get_contigs() ) {
            my @bases = $contig->get_bases();
            
            $length += scalar @bases;

        }
        return $length;
    }


        

    sub get_contig_names {

        my ($self) = @_;

        return keys %{$contigs_of{ident $self}};
    }

    sub get_contigs {

        my ($self) = @_;

        return values %{$contigs_of{ident $self}};
    }

    #calc once
    sub get_mean_contig_depth {

        my ($self) = @_;

        my @depths;

        foreach my $contig ( $self->get_contigs() ) {

            push @depths, $contig->get_coverage_depth();
        }

        my $mean = mean ( @depths );
        my $stddev = stddev ( @depths );

        return $mean, $stddev;
    }

    #not worth a hash but costly, calc once and store in var
    sub get_coverage_depth {

        my ($self) = @_;

        my $total_cov;
        my $total_bases;
        foreach my $contig ( $self->get_contigs() ) {

            foreach ( $contig->get_bases() ) {
                
                $total_cov += $_;
                $total_bases++;
            }
        }

        return $total_cov / $total_bases;
    }



    #once again, not worth hash but costly
    sub get_coding_density {

        my ($self) = @_;

        my $coding_bases;
        my $total_bases;

        foreach my $contig ( $self->get_contigs() ) {
            
            $total_bases += $contig->get_length();

            foreach my $gene ( $contig->get_genes() ) {

                $coding_bases += $gene->get_length;
            }
        }

        return $coding_bases / $total_bases;
    }


    sub get_median_coverage_depth {
        
        my ($self) = @_;

        my @bases;
        foreach my $contig ( $self->get_contigs() ) {

            foreach ( $contig->get_bases() ) {
            
                push @bases, $_;


            }
        }

        return median @bases;
    }       



    #                  #
    ##                ##
    ###Setter Methods###
    ##                ##
    #                  #


    #                    #
    ##                  ##
    ###Internal Methods###
    ##                  ##
    #                    #

    #calculates the average of the depths of the contigs (contig size independent) 
    sub _calc_mean_contig_depth {
        my ($self) = @_;

        my @depths;

        foreach my $contig ( $self->get_contig_refs() ) {

            push @depths, $contig->get_coverage_depth();
        }

        my $mean = mean ( @depths );
        my $stddev = stddev ( @depths );

        return $mean, $stddev;
    }

    


    sub _parse_depth {
        my ($self, $depth_file) = @_;

        open my $IN, "<", $depth_file or die "Could not open $depth_file\n";

        my $current_con;
        my @bases;
        while ( my $line = readline $IN ) {
            
            chomp $line;
            
            my ($contig, $base_pos, $coverage) = split "\t", $line;

            #set the initial value of current_con
            $current_con = $contig if ( !defined $current_con );

            #need to add the base values to the appropriate contig
            #would that be better done by calculating per contig here and passing an array to add bases or by each single base to the appropriate contig's add base?
            
            #if the contig doesn't match then end the current contig and send the array of bases and reset
            if ( ($contig ne $current_con) and @bases ) {
                
                #get the object reference 
                my $ref = $self->_get_contig_ref( $current_con );
                
                $ref->_add_bases(\@bases);
                
                @bases = ();

                $current_con = $contig;
            }
            
           
            push @bases, $coverage;
        }

        #add the bases to the last contig after the file ends
        my $ref = $self->_get_contig_ref( $current_con );
        $ref->_add_bases(\@bases);

        close $IN;

        return;
    }

    sub _get_contig_ref { 

        my ($self, $key) = @_;
        
        return $contigs_of{ident $self}->{$key} or "Die could not find a contig by the name of $key\n";
    }

    sub _parse_genbank {
        my ($self, $gbk_file) = @_;

        use Bio::SeqIO;

        my $seqio_object = Bio::SeqIO->new( -file => $gbk_file ); 
        my $contig_num = 1;    
        while ( my $seq_object = $seqio_object->next_seq() ) {
            
            #get the contig
            my $id = $seq_object->display_id();
            
            #add the contig if it isn't already
            my $contig = $self->_add_contig($id);
            #print "$id\n";


            my $gene_num = 1;
            #get the tag objects
            for my $feat_object ($seq_object->get_SeqFeatures) {
        
                #only get info from the CDS tag
                if ($feat_object->primary_tag() eq "CDS" ) {
                    
                    my ($gene_tag) = $feat_object->get_tag_values('locus_tag') if ($feat_object->has_tag('locus_tag' ));
                    my $name = "C$contig_num-G$gene_num"; 
                    my $gene = $contig->_add_gene($gene_tag);
                    
                    $gene->set_name ( $name );
                    $gene->set_start( $feat_object->location()->start() );
                    $gene->set_end( $feat_object->location()->end() );
                    

                    #set optional parameters
                    $gene->set_mol_type( $feat_object->get_tag_values("mol_type") ) if ($feat_object->has_tag("mol_type"));
                    $gene->set_product( $feat_object->get_tag_values("product") )   if ($feat_object->has_tag("product"));
                $gene_num++;
                }
            }
            $contig_num++;
        } 

    }



    sub _add_contig {
        my ($self, $contig_name) = @_;

        #print "$contig_name\n";
        #if the contig is already made, return
        return $contigs_of{ident $self}{$contig_name} if ( defined $contigs_of{ident $self}{$contig_name} );
        
        #make the new contig object if needed
        my $contig = PlasmidFinder::Contig->new( {'name' => $contig_name} );

        #add the object reference to the contig_name key 
        $contigs_of{ident $self}{$contig_name} = $contig;

        return $contig;
    }

    sub DESTROY {
        my ($self) = @_;

        #forced recycling of contigs and genes associated with Genome
        #foreach ( $self->get_contigs() ) {
        #     foreach ($_->get_genes() ) {
        #        $_->DESTROY();
        #    }
        #    $_->DESTROY();
        #}

        delete $name_of{ident $self};
        delete $size_of{ident $self};
        delete $contigs_of{ident $self};
        #delete $mean_contig_depth_of{ident $self};
        delete $coding_density_of{ident $self};
        delete $coverage_depth_of{ident $self};

    }



}
1;
