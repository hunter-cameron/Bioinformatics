package PlasmidFinder::Gene;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;
use Scalar::Util qw(weaken);
{
    #Class Attributes
    my %gene_tag_of;
    my %name_of;
    my %size_of;
    my %start_of;
    my %end_of;
    my %product_of;
    my %mol_type_of;
    my %length_of;
    my %contig_of;
    my %coverage_depth_of;

    sub new {
        my ($class, $args_href) = @_;
        
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{gene_tag}, $args_href->{contig} ) {
            die "Cannot create new Gene. Essential parameter undefined.\n";
        }

        #make new object
        my $new_obj = bless \do{my $anon_scalar}, $class;


        #set parameters for object
        $gene_tag_of{ident $new_obj} = $args_href->{gene_tag};
        $contig_of{ident $new_obj} = $args_href->{contig};

        #weaken the backreference to contig so contig can be garbage collected
        weaken $contig_of{ident $new_obj};
        return $new_obj;
    }


    #                  #
    ##                ##
    ###Getter Methods###
    ##                ##
    #                  #




    sub get_name { 
        my ($self) = @_;
        if ( defined $name_of{ident $self}) {
            return $name_of{ident $self};
        }
        else {
            return 'null';
        }
    }
    
    sub get_gene_tag { 
        my ($self) = @_;
        if ( defined $gene_tag_of{ident $self}) {
            return $gene_tag_of{ident $self};
        }
        else {
            return 'null';
        }
    }

    sub get_contig {
        my ($self) = @_;

        return $contig_of{ident $self};
    }


    sub get_start {
        my ($self) = @_;
        
        if ( defined $start_of{ident $self} ) {
            return $start_of{ident $self};
        }
        else {
            return 'null';
        }
    }

    sub get_end {
        my ($self) = @_;
        
        if ( defined $end_of{ident $self} ) {
            return $end_of{ident $self};
        }
        else {
            return 'null';
        }
    }

    sub get_length {
        my ($self) = @_;
        
        if ( defined $length_of{ident $self} ) {
            return $length_of{ident $self};
        }
        else {

            $length_of{ident $self} = $self->get_end() - $self->get_start();
            return $length_of{ident $self};
        }
    }
        

    sub get_product {
        my ($self) = @_;
        
        if ( defined $product_of{ident $self} ) {
            return $product_of{ident $self};
        }
        else {
            return 'null';
        }
    }  

    sub get_mol_type {
        my ($self) = @_;
        
        if ( defined $mol_type_of{ident $self} ) {
            return $mol_type_of{ident $self};
        }
        else {
            return 'null';
        }
    }

    sub get_coverage_depth {
        my ($self) = @_;
        
        if ( defined $coverage_depth_of{ident $self} ) {
            return $coverage_depth_of{ident $self};
        }
        else {
           
            ($coverage_depth_of{ident $self}, my $median) = $self->get_contig()->_get_gene_coverage_depth($self);

            return ($coverage_depth_of{ident $self}, $median);
        }
    }



    #                  #
    ##                ##
    ###Setter methods###
    ##                ##
    #                  #

    sub set_name {
        my ($self, $name) = @_;

        $name_of{ident $self} = $name;
        return;
    }


    sub set_start {
        my ($self, $start) = @_;

        $start_of{ident $self} = $start;
        return;
    }

    sub set_end {
        my ($self, $end) = @_;

        $end_of{ident $self} = $end;
    }

    sub set_product {
        my ($self, $product) = @_;
        
        $product_of{ident $self} = $product;
    }

    sub set_mol_type {
        my ($self, $mol_type) = @_;

        $mol_type_of{ident $self} = $mol_type;
    }
    
    
    
    
    
    #                    #
    ##                  ##
    ###Internal methods###
    ##                  ##
    #                    #


   sub DESTROY {
       my ($self) = @_;

        delete $gene_tag_of{ident $self};
        delete $name_of{ident $self};
        delete $size_of{ident $self};
        delete $start_of{ident $self};
        delete $end_of{ident $self};
        delete $product_of{ident $self};
        delete $mol_type_of{ident $self};
        delete $length_of{ident $self};
        delete $contig_of{ident $self};
        delete $coverage_depth_of{ident $self};
    }





}
1;
