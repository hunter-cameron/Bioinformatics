package PlasmidFinder::Contig;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;
use Statistics::Basic qw(mean stddev median);
use Carp qw(croak);


use PlasmidFinder::Gene;

use Data::Dumper;

{
    #Class Attributes
    my %name_of;
    my %size_of;
    my %bases_of;
    my %length_of;
    my %genes_of;
    my %coding_density_of;
    my %coverage_depth_of;

    sub new {
        my ($class, $args_href) = @_;
        
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{name} ) {
            die "Cannot create new Contig. Essential parameter undefined.\n";
        }

        #make new object
        my $new_obj = bless \do{my $anon_scalar}, $class;


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

    sub get_coverage_depth {
        
        my ($self) = @_;

        #return if already calculated
        if ( defined $coverage_depth_of{ident $self} ) {

            return $coverage_depth_of{ident $self};
        }

        #if not already calculated; calculate
        else {

            return $self->_calc_coverage_depth();
        }
    }

    sub get_median_coverage_depth {

        my ($self) = @_;

        return median $bases_of{ident $self};
    }


    sub get_coding_density {

        my ($self) = @_;

        #return if already calculated
        if ( defined $coding_density_of{ident $self} ) {

            return $coding_density_of{ident $self};
        }

        #if not already calculated; calculate
        else {

            return $self->_calc_coding_density();
        }
    }


    sub get_length {

        my ($self) = @_;
        
       #return if already calculated
        if ( defined $length_of{ident $self} ) {

            return $length_of{ident $self};
        }

        #if not already calculated; calculate
        else {

            return $self->_calc_length();
        }
    } 

    
    sub get_genes {
        
        my ($self) = @_;

        return values %{$genes_of{ident $self}};

    
    }

    sub get_bases {
        my ($self) = @_;

        return @{$bases_of{ident $self}};
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
    
    #This length is lower than the true length because it is derrived from the depth results which omit either cov = 0 or low quality reads; real length could be parsed from genbank?
    sub _calc_length {

        my ($self) = @_;
        
        return scalar @{$bases_of{ident $self}} or die "Could not calculate length\n";

    }

    sub _calc_coding_density {

        my ($self) = @_;

        my $coding_bases = 0;
        
        
        foreach my $ref ( $self->get_genes() ) {
            
            $coding_bases += $ref->get_length();
        }
        $coding_density_of{ident $self} = $coding_bases / $self->get_length();

        return $coding_density_of{ident $self};
    }



    sub _calc_coverage_depth {

        my ($self) = @_;

        #print Dumper %bases_of;
        if ( !defined $bases_of{ident $self} ) {
            print Dumper $bases_of{ident $self},
            die;
        }
        
        $coverage_depth_of{ident $self} = mean $bases_of{ident $self};
        return $coverage_depth_of{ident $self};
    }


    sub _add_gene {
        my ($self, $gene_tag) = @_;

        #if the gene is already made, return
        die "Already Made $gene_tag" if ( defined $genes_of{ident $self}{$gene_tag} );
        #return $genes_of{ident $self}{$gene_name} if ( defined $genes_of{ident $self}{$gene_name} );
        #make the new gene object if needed
        my $gene = PlasmidFinder::Gene->new( {'gene_tag' => $gene_tag, 'contig' => $self} );

        #add the object reference to the contig_name key 
        $genes_of{ident $self}{$gene_tag} = $gene;

        return $gene;
    }

    sub _get_gene_coverage_depth {
        my ($self, $gene) = @_;

        #get the start and end of gene, -1 to translate to base 0
        my $start = $gene->get_start() - 1;
        my $end = $gene->get_end() - 1;

        my $coverage = mean @{$bases_of{ident $self}}[$start..$end];
        my $median = median @{$bases_of{ident $self}}[$start..$end]; 
        return ($coverage, $median);
    }


    sub _add_bases {
        my ($self, $bases) = @_;
      
        #these lines make a dereferenced copy of bases and then send the reference to this copy to the hash; needs to be done or there will be multiple references to the same array of bases. 
        my @bases = @{$bases}; 
        $bases_of{ident $self} = \@bases;
        
        return;
    }

    sub DESTROY {
        my ($self) = @_;

        delete $name_of{ident $self};
        delete $size_of{ident $self};
        delete $bases_of{ident $self};
        delete $length_of{ident $self};
        delete $genes_of{ident $self};
        delete $coding_density_of{ident $self};
        delete $coverage_depth_of{ident $self};
    }


}
1;
