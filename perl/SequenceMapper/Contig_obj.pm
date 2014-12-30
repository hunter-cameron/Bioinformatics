package SequenceMapper::Contig_obj;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;
use Data::Dumper;

use SequenceMapper::Alignment_obj;

{
    #Class Attributes
    my %length_of;
    my %alignments_of;
    my %start_of;
    my %end_of;
    
    sub new {
        my ($class, $args_href) = @_;
        #print "$class\n";
        #print "$args_href->{'name'}\n";    
        #check if essential parameters defined, could use to give better error report
        #if ( any {! defined} $args_href->{start}, $args_href->{end} and ! defined $args_href->{from_alignment} and ! defined $args_href->{from_build}) {
        if (! defined $args_href->{from_build}) {
            print $args_href->{from_build}, "\n";
            die "Cannot create new Sequence object. Parameter undefined.\n";
        }


        #bless takes two args, a reference to the variable and a string containing th ename of the class
        #/do{my $anon_scalar} = reference to a scalar that only has scope within the do statement. So the object doesn't "exist" after the end of the do but it still kept alive because of the blessed reference to it
        #could be done with the anon_scalar() method from class::std::utils
        my $new_obj = bless \do{my $anon_scalar}, $class;
        
        if ( defined $args_href->{from_alignment} ) {
           $new_obj->_new_from_alignment($args_href->{from_alignment});
        }

        elsif ( defined $args_href->{from_build} ) {
            $new_obj->_new_from_build($args_href->{from_build});
        }

        else {

        #set parameters for object; length() returns a unique id for the object
        $start_of{ident $new_obj} = $args_href->{start};
        $end_of{ident $new_obj} = $args_href->{end};

        #store a reference to a anon has within an anon array for the alignments hash
        #array to order the alignments along the contig; hash to store alignment ID, start, and end
        $alignments_of{ident $new_obj} = [{}]; 
        }
        
        return $new_obj;
    }



    

    sub _new_from_build {
        my ($self, $contig_array) = @_;

        $start_of{ident $self} = $contig_array->[0]{'start'} or die Dumper $contig_array, "Missing Start\n";
        $end_of{ident $self} = $contig_array->[-1]{'end'} or die "missing end\n";
        $alignments_of{ident $self} = $contig_array;
        $self->set_length();
    }


    sub _new_from_alignment {
        my ($self, $alignment) = @_;
        
        #set initial contig params to the alignment (since the alignment makes up the whole contig)
        $start_of{ident $self} = $alignment->get_start();
        $end_of{ident $self} = $alignment->get_end();
        
        

        #add the new alignment hash to the array
        #would making dual calls to alignment be better or one call in beginning to set a scalar?
        push @{$alignments_of{ident $self}}, {'id' => $alignment,
                                              'start' => $alignment->get_start(),
                                              'end' => $alignment->get_end(),};
    }

    
    sub set_length {
        my ($self) = @_;

        my $length = $end_of{ident $self} - $start_of{ident $self};

        $length_of{ident $self} = $length;
    }

    sub get_start {
        my ($self) = @_;

        return $start_of{ident $self};
    }

    sub get_end {
        my ($self) = @_;

        return $end_of{ident $self};
    }

    sub get_length {
        my ($self) = @_;
        if ( ! defined $length_of{ident $self} ) {
            return "Length of object is not defined. Make sure object was created. \n";
        }
        else {
            return $length_of{ident $self};
        }
    }



    #something needs to be done with this to return something more sensible than an anon array of anon hashes
    #perhaps return an array of alignment references? - but then misleading info about the start/ends
    sub get_alignments {
        my ($self) = @_;
        
        #Do I want to return an array? or the reference?
        #return the ref because it is better on memory
        return $alignments_of{ident $self};
    }


    

    



    #delete all attributes of the object; won't destroy the scalar object through b/c there is still ref to it in the main program?
    sub DESTROY {
        my ($self) = @_;

        delete $length_of{ident $self};
        delete $alignments_of{ident $self};
        delete $start_of{ident $self};
        delete $end_of{ident $self};

    }







}








1;
