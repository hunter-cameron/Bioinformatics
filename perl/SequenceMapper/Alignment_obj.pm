package SequenceMapper::Alignment_obj;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;

{
    #Class Attributes
    my %name_of;
    my %start_pos;
    my %end_pos;
    my %qstart_of;
    my %qend_of;
    my %perc_iden;
    my %subj_of;

    sub new {
        my ($class, $args_href) = @_;
        
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{name}, $args_href->{s_start}, $args_href->{s_end} ) {
            die "Cannot create new Alignment. Essential parameter undefined.\n";
        }


        #bless takes two args, a reference to the variable and a string containing th ename of the class
        #/do{my $anon_scalar} = reference to a scalar that only has scope within the do statement. So the object doesn't "exist" after the end of the do but it still kept alive because of the blessed reference to it
        #could be done with the anon_scalar() method from class::std::utils
        my $new_obj = bless \do{my $anon_scalar}, $class;

        #set parameters for object; length() returns a unique id for the object
        $name_of{ident $new_obj} = $args_href->{name};

        #check for reverse alignments 
        if ( $args_href->{s_start} > $args_href->{s_end} ) {
            #switch the start and end parameters for reverse alignments
            $start_pos{ident $new_obj} = $args_href->{s_end};
            $end_pos{ident $new_obj} = $args_href->{s_start};
        }
        else {
            $start_pos{ident $new_obj} = $args_href->{s_start};
            $end_pos{ident $new_obj} = $args_href->{s_end};
        }
        
        #set optional parameters
        $perc_iden{ident $new_obj} = $args_href->{perc_iden} if ( defined $args_href->{perc_iden} );
        $subj_of{ident $new_obj} = $args_href->{subj} if ( defined $args_href->{subj} ); 
        $qstart_of{ident $new_obj} = $args_href->{q_start} if ( defined $args_href->{q_start} );
        $qend_of{ident $new_obj} = $args_href->{q_end} if ( defined $args_href->{q_end} );

        return $new_obj;
    }


    sub get_name { 
        my ($self) = @_;
        if ( ! defined $name_of{ident $self}) {
            return "Header of object is not defined. Make sure object was created. \n";
        }
        else {
            return $name_of{ident $self};
        }
    }

    sub get_subject {
        my ($self) = @_;

        return $subj_of{ident $self};
    }

    sub get_length {
        my ($self) = @_;
        
        my $start = $start_pos{ident $self};
        my $end = $end_pos{ident $self};

        my $length = $end - $start;

        return $length;

    }

    sub get_start {
        my ($self) = @_;

        return $start_pos{ident $self};
    }

    sub get_end {
        my ($self) = @_;

        return $end_pos{ident $self};
    }



    sub get_qstart {
        my ($self) = @_;

        return $qstart_of{ident $self};

    }

    sub get_qend {
        my ($self) = @_;

        return $qend_of{ident $self};

    }

}
1;
