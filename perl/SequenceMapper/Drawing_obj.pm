package SequenceMapper::Drawing_obj;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;
use GD::Simple;

{
    #Class Attributes
    my %name_of;
    my %canvas_of;
    my %end_pos;
    my %perc_iden;

    my $Y = 480;
    my $X = 1000;
    
    sub new {
        my ($class, $args_href) = @_;
        
        #check if essential parameters defined, could use to give better error report
        #if ( any {! defined} $args_href->{name}, ) {
        #    die "Cannot create new Drawing. Essential parameter undefined.\n";
        #}


        #bless takes two args, a reference to the variable and a string containing th ename of the class
        #/do{my $anon_scalar} = reference to a scalar that only has scope within the do statement. So the object doesn't "exist" after the end of the do but it still kept alive because of the blessed reference to it
        #could be done with the anon_scalar() method from class::std::utils
        my $new_obj = bless \do{my $anon_scalar}, $class;

        #set parameters for object; length() returns a unique id for the object
        $name_of{ident $new_obj} = $args_href->{name} if ( defined $args_href->{'name'} );

        $canvas_of{ident $new_obj} = GD::Simple->new($X, $Y);

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

    sub draw_contigs {
        my ($self, $sequence) = @_;

        my $canvas = $canvas_of{ident $self};
        
        #in to top left, write the sequence name
        $canvas->bgcolor('black');
        $canvas->moveTo(50, 50);
        $canvas->string($sequence->get_name());

        #graph a representation of the sequence w/ start and end labels
        $canvas->bgcolor('black');
        $canvas->rectangle(50, 200, $X - 50, 220);
        $canvas->moveTo(50, 180);
        $canvas->string("1");
        $canvas->moveTo($X - 50, 180);
        $canvas->string($sequence->get_length());

        #set scale factor. X - 100 = the number of pixels allowed for drawing
        my $scale = ($X - 100) / $sequence->get_length();


        $canvas->moveTo(20, 270);
        $canvas->string("Contigs");
        my $color = 'red';
        foreach my $contig ( @{$sequence->get_contigs()} ) {
            my $start_plot = int ($contig->get_start() * $scale) + 50;
            my $end_plot = int ($contig->get_end() * $scale) + 50;
            
            $canvas->bgcolor($color);
            $canvas->rectangle($start_plot, 250, $end_plot, 270);

            $self->draw_alignments($contig, $scale);

            $color = _color_switch($color);
        }

        open my $OUT, ">", "testpic.png";

        print $OUT $canvas->png;

        close $OUT;

    }

     sub draw_alignments {
        my ($self, $contig, $scale) = @_;

        my $canvas = $canvas_of{ident $self};
        
        $canvas->moveTo(15, 310);
        $canvas->string("Alignments");

        my $color = 'red';
        foreach my $alignment ( @{$contig->get_alignments()} ) {
            my $start_plot = int ($alignment->{'id'}->get_start() * $scale) + 50;
            my $end_plot = int ($alignment->{'id'}->get_end() * $scale) + 50;
            
            $canvas->bgcolor($color);
            $canvas->rectangle($start_plot, 290, $end_plot, 310);

            $color = _color_switch($color);
        }


    }
    sub _color_switch {

        my ($color) = @_;

        return 'red' if ( $color eq 'blue' );
        return 'blue' if ( $color eq 'red' );

    }


        



















   }

1;
