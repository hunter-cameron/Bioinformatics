package SequenceMapper::Sequence_obj;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Data::Dumper;
use Class::Std::Utils;

use SequenceMapper::Alignment_obj;
use SequenceMapper::Contig_obj;

{
    #Class Attributes
    my %name_of;
    my %length_of;
    my %alignments_of;
    my %contigs_of;
    my %coverage_of;
    
    sub new {
        my ($class, $args_href) = @_;
        #print "$class\n";
        #print "$args_href->{'name'}\n";    
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{name}, $args_href->{length} ) {
            die "Cannot create new Sequence object. Parameter undefined.\n";
        }


        #bless takes two args, a reference to the variable and a string containing th ename of the class
        #/do{my $anon_scalar} = reference to a scalar that only has scope within the do statement. So the object doesn't "exist" after the end of the do but it still kept alive because of the blessed reference to it
        #could be done with the anon_scalar() method from class::std::utils
        my $new_obj = bless \do{my $anon_scalar}, $class;

        #set parameters for object; length() returns a unique id for the object
        $name_of{ident $new_obj} = $args_href->{name};
        $length_of{ident $new_obj} = $args_href->{length};
        #print "$name_of{ident $new_obj}\n"; 
        #store a reference to an empty array for the alignments hash
        $alignments_of{ident $new_obj} = []; 
        #print "success!\n"; 
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

    sub get_length {
        my ($self) = @_;
        if ( ! defined $length_of{ident $self} ) {
            return "Length of object is not defined. Make sure object was created. \n";
        }
        else {
            return $length_of{ident $self};
        }
    }

    sub add_alignment {
        my ($self, $alignment) = @_;
        #my $alignments_of; 
        #add the new alignment object reference to the array
        push @{$alignments_of{ident $self}}, $alignment;
    }


    sub get_alignments {
        my ($self) = @_;
        
        #Do I want to return an array? or the reference?
        #return the ref because it is better on memory
        return $alignments_of{ident $self};
    }

    sub get_contigs {
        my ($self) = @_;

        return $contigs_of{ident $self};
        }


    sub build_contigs {
        my ($self) = @_;
        
        #sort by alignment length
        my @sorted_aligns = sort {$a->get_length() <=> $b->get_length()}      @{$self->get_alignments};

        #set first and last positions for the sequence (use intuitive scale-starting with one)
        my @sequence;
        my $last = $self->get_length();

        #STEP 1: add each alignment giving precidence to the longer alignments (b/c of the sort)
        foreach my $alignment (@sorted_aligns) {

            for ( my $i = $alignment->get_start(); $i <= $alignment->get_end(); $i++ ) {
            
                $sequence[$i] = $alignment if ( ! defined $sequence[$i] );
                
            }
        
            print "$alignment\n";
        }
        print "seq elements =", scalar @sequence, " \n";


        #STEP 2: Parse the sequence for contigs

        #data type to hold contig alignments = anon hash within anon array
        my $contig = [];
        for ( my $i = 0; $i < @sequence; $i++ ) {
            if ( ! defined $sequence[$i] ) {
                next;
            }
            
            
            my $id = $sequence[$i];
            my $start = $i;

            while ($sequence[$i] = $sequence[$i+1]) {
                $i++;
                #print "$sequence[$i]\t$sequence[$i+1]\n";
            }

            my $end = $i;

            #add the portion of the contig to the whole thing (as a hash ref)
            push @{$contig}, {'id' => $id, 'start' => $start, 'end' => $end};
            #end the contig if there is an undefined index next
            if ( ! defined $sequence[$i+1] ) {
                $self->add_contig($contig);
                $contig = [];
            }
        }
    }

    sub add_contig {
        my ($self, $contig) = @_;
        my $new_contig = SequenceMapper::Contig_obj->new({'from_build' => $contig});

        push @{$contigs_of{ident $self}}, $new_contig;
    }

    sub perc_cov {
        my ($self) = @_;

        my $total_coverage = 0;

        foreach my $contig (@{$contigs_of{ident $self}}) {
            $total_coverage += $contig->get_length;
        }

        my $percent_cov = ($total_coverage / $self->get_length) * 100;

        return $percent_cov;
    }


    #may be handled better by a map object(that can keep track of stats, etc)
    #or perhaps each contig should be a contig object with stats
    #either way, this is very costly, don't want to allow the user to do more than once.
    sub percent_coverage {
        my ($self) = @_;
       
        #check if it has already been calculated first
        return $coverage_of{ident $self} if ( defined $coverage_of{ident $self} );



        #arrays refs to store start and end value pairs
        my @start = ();
        my @end = ();
        
        foreach my $alignment (@{$self->get_alignments()}) {
            print $alignment->get_name(), "\tstart= ", $alignment->get_start(), "\tend= ", $alignment->get_end(),"\n";
            print "Beginning ", join "\t", @start, "\n";
            print "Beginning ", join "\t", @end, "\n\n";
            

            #declare these to avoid repeated calls
            my $aln_start = $alignment->get_start();
            my $aln_end = $alignment->get_end();
           

            #this fires up only on the first alignment
            if ( ! @start ) {
                push @start, $aln_start;
                push @end, $aln_end;
                next;
            }

                
            my $cur_contig = 0;
            my $collapse = 0;
            my $new_min = $aln_start;
            my $new_max = $aln_end;
        
            for (; $cur_contig < @start; $cur_contig++ ) {
                print "$cur_contig = contig \t";
                print "$start[$cur_contig] :: $end[$cur_contig]\n";
            
                
                if ( $aln_start <= $start[$cur_contig] ) {
                
                     
                    for (my $i = $cur_contig; $i < @start; $i++ ) {
                        if ( $new_max < $start[$i] ) {
                            print "cur = $cur_contig\n";
                            print "i = $i\n";
                            last;
                        }
                        else {
                            print "cur = $cur_contig\n";
                            $new_max = $end[$i] if ( $end[$i] > $new_max );
                            $collapse++;
                        }
                    }

                    #end the loop here (keeps the contig value right) might be able to remove main loop if use a loop like this one in the second half
                    last;
                }

                elsif ( $aln_start > $start[$cur_contig] )  {

                    #next if alignment ends later than contig; if last contig, it will make a new entry
                    next if ( $aln_start > $end[$cur_contig] );
                    
                    if ( $aln_start <= $end[$cur_contig] ) {

                        $new_min = $start[$cur_contig];
                        
                        
                        for ( my $i = $cur_contig; $i < @start; $i++ ) {
                            if ( $new_max < $start[$i] ) {
                                last;
                            }
                            else {
                                $new_max = $end[$i] if ( $end[$i] > $new_max );
                                $collapse++;
                            
                            }
                            #should change loop to use cur_contig or remove the main loop instead of this hack
                        } 


                    
                    }
                    #last to end the loop -> may be able to collapse main loop
                    last;

                }
                print "end cur = $cur_contig\n";

            
            } #end contig loop


            #splice the contigs 
            print "cur = $cur_contig\n";
            splice( @start, $cur_contig, $collapse, $new_min );
            splice( @end, $cur_contig, $collapse, $new_max );
            
            print "Splice = start, ", $cur_contig, ", $collapse, $new_min\n";
            print "\nEnd ", join "\t", @start, "\n";
            print "End ", join "\t", @end, "\n\n";



        } #end alignment


        #calculate percent
        my $total_cov = 0;
        
        for ( my $i = 0; $i < @start; $i++ ) {
            $total_cov += ($end[$i] - $start[$i]);
        }
        print "$total_cov ", $self->get_length, "\n"; 
        my $perc_cov = $total_cov / $self->get_length();

        #store the value to prevent multiple calculations
        $coverage_of{ident $self} = $perc_cov;

        return $perc_cov;
        

    }

    sub Build_contigs {
        my ($self) = @_;

    }






    #delete all attributes of the object; won't destroy the object through b/c there is still ref to it in the main program?
    sub DESTROY {
        my ($self) = @_;

        delete $name_of{ident $self};
        delete $length_of{ident $self};
        delete $alignments_of{ident $self};
        delete $coverage_of{ident $self};

    }







}








1;
