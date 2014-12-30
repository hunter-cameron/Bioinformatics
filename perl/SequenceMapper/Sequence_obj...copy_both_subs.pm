package SequenceMapper::Sequence_obj;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;

use SequenceMapper::Alignment_obj;

{
    #Class Attributes
    my %name_of;
    my %length_of;
    my %alignments_of;
    
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


    

    

    #may be handled better by a map object(that can keep track of stats, etc)
    #or perhaps each contig should be a contig object with stats
    #either way, this is very costly, don't want to allow the user to do more than once.
    sub percent_coverage {
        my ($self) = @_;
        
        #arrays to store start and end value pairs
        my @start = ();
        my @end = ();
        
        #my $tester = 0; #works with only one iteration (b/c of the initial thing)
        foreach my $alignment (@{$self->get_alignments()}) {
            print $alignment->get_name(), "\tstart= ", $alignment->get_start(), "\tend= ", $alignment->get_end(),"\n";
            print "Beginning ", join "\t", @start, "\n";
            print "Beginning ", join "\t", @end, "\n\n";
            #last if ($tester == 2);
            #$tester++;
            
            #print $alignment->get_name(), "\n";
            #next;

            #declare these to avoid repeated calls
            my $aln_start = $alignment->get_start();
            my $aln_end = $alignment->get_end();
            
            #this fires up only on the first alignment
            if ( ! @start ) {
                print "$aln_start to start\n";
                push @start, $aln_start;
                push @end, $aln_end;
                next;
            }
            else {
                print "Running for loop\n";
            }

                

            
            for ( my $contig = 0; $contig < @start; $contig++ ) {
                #print $alignment->get_name(), "\t";
                print "$contig = contig \t";
                print "$start[$contig] :: $end[$contig]\n";
                #
                ##if start is before the contig
                #
               
                if ( $aln_start < $start[$contig] ) {
                #if ( $aln_start < $start[$contig] ) {
                    #if the end is less than the start of the contig, add a new contig before that one
                    #(alignment is a separate contig before the current first)
                    if ( $aln_end < $start[$contig]) {
                        print "Here 1\n";
                        splice @start, $contig, 0, $aln_start;
                        splice @end, $contig, 0 ,$aln_end;
                    }

                    #if the end is encompassed in the contig; update the start
                    #(alignment extended the start of the contig)
                    elsif ( $aln_end <= $end[$contig] ) {
                        print "Here 2 \n";
                        $start[$contig] = $aln_start;
                    }

                    #if the end is greater than the contig; update the start and the end
                    #(contig is encompassed within the alignment)
                    elsif ( $aln_end > $end[$contig] ) {
                        print "Here 3\n";
                        $start[$contig] = $aln_start;
                        $end[$contig] = $aln_end;

                        #better way to do this? work with refs the whole time?
                        my ($ref_start, $ref_end) = check_collapse($aln_start, $aln_end, \@start, \@end, $contig);
                        @start = @{$ref_start};
                        @end = @{$ref_end};
                    }

                }

                #
                ##if the start is in the middle of the contig
                #
                elsif ( $aln_start >= $start[$contig] and $aln_start <= $end[$contig] ) {
                
                    #if the end is larger, update the end
                    #(alignment extended the end of the contig)
                    if ( $aln_end > $end[$contig] ) {
                        print "Here 4\n";
                        #$end[$contig] = $aln_end;
    
                        #better way to do this? work with refs the whole time?
                        my ($ref_start, $ref_end) = check_collapse($aln_start, $aln_end, \@start, \@end, $contig);
                        @start = @{$ref_start};
                        @end = @{$ref_end};
    
                    }

                }
             
                #
                ##if the start is after the contig
                #
                elsif ( $aln_start > $end[$contig] ) {
                    #check if this is the last contig; 
                    if ( $contig == @start - 1 ) {  #if so, make this alignment a new contig
                        print "Here 5\n";
                        push @start, $aln_start;
                        push @end, $aln_end;
                    }
                    else {      #otherwise go to the next contig
                        print "Here 6\n";
                        next;
                    }
                }

                #if it makes it here (i.e. did not get caught by  the next loop above, then alignment was mapped so go to next
                print "Here end\n";
                last;
                
            } #end for contigs
            
                print "\nEnd ", join "\t", @start, "\n";
                print "End ", join "\t", @end, "\n\n";
        } #end foreach alignment
        

        #calculate percent
        my $total_cov = 0;
        
        for ( my $i = 0; $i < @start; $i++ ) {
            $total_cov += ($end[$i] - $start[$i]);
        }
        print "$total_cov ", $self->get_length, "\n"; 
        my $perc_cov = $total_cov / $self->get_length();

        return $perc_cov;
    }
    
    #could I do all this with splice? the whole alignment assembly? Much fewer if statements?
    #basically, this sub with a couple extra qualifications could do all the work of the nested ifs; just replace every time in this fashion.
    sub check_collapse {
        my ($aln_start, $aln_end, $start_ref, $end_ref, $cur_contig) = @_;
        
        my $collapse_indx = 0;
        my $new_min = $aln_end;
        my $new_max = $aln_end;
        
        
        #start loop at the next contig
        for ( my $i = $cur_contig; $i < @{$start_ref}; $i++ ) {
            #no need to collapse if the end is less than the start of the next
            $new_min = $start_ref->[$i] if ( $start_ref->[$i] < $new_min );

            if ( $aln_end < $start_ref->[$i] ) {
                last;
            }
            #if end is greater than or equal to the start of the next, need to collapse
            elsif ( $aln_end >= $start_ref->[$i] ) {
                $collapse_indx += 1;
                $new_max = $end_ref->[$i] if ( $end_ref->[$i] > $new_max );
            }

        }

        #collapse contigs if necessary
        #collapse current contig and all others that overlap and replace with new contig entry
 
        #splice the contigs correctly
        splice( @{$start_ref}, $cur_contig, $collapse_indx, $new_min );
        splice( @{$end_ref}, $cur_contig, $collapse_indx, $new_max );
        
        #return a new array (though with the same name)
        return $start_ref, $end_ref;
        

    }

    #could I do all this with splice? the whole alignment assembly? Much fewer if statements?
    #basically, this sub with a couple extra qualifications could do all the work of the nested ifs; just replace every time in this fashion.
    sub new_perc {
        my ($self) = @_;
       
        #arrays refs to store start and end value pairs
        my @start = ();
        my @end = ();
        
        #my $tester = 0; #works with only one iteration (b/c of the initial thing)
        foreach my $alignment (@{$self->get_alignments()}) {
            print $alignment->get_name(), "\tstart= ", $alignment->get_start(), "\tend= ", $alignment->get_end(),"\n";
            print "Beginning ", join "\t", @start, "\n";
            print "Beginning ", join "\t", @end, "\n\n";
            #last if ($tester == 2);
            #$tester++;
            
            #print $alignment->get_name(), "\n";
            #next;

            #declare these to avoid repeated calls
            my $aln_start = $alignment->get_start();
            my $aln_end = $alignment->get_end();
           
          
                    
   

            #this fires up only on the first alignment
            if ( ! @start ) {
                push @start, $aln_start;
                push @end, $aln_end;
                next;
            }
            else {
                print "Running for loop\n";
            }


                
            #start loop at the next contig
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
                            print "cur = $cur_contig\n";

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
                            #$cur_contig = $i;
                        } 


                    
                    }
                    #last to end the loop -> may be able to collapse main loop
                    last;

                }
                print "end cur = $cur_contig\n";

            
            } #end contig loop


            #splice the contigs correctly
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

        return $perc_cov;
        

    }









}








1;
