package SequenceMapper::Blast_Parser_obj;


use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);
use Class::Std::Utils;

use Directory_obj;
use SequenceMapper::Alignment_obj;
use SequenceMapper::Sequence_obj;
use SequenceMapper::Drawing_obj;

use Data::Dumper;


{

    #class attributesi
    my %db_of;
    my %query_of;
    my %subjects_of;

    sub new {
        my ($class, $args_href) = @_;
        #print "$class\n";
        #print "$args_href->{'name'}\n";    
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{db}, $args_href->{query} ) {
            die "Cannot create new Blast Parser object. Parameter undefined.\n";
        }

        #$bless a new object
        my $new_obj = bless \do{my $anon_scalar}, $class;

        #set parameters for object; length() returns a unique id for the object
        $db_of{ident $new_obj} = $args_href->{db};
        $query_of{ident $new_obj} = $args_href->{query};
        $subjects_of{ident $new_obj} = [];
        
        #perform the BLAST and read in the subjects, return 0 if failed. 
        if ( $new_obj->_parse_by_subject() ) {
        return $new_obj;
        }
        else {
            return 0;
        }
    }

    sub get_subjects {
        my ($self) = @_;
        #return an actual array here for easy use in for loops. Is this really what I want to do?
        return @{$subjects_of{ident $self}};
    }

    sub add_subject {
        my ($self, $subject) = @_;
        #print "Subj to add $subject\n";
        #no protection against users adding subjects (no validation)
        push @{$subjects_of{ident $self}}, $subject;
    }
    
    sub _parse_by_subject {
        my ($self) = @_;
    
        my $database = $db_of{ident $self};
        my $query = $query_of{ident $self};
    
        #call to blast, return failed if unsuccessful
        if ( system( "blastn -db $database -query $query -out 'temp_results.out' -outfmt '6 std stitle'" ) != 0 ) {
            print "Blast could not be done successfully:\n blastn -db $database -query $query -out 'temp_results.out' -outfmt '6 std stitle'\n";
            return 0;
        }

        #process each line of the blast results
        open my $IN, "<", 'temp_results.out';
        while ( my $line = readline $IN ) {
            chomp $line;
    
            #print "$.\n";
            my @ln_elements = split "\t", $line;

            my %ln_elements = ( 'qseqid'    =>  $ln_elements[0],
                                'sseqid'    =>  $ln_elements[1],
                                'pident'    =>  $ln_elements[2],
                                'length'    =>  $ln_elements[3],
                                'mismatch'  =>  $ln_elements[4],
                                'gapopen'   =>  $ln_elements[5],
                                'qstart'    =>  $ln_elements[6],
                                'qend'      =>  $ln_elements[7],
                                'sstart'    =>  $ln_elements[8],
                                'send'      =>  $ln_elements[9],
                                'evalue'    =>  $ln_elements[10],
                                'bitscore'  =>  $ln_elements[11],
                                'stitle'    =>  $ln_elements[12],
                                );

            #print "$ln_elements{'sseqid'}\n";
            _check_subj($self, \%ln_elements);
    
            #last;
        }

        #close filehandle and delete the temporary file
        close $IN;
        unlink 'temp_results.out';

    }



    sub _get_slength {
        my ($args_href) = @_;
        #print "$args_href->{'stitle'}\n";
        my $slength;
        #match len_ and some digits
        if ( $args_href->{'stitle'} =~ m/len_([\d]+)/ ) {
        $slength = $1;
        }
        else {
            die "Couldn't get the length of the subject";
        }
    
        return $slength;
    }

    sub _check_subj {
        my ($self, $args_href) = @_;
        #print "Checking Subjects...\n";
        
        foreach my $subj ($self->get_subjects()) {
            #print Dumper $subj;
            #check if the subject already exists, if so add the reference
            #print $subj->get_name() . " <> " . $args_href->{'sseqid'} . "\n";
            if ( $subj->get_name() eq $args_href->{'sseqid'} ) {
                my $alignment = _new_alignment($args_href, $subj);
                $subj->add_alignment($alignment);
                return;
            }
        }

        #if the subject doesn't exist, make it and then add the reference
        my $subj = _new_subj($args_href);
        my $alignment = _new_alignment($args_href);
        $subj->add_alignment($alignment);
        
        #print "Here is the subject $subj\n";
        #add the subject to the array
        $self->add_subject($subj);

        return;
    }

    sub _new_subj {
        my ($args_href) = @_;

        #get the length of the subject
        my $slength = _get_slength($args_href);
        #print "$args_href->{'sseqid'}\n";
        #print "$slength\n";
        #make a new sequence object
        my $subj = SequenceMapper::Sequence_obj->new({'name' => $args_href->{'sseqid'}, 'length' => $slength,});

        return $subj;
    }
 
    sub _new_alignment {
        my ($args_href, $subject) = @_;
    
        #make a new alignment object
        my $alignment = SequenceMapper::Alignment_obj->new({    'name' => $args_href->{'qseqid'}, 
                                                                's_start' => $args_href->{'sstart'}, 
                                                                's_end' => $args_href->{'send'},
                                                                'q_start' => $args_href->{'qstart'},
                                                                'q_end' => $args_href->{'qend'},
                                                                'perc_iden' => $args_href->{'pident'}, 
                                                                'subj' => $args_href->{'sseqid'},
                                                                });

        return $alignment;
    }


    sub DESTROY {
        my ($self) = @_;

        delete $db_of{ident $self};
        delete $query_of{ident $self};
        delete $subjects_of{ident $self};

    }

}

1;

