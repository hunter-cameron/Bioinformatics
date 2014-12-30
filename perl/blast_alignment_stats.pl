#!/usr/bin/perl

###Gets some alignment stats for a single fasta or a folder of fastas. Alignment stats reflect similarity between query sequences. 


use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);

use Directory_obj;
use SequenceMapper::Blast_Parser_obj;
use SequenceMapper::Alignment_obj;
use SequenceMapper::Sequence_obj;

use Data::Dumper;

my $USAGE = "Usage: blast_alignment_stats.pl --db <database> --query <query.fasta>\n";


my $alignments = [];

#Set up options
my $database;
my $query;
my $min_length = 500;
my $range = 100;

GetOptions(     'db=s' => \$database,
                'query=s' => \$query,
          );

if ( any {!defined} $database, $query ) {
    die "$USAGE";
}


my $directory_obj = Directory_obj->new({'directory' => $query,});

while ( my $file = $directory_obj->next_file() ) {
    my $filename = $directory_obj->get_filename();
    main($database, $file, $filename);
    #last; 
}

print_results($alignments);

#print Dumper $s_data;


    
sub main {
    my ($db, $query, $filename) = @_;

    my $parser = SequenceMapper::Blast_Parser_obj->new({'db' => $database, 'query' => $query});

    #return if parser failed; means blast file was bad probably
    return if ( $parser == 0 );

    #return an array of all the subject objects for the file
    my @subjects = $parser->get_subjects();

    foreach my $subject (@subjects) {

        foreach my $align ( @{$subject->get_alignments()} ) {
            #remove matches below a threshhold
            next if ( $align->get_length() < $min_length );

            $alignments = process_alignment ( $alignments, $align, $filename );

        }
        #last;
    }



} #end main

sub process_alignment {

    my ($aln_ref, $align, $filename) = @_;
    
    my $name = "$filename | " . $align->get_name();
    my $q_start = $align->get_qstart();
    my $q_end = $align->get_qend();
    my $subject = $align->get_subject();
    my $s_start = $align->get_start();
    my $s_end = $align->get_end();
    
    #print Dumper $aln_ref;
    
    foreach my $existing ( @{$aln_ref} ) {
        
        #step 1: check for the same alignment
        for ( my $i = 0; $i < @{$existing->{aliases}}; $i++ ) {

            if ( $name eq $existing->{aliases}[$i]
                and $q_start ~~ [$existing->{q_starts}[$i] - $range..$existing->{q_starts}[$i] + $range]
                and $q_end ~~ [$existing->{q_ends}[$i] - $range..$existing->{q_ends}[$i] + $range] ) {


                #alignments match check if subjects do 
                for ( my $s = 0; $s < @{$existing->{subjects}}; $s++ ) {

                    if ( $subject eq $existing->{subjects}[$s]
                        and $s_start ~~ [$existing->{s_starts}[$s] - $range..$existing->{s_starts}[$s] + $range]
                        and $s_end ~~ [$existing->{s_ends}[$s] - $range..$existing->{s_starts}[$s] + $range]) {

                        #subjects match as well; do nothing, shouldn't happen often
                        #only when blast aligns two reads within $range bp of each other
                        return $aln_ref;
                    }

                    else {      #subjects do not match add a new one
                        

                        push @{$existing->{subjects}}, $subject;
                        push @{$existing->{s_starts}}, $s_start;
                        push @{$existing->{s_ends}}, $s_end;
                        return $aln_ref;
                    }
                }

            }
        }
        
        #does not match an alignment; check if it matches subject

        for ( my $s = 0; $s < @{$existing->{subjects}}; $s++ ) {

            if ( $subject eq $existing->{subjects}[$s]
                and $s_start ~~ [$existing->{s_starts}[$s] - $range..$existing->{s_starts}[$s] + $range]
                and $s_end ~~ [$existing->{s_ends}[$s] - $range..$existing->{s_starts}[$s] + $range]) {
                        
                
                #subjects match; alignments couldn't match because already checked so add new alignment

                push @{$existing->{aliases}}, $name;
                push @{$existing->{q_starts}}, $q_start;
                push @{$existing->{q_ends}}, $q_end;
                return $aln_ref;
            }
        }
    }

    #if it gets here, it has no matches so add a new alignment
    push @{$aln_ref}, {   aliases =>  [$name],
                            q_starts => [$q_start],
                            q_ends =>   [$q_end],
                            subjects => [$subject],
                            s_starts => [$s_start],
                            s_ends =>   [$s_end],
                        };

    return $aln_ref;

}

         
sub print_results {
    my ($aln_ref) = @_;

    open my $OUT, ">", "alignment_stats.csv";


    print $OUT "Alignment ID,# Subjects Matched,# Contig Ranges Matched,Subjects (name => start..end),Contig Ranges(file | contig => start..end)\n";


    my $aln_id = 1;
    foreach my $alignment (@{$aln_ref}) {

        #build subjects field
        my @subjects;
        for ( my $i = 0; $i < @{$alignment->{subjects}}; $i++ ) {
            
            push @subjects, "$alignment->{subjects}[$i] => $alignment->{s_starts}[$i]..$alignment->{s_ends}[$i]";
        }

        #build contig ranges field
        my @contig_ranges;
        for ( my $i = 0; $i < @{$alignment->{aliases}}; $i++ ) {

            push @contig_ranges, "$alignment->{aliases}[$i] => $alignment->{q_starts}[$i]..$alignment->{q_ends}[$i]";
        }

        #print results
        print $OUT join ",", (   "Alignment $aln_id",
                                scalar @{$alignment->{subjects}},
                                scalar @{$alignment->{aliases}},
                                join( "&&", @subjects),
                                join( "&&", @contig_ranges),
                            ), "\n";
        $aln_id++;

    }
}
            
            
