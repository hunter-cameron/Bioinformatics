#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);

use Directory_obj;
use SequenceMapper::Alignment_obj;
use SequenceMapper::Sequence_obj;
use SequenceMapper::Drawing_obj;

my $USAGE = "Usage: map_plasmid_by_contig.pl --db <database> --query <query.fasta>\n";


#Set up options
my $database;
my $query;


my @subjects;

GetOptions(     'db=s' => \$database,
                'query=s' => \$query,
          );

if ( any {!defined} $database, $query ) {
    die "$USAGE";
}


my $directory_obj = Directory_obj->new({'directory' => $query,});

while ( my $file = $directory_obj->next_file() ) {
    my $filename = $directory_obj->get_filename();
    main($file, $filename);

    #clean up at the end or each loop
    foreach my $subject (@subjects) {
        $subject->DESTROY();
    }

    @subjects = ();
    unlink "temp_results.out";
}



    
sub main {
    my ($query, $filename) = @_;
    print "Query = $query\n";
    #call to blast, die if unsuccessful
    system( "blastn -db $database -query $query -out 'temp_results.out' -outfmt '6 std stitle'" ) == 0 
        or die "Blast could not be done successfully\n";


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
        check_subj(\%ln_elements);
    
        #last;
    }

    close $IN;
    print_results($filename);    
   

} #end main

sub print_results {
    my ($filename) = @_;
    
    my %cov;
    open my $OUT, ">", "$filename.txt";
    
    foreach my $subject (@subjects) {
        $subject->build_contigs();
        
    }

   
    #gi|410692015|ref|NC_019324.1| is different between the two algorithms, build has a greater coverage
    #sort the output in order of greatest coverage
    foreach my $subject (sort {$b->perc_cov() <=> $a->perc_cov()} @subjects ) {
        $cov{$subject} = $subject->perc_cov();
        print $OUT $subject->get_name(), "\t Cov = $cov{$subject} \n";
        foreach my $alignment ( @{$subject->get_alignments()} ) {
            #print $OUT "\t" . $alignment->get_name() . "\n";
        }
        my $png = SequenceMapper::Drawing_obj->new();

        $png->draw_contigs($subject);

        last; #print only the best complete w/ alignments.
    }

    close $OUT;
}

sub get_slength {
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

sub check_subj {
    my ($args_href) = @_;
    #print "$args_href->{'sseqid'}\n";
    foreach my $subj (@subjects) {
        #check if the subject already exists, if so add the reference
        #print $subj->get_name() . " <> " . $args_href->{'sseqid'} . "\n";
        if ( $subj->get_name() eq $args_href->{'sseqid'} ) {
            my $alignment = new_alignment($args_href);
            $subj->add_alignment($alignment);
            return;
        }
    }
    #print "Making new subject\n";
    #if the subject doesn't exist, make it and then add the reference
    my $subj = new_subj($args_href);
    my $alignment = new_alignment($args_href);
    $subj->add_alignment($alignment);

    #add the subject to the array
    push @subjects, $subj;

    return;
}

sub new_subj {
    my ($args_href) = @_;

    #get the length of the subject
    my $slength = get_slength($args_href);
    #print "$args_href->{'sseqid'}\n";
    #print "$slength\n";
    #make a new sequence object
    my $subj = SequenceMapper::Sequence_obj->new({'name' => $args_href->{'sseqid'}, 'length' => $slength,});

    return $subj;
}
 
sub new_alignment {
    my ($args_href) = @_;
    
    #make a new alignment object
    my $alignment = SequenceMapper::Alignment_obj->new({'name' => $args_href->{'qseqid'}, 's_start' => $args_href->{'sstart'}, 's_end' => $args_href->{'send'}, 'perc_iden' => $args_href->{'pident'}});

    return $alignment;
}

