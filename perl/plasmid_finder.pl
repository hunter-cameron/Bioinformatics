#!/usr/bin/perl

###A new algorithm to find plasmids that takes into account cBar, length, coverage, and homology with other plasmids. 

#perhaps it shouldn't take into account homology unitl it is clear that chromosomal sequences don''t share great homology with them. 

#Possibly also use homology between files? Ie. homology to high schoring plasmids

#also perhaps use kmer profiles within the sequences - expect plasmids to be different than chromosomes


use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);
use Statistics::Basic qw(stddev mean);
use Cwd qw(getcwd);

use BioUtils::FastaIO;
use BioUtils::FastaSeq;


use Directory_obj;
use SequenceMapper::Blast_Parser_obj;
use SequenceMapper::Alignment_obj;
use SequenceMapper::Sequence_obj;

use Data::Dumper;

my $USAGE = "Usage: plasmid_finder.pl -db <database> -query <query.fasta>\n";


#set options for parameter scoring
my $CBAR_VALUE = 10;        #cBar plasmid
my $COVERAGE_VALUE = 10;    #each stdev above mean
my $LENGTH_VALUE = 1;       #which percent of length does it fit in? by 20%; eg: top 20% = 5 * $LENGTH_VAULE
                            #below 1000 = 0
my $CHROMOSOME_VALUE = -5;  #95% cov by and 90% ident to a single sequence in NT (most likely a chromosome)
my $THRESHHOLD = 13;        #minimum score to be classified as plasmid



#Set up options
my $database;
my $query;

GetOptions(     'db=s' => \$database,
                'query=s' => \$query,
          );


#die if undefined
if ( any {!defined} $database, $query ) {
    die "$USAGE";
}



#shell to hold info for multiple files
my $DATA = {};


#make a directory object
my $directory_obj = Directory_obj->new({'directory' => $query,});


#loop through all the files in the directory (or just the single file)
while ( my $file = $directory_obj->next_file() ) {

    #get the filename (no path or ext)
    my $filename = $directory_obj->get_filename();

    #call the main program
    main($database, $file, $filename);
}







sub main {

    my ($db, $file, $filename) = @_;

    #make a new reference for convenience
    my $file_ref = $DATA->{$filename}; 


    #make BioUtils object to read in sequences 
    my $in = BioUtils::FastaIO->new({stream_type => '<', file => $file});

    #loop through all sequences and get the length
    while ( my $seq_obj = $in->get_next_seq() ) {
        my $header = $seq_obj->get_header();
        my $seq = $seq_obj->get_seq();

        $file_ref->{$header}{length} = length $seq;

        }


    #calls to the subs to score different aspects

    $file_ref = score_length($file_ref);
    $file_ref = score_coverage($file_ref);
    $file_ref = score_cbar_result($file, $file_ref);
    $file_ref = score_chromosome_blast($file, $file_ref);   #takes vastly longer than the other steps; can be omitted
    #score_blast_homology($db, $file, $file_ref);



    #make predictions by adding up scores and print results
    $file_ref = make_predictions($file_ref);
    print_table($filename, $file_ref);

}





sub make_predictions {

    my ($file_ref) = @_;
    
    #loop through each contig
    foreach my $contig ( keys %{$file_ref} ) {

        #make contig a reference
        $contig = $file_ref->{$contig};


        #calculate score
        my $score = 0;
        $score += $contig->{cbar_score} * $CBAR_VALUE;
        $score += $contig->{coverage_score} * $COVERAGE_VALUE;
        $score += $contig->{length_score} * $LENGTH_VALUE;
        $score += $contig->{chromosome_score} * $CHROMOSOME_VALUE;  
        $contig->{score} = $score;

        #make predictions
        if ( $score >= $THRESHHOLD ) {
            $contig->{prediction} = "Plasmid";
        }
        else {
            $contig->{prediction} = "Chromosome";
        }
    }

    #return additions to hash
    return $file_ref;
}




#Scores the contigs based on their relative lengths (relative to account for small genomes/poor contig length) 

sub score_length {
    my ( $file_ref ) = @_;

    #counter to keep track of place to determine which percentage 
    my $counter = 1;

    #loop through the contigs sorted by length
    foreach my $contig ( sort {$file_ref->{$a}{length} <=> $file_ref->{$b}{length}} keys %{$file_ref} ) {
        my ($total_contigs) = scalar keys %{$file_ref};

        #change the string to a proper reference. 
        $contig = $file_ref->{$contig};

        

        #if contig isn't at least 1000 bp long, give it a score of 0
        if ( $contig->{length} <= 1000 ) {
            $contig->{length_score} = 0;

        }
        elsif ( int ($counter / $total_contigs * 100) ~~ [0..20]) {      #bottom 20% gets a 1 and so on...
            $contig->{length_score} = 1;
        }
        elsif ( int ($counter / $total_contigs * 100) ~~ [21..40]) {
            $contig->{length_score} = 2;
        }
        elsif ( int ($counter / $total_contigs * 100) ~~ [41..60]) {
            $contig->{length_score} = 3;
        }
        elsif ( int ($counter / $total_contigs * 100) ~~ [61..80]) {
            $contig->{length_score} = 4;
        }
        elsif ( int ($counter / $total_contigs * 100) ~~ [81..100]) {
            $contig->{length_score} = 5;
        }
        
        $counter++;
    }

    #return additions to the hash
    return $file_ref;
}




#scores contigs based on standard deviations from the mean of coverage

sub score_coverage {

    my ($file_ref) = @_;

    #an array to hold the contigs
    my @coverages;
    
    #parse the coverage out of the contig name
    foreach my $contig ( keys %{$file_ref} ) {

        my $ref = $file_ref->{$contig};

        #match the _cov_ number
        my $coverage;
        if ( $contig =~ m/_cov_([\d.]*)/ ) {     #match _cov_ and then some digits or a decimal
            $coverage = $1;
        }

        if ( ! defined $coverage ) {

            die "Warning: coverage for $contig not found \n";
        }

        $ref->{coverage} = $coverage;
        push @coverages, $coverage;
    }



    #calculate mean and stddev
    my $mean = mean(@coverages);
    my $stddev = stddev(@coverages);
    

    #calculate the score. score is a measure of deviations from the mean. 
    foreach my $contig ( keys %{$file_ref} ) {

        #turn contig into a ref
        my $contig = $file_ref->{$contig};
        
        my $score = (($contig->{coverage} - $mean) / $stddev); 
        
        #adds either a .5 or a -.5 depending on whether score is positive of neagive to round rather than floor
        $score = int ( $score + $score/abs($score*2) ); 
        
        $contig->{coverage_score} = $score;
    
    }

    #returns the additions to the hash
    return $file_ref;
}




#scores the contigs based on whether or not cBar says they are a plasmid (~70 acc?)

sub score_cbar_result {
    my ($file, $file_ref) = @_;

    my $tempout = "cBar_temp.txt";

    #call cBar - needs to be updated to be used other than from my home
    system( "perl ~/scripts/plasgraph/cBar.1.2/cBar.pl $file $tempout" ) == 0 or die "Could not run cBar\n";

    #read in the cBar output
    open my $IN, "<", $tempout;

    while ( my $line = readline $IN ) {


        #skip header line
        next if ($. == 1);

        chomp $line;

        #parse out cBar's prediction
        my ($contig, $length, $prediction) = split "\t", $line;

        if ( $prediction eq "Plasmid" ) {
            $file_ref->{$contig}{cbar_score} = 1;
        }
        elsif ( $prediction eq "Chromosome" ) {
            $file_ref->{$contig}{cbar_score} = 0;
        }
    }

    close $IN;

    #return additons to the hash
    return $file_ref;
}




#blasts to the NT database and classifies things with 95% coverage and 90% identity as chromosomes - it would be better to use a chromosome-only database, but our plasmids (so far) have not had high identity to anything - which might make sense, at least until more plasmids get sequenced.
        
sub score_chromosome_blast {
    my ($file, $file_ref) = @_;
   
    #get working path
    my $path = getcwd();

    #make an output for the chromosome_reversal program
    my $tempout = "$path/chromosome.tempout";

    system( "chromosome_reversal_forked.pl -out $tempout -db /nas02/data/blast/nt -query $file") == 0 or die "Could not call chromosomal_reversal.\n";

    #read output
    open my $IN, "<", $tempout or die "Could not open $tempout for reading\n";

    while ( my $line = readline $IN ) {
        chomp $line;

        #skip header
        next if ( $. == 1 );
        my ($contig, $coverage, $identity, $prediction ) = split "\t", $line;

        if ( $prediction eq "Chromosome" ) {
            $file_ref->{$contig}{chromosome_score} = 1;
        }
        if ( $prediction eq "Plasmid" ) {
            $file_ref->{$contig}{chromosome_score} = 0;
        }
    }

    close $IN;

    #return additions to the hash
    return $file_ref;

}




#prints the results 

sub print_table {

    my ($filename, $file_ref) = @_;

    open my $OUT, ">", "$filename\_table.tab";
    
    #print header
    print $OUT "Contig\tPrediction\tScore\tLength Score\tCoverage Score\tcBar Plasmid Prediction\tChromosome Prediction\n";

    #loop through contigs and print info
    foreach my $contig ( keys %{$file_ref} ) {

        my $ref = $file_ref->{$contig};

        print $OUT join "\t", ( $contig,
                                $ref->{prediction},
                                $ref->{score},
                                $ref->{length_score},
                                $ref->{coverage_score},
                                $ref->{cbar_score},
                                $ref->{chromosome_score},
                                ), "\n";

    }
}


__END__

    
sub score_blast_homology {
    my ($db, $file, $filename) = @_;

    my $parser = SequenceMapper::Blast_Parser_obj->new({'db' => $db, 'query' => $file});

    #return if parser failed; means blast file was bad probably
    return if ( $parser == 0 );

    #return an array of all the subject objects for the file
    my @subjects = $parser->get_subjects();

    foreach my $subject (@subjects) {

        foreach my $align ( @{$subject->get_alignments()} ) {
            #remove matches below a threshhold
            next if ( $align->get_length() < $min_length );

            #get name and increment

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
            
            
