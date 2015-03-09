#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use Getopt::Std;
use File::Basename;
use Cwd 'abs_path';


BEGIN {
	my $cwd = getcwd;
}




my %opts = ();
getopts("hi:c:o:s:a:", \%opts);

if ($opts{h}) {

    print "\n";
    print "\tdemultiplex.pl\n";
    print "\n\tDESCRIPTION\n";
    print "\t\tA script to perform customizable demultiplexing.\n";
    
    print "\n\tARGUMENTS\n";
    print "\t\t-i =\tNucleotide FASTA or FASTQ file of sequences to be demultiplexed. [required]\n";
    print "\t\t-o =\tOutput nucleotide FASTA or FASTQ basename of file where demultiplexed sequences should be printed. [required]\n";
    print "\t\t-c =\tConfiguration file containing demultiplexing parameters. [required]\n";
    print "\t\t-s =\tSample sheet specifying valid sample IDs. [required]\n";
    print "\t\t-a =\tAppend character that links the extracted barcode and sequence header [optional, default=\"_\"]\n";
    
    print "\n\tNOTES\n";
    print "\t\t- Assumes fasta and fastq sequences take one line each, and that each entry immediately follows the one before.\n";
    print "\t\t- Assumes the match indices for the nucleotide barcode start from 1, NOT 0.\n";
    die "\n";
}


# Read and check parameters.
my $params = read_params(\%opts);

# Update parameters to include sample sheet, which maps Sample ID to output file.
$params = read_sample_sheet($params);

# Perform the demultiplexing.
my $seq;
open my $seqFH, $params -> {inputFile};
while($seq = get_next($seqFH,$params)){
	$seq = demultiplex_sequence($seq,$params);
}
close($seqFH);

# Close all file handles.
my $ssref = $params -> {sampleSheet};
while( my( $key, $value ) = each  %$ssref){
    close($value);
}


# Subroutines.
sub read_sample_sheet {
	my $params = shift;

	# Create the sample sheet hash.
	# This hash maps SampleID to an output file handle.
	open FH, $params -> {sampleSheetFile};
	my $idSplit = $params -> {config} -> {idsplit} eq "true";
	my (@splits, %sampleSheet, %conv_code, $fh);
	open $fh, '>', $params->{outputBasename} . "_determined" . "." . $params->{fileType} unless $idSplit;
	while (<FH>) {
		chomp($_);
	    
        # HC - to write the files using Tatiana's sample ids.
        my ($id, $sample) = split "\t", $_;

		if ($idSplit) {
            open $sampleSheet{$id}, '>', $params->{outputBasename} . "_" .  $sample . "." . $params->{fileType};
            $conv_code{$id} = $sample;
            #my $out = $sampleSheet{$id};
            #print $out "# Seqs with bc $id \n";
		} else {
			$sampleSheet{$id} = $fh;
            $conv_code{$id} = $sample;
		}
	}
	
	#print "$_\n" for keys %sampleSheet;
	
	# Create file handles for undetermined and partially undetermined reads.
	open my $fh_u, '>', $params->{outputBasename} . "_undetermined". "." . $params->{fileType};
	open my $fh_pu, '>', $params->{outputBasename} . "_partially_undetermined". "." . $params->{fileType};
	$sampleSheet{undetermined} = $fh_u;
	$sampleSheet{partially_undetermined} = $fh_pu; 
	
	$params -> {sampleSheet} = \%sampleSheet;
    $params -> {conv_code} = \%conv_code;
	return($params);
	
}

sub demultiplex_sequence {
	my $seq = shift;
	my $params = shift;
	
	my (@matches, $barcode, $trimmedSeq, $sid, $qid, $qsOffset, $qsLength, @splits);
	my $regex = $params -> {config} -> {regex};
	my $bcArrRef = $params -> {config} -> {barcode};
	my $appendChar = $params -> {appendChar};
	if (@matches = ($seq->{sequence} =~ m/$regex/smg)) {
		
		## Form the barcode.
        #foreach (@$bcArrRef) {
        #	# Treat the specified barcode index
        # 	$barcode  = $barcode . $appendChar . $matches[$_-1];    # hc added append char	
        #}
	    	

        # HC Custom barcode
        my @bc;
        foreach my $indx ((0, 1, 3)) {
            if ($indx == 0 or $indx == 3) {    # these indx are frameshifts
                if (length($matches[$indx]) <= 2) {
                    push(@bc, "1,2,3");
                }
                else {
                    push(@bc, "4,5,6");
                }
            }
            else {
                push(@bc, $matches[$indx]);
            }
        }

        $barcode = join($appendChar, @bc);

		
		## Trim the sequence (and quality scores, if fastq) if requested.
		unless ($params -> {config} -> {trimto} == 0) {
			$trimmedSeq = @matches[$params->{config}->{trimto}-1];
		} else {
			$trimmedSeq = $seq -> {sequence};
		}
		
		## Update the $seq object
		# First update the sequenceID
        
		$sid = $seq -> {sequenceID};
		$seq -> {sequenceID} = substr($sid,0,1) . $barcode . $appendChar . substr($sid,1);
		
		# Next update the sequence.
		$seq -> {sequence} = $trimmedSeq;
		
		# Update fastq specific fields
		if ($params -> {fileType} eq "fastq") {
			$qid = $seq -> {qualityScoreID};
			$seq -> {qualityScoreID} = substr($qid,0,1) . $barcode . $appendChar . substr($qid,1);
			
			$qsOffset = $-[$params->{config}->{trimto}];
			$qsLength = $+[$params->{config}->{trimto}] - $qsOffset;
			$seq -> {qualityScore} = substr($seq->{qualityScore}, $qsOffset, $qsLength);
		}
		
		## Print the sequence to the appropriate file.
		# The sample ID of the sequence is everything between the defline character ('>' for fasta, '@' for fastq)
		# and the first underscore. This sample ID used to determine where the sequence will be printed.
		# This sample ID may not be one in the sample sheet, in which case this sequence will be classified
		# as partially undetermined.
		
        
        @splits = split("_", $seq->{sequenceID});
        my $sampleID = substr($splits[0],1);
        my $id_update;
        if ( exists $params->{conv_code}{$sampleID} ) {

            $id_update = $params->{conv_code}{$sampleID};
            $seq->{sequenceID} =~ s/$sampleID/$id_update/;
            $seq->{qualityScoreID} = $seq->{sequenceID};
        }

		$seq -> {sampleID} = substr($splits[0],1) if exists($params -> {sampleSheet} -> {substr($splits[0],1)});
		$seq -> {sampleID} = "partially_undetermined" unless exists($params -> {sampleSheet} -> {substr($splits[0],1)});
		
		#my $ss = $params -> {sampleSheet};
		#print "$_\n" for keys %$ss;
		
		
	} else {
		# If the sequence doesn't fit the general regex in the first place, then it's automatically undetermined.
		$seq -> {sampleID} = "undetermined";
	}
	
	print_seq($seq, $params->{sampleSheet}->{$seq->{sampleID}}, $params);
}

sub print_seq {
	my $seq = shift;
	my $fh = shift;
	my $params = shift;
	
	if (defined($fh)) {
		print $fh $seq->{sequenceID} . "\n";
		print $fh $seq->{sequence} . "\n";
		print $fh $seq->{qualityScoreID} . "\n" if $params -> {fileType} eq "fastq";
		print $fh $seq->{qualityScore} . "\n" if $params -> {fileType} eq "fastq";
	} else {
		print $seq->{sequenceID} . "\n";
		print $seq->{sequence} . "\n";
		print $seq->{qualityScoreID} . "\n" if $params -> {fileType} eq "fastq";
		print $seq->{qualityScore} . "\n" if $params -> {fileType} eq "fastq";
	}
	
}

sub get_next {
	my $fh = shift;
	my $params = shift;
	
	my ($tmp,%seq);
	$tmp = <$fh>;
	if ($tmp){
		$seq{sequenceID} = $tmp;
		$seq{sequence} = <$fh>;
		chomp($seq{sequenceID});
		chomp($seq{sequence});
		
		if ($params->{fileType} eq "fastq") {
			$seq{qualityScoreID} = <$fh>;
			$seq{qualityScore} = <$fh>;
			chomp($seq{qualityScoreID});
			chomp($seq{qualityScore});
		}
		
		return \%seq;
	} else {
		return 0;
	}	
}

sub read_params {
	my $opts = shift;
	my $inputFile = $opts -> {i};
	my $outputBasename = $opts -> {o};
	my $sampleSheet = $opts -> {s};
	my $config = read_config($opts -> {c});
	my $appendChar = $opts -> {a};
	
	## Basic parameter checks ##
	# Required parameters
	die "Input FASTA or FASTQ file required" if !defined($inputFile);
	die "Configuration file required" if !defined($config);
	die "Sample sheet required" if !defined($sampleSheet);
	die "Missing field \"regex\" in configuration file" if !defined($config->{regex});
	die "Missing field \"barcode\" in configuration file" if !defined($config->{barcode});
	die "Missing field \"trimto\" in configuration file" if !defined($config->{trimto});
	$appendChar = "_" if !defined($appendChar);
	
	# Further process config parameters
	my @bcArr = split(",", $config->{barcode});
	$config -> {barcode} = \@bcArr;
	
	
	# Wrap everything into a single hash.
	my %params;
	$params{sampleSheetFile} = $sampleSheet;
	$params{outputBasename} = $outputBasename;
	$params{inputFile} = $inputFile;
	$params{config} = $config;
	$params{fileType} = file_type($inputFile);
	$params{appendChar} = $appendChar;
	return \%params;
}

sub file_type {
	# The fastq check is very simple. If the first character of the first line is an '@' then its fastq. Otherwise we assume its fasta.
	# We will want to make this check more robust in the future.
	my $file = shift;
	open FH, $file;
	my $firstLine = <FH>;
	my $fileType = (substr($firstLine,0,1) eq '@') ? "fastq" : "fasta";
	close(FH);
	return $fileType;
}

sub read_config {
	my $configFile = shift;
	open C, $configFile;
	my (%config, @splits);
	while (<C>) {
		chomp($_);
		@splits  = split("\t", $_);
		$config{$splits[0]} = $splits[1];
	}
	close(C);	
	return \%config;
}


