#!/usr/bin/perl

#Takes a cBar output file and the corresponding fasta and makes a fasta out of the sequences cBar determined were plasmid sequences. 




use strict; 
use warnings;
use IO::Handle;

#intake arg for folder
(my $fasta_file, my $p_list) = @ARGV;

my $fasta_plasmid = "plasmids.fasta";

#read in all the files to an array
open(my $LIST, "<", $p_list) or die ("Could not open the list of plasmid sequences");
chomp(my @plasmids = <$LIST>);
close($LIST);

open(my $FASTA, "<", $fasta_file) or die("Could not open the fasta file.");
open(my $COMB, ">>", $fasta_plasmid); #problem if file already exists
my @headers;
my @line_nums;

#set printing to false
my $chk_print = 0;
my $pindex = 0;
#begin reading the fasta
while (my $line = <$FASTA>) {
    if ($line m/>/) {
        foreach my $plasmid (@plasmids) {
	    if (
    # if already writing stop on the >; will this write back to back plasmids?
    if ($line !~ m/$plasmids[$pindex]/ and $line =~ m/>/) {
        $chk_print = 0;
    }
    #would it write back to back if this was two ifs?
    #if matches, then start printing and increment index. This method should work since cBar keeps it in order. Would be better to rewrite it in case using a sorted cBar file. 
    elsif ($line =~ m/$plasmids[$pindex]/) {
        $chk_print = 1;
        $pindex = $pindex +1;
    }
    
    #print the line if needed
    if ($chk_print == 1) {
        print $COMB "$line";
    }



