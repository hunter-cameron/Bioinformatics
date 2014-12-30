#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $tline = "";
my @tlist = ();
my $tempi = 0;
my $tempj = 0;
my $tempk = 0;

my $pline = "";
my @plist = ();

my $tid = 0;
my %tidx = ();

@ARGV==3 or die "".
"Error in syntax!\n".
"    ./genKMP.fasta.pl <input.fasta> <KMer> <output.txt>\n";

my $egFileI = $ARGV[0];
my $egKMer  = $ARGV[1];
my $egFileO = $ARGV[2];

my $egDirNow = $0; $egDirNow=~s/[\r\n]//g;
@tlist = split(/\//, $egDirNow);
$egDirNow=~s/\/$tlist[@tlist-1]$//g;

print "Formatting the sequence file ... ";
open(efO, ">$egFileI.formatted.fasta") or die "Error when formatting!\n";
my $egSeqObj = new Bio::SeqIO( -file => "$egFileI", -format => "fasta" );
my $tNumP = 0;
my %egPID2Length = ();
while( my $egSeqO = $egSeqObj->next_seq )
{
	my $egID = $egSeqO->id;
	my $egSeq = $egSeqO->seq; $egSeq=~s/[\s\d\.\-\r\n]//g;

	if( (defined $egSeq) and ($egSeq ne "") )
	{
		$egPID2Length{"$egID"} = length($egSeq);
		print efO ">".length($egSeq)."|$egID\n\U$egSeq\E\n";
		$tNumP++;
	}
}
print efO "\n";
close(efO);
print " [done] [NucleotideNumber:$tNumP]\n";

print "Generating the KMer profiles ... ";
my $egCmdLine = "$egDirNow/KMerProfiler -index \"$egFileO.K-$egKMer.kmp-idx\" -i \"$egFileI.formatted.fasta\" -o \"$egFileO\" -k $egKMer";
if( system("$egCmdLine") )
{
	die "-Error when processing this FASTA file!\n---[$egCmdLine]---\n";
}
print " [done]\n";



sub efFormatSeq
{
	my($fSeq, $tBlock, $tNumber) = @_;
	# $tBlock = 10
	# $tNumber = 5
	my $tPos = 0;
	my $tRSeq = "";
	while(1)
	{
		my $tStart = $tBlock*$tPos;
		my $tEnd = $tBlock*($tPos+1)-1;
		if( $tStart>=length($fSeq) )
		{
			last;
		}
		
		if( $tEnd>=length($fSeq) )
		{
			$tEnd = length($fSeq)-1;
		}

		$tRSeq = $tRSeq.substr($fSeq, $tStart, abs($tEnd-$tStart)+1)." ";
		if( ($tPos>0) and ($tPos%$tNumber==0) )
		{
			$tRSeq = $tRSeq."\n     ";
		}
		$tPos++;
	}
	$tRSeq=~s/^ +//g;
	$tRSeq=~s/ +$//g;
	$tRSeq = "     ".$tRSeq;
	return $tRSeq;
}


