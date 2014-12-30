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

my $egDirScript = "/home/xynlab/zhouff/work/script_app";
my $egDirGenomes = "/db/ncbi-genomes-bacteria/Bacteria092008";

my $egDate = `date \"+%F %T\"`; $egDate=~s/[\r\n]//g;

@ARGV==2 or die "".
"Error in syntax!\n".
"--------------------------------------------------------------------------\n".
"cBar: an algorithm for accurate identification of chromosomal sequence fragments from metagenome data\n".
"                                Fengfeng Zhou (FengfengZhou\@gmail.com)\n\n".
"     ./cBar.pl <input.fasta> <output.prediction.txt>\n";

my $egFileI = $ARGV[0];
my $egFileO = $ARGV[1];

my $egDirNow = $0; $egDirNow=~s/[\r\n]//g;
@tlist = split(/\//, $egDirNow);
$egDirNow=~s/\/$tlist[@tlist-1]$//g;

my $egScriptFASTA = "genKMP.fasta.pl";
my $egScriptARFF  = "genKMP.arff.pl";

my $egDirTemp = "$egDirNow/temp";

if( !(-d $egDirTemp) )
{
	die "-Error! Temporary directory does not exist!\n---[$egDirTemp]----\n";
}

my @egName = ();
my @egSeqLength = ();
my @egPrediction = ();
my $egNum = 0;

@tlist = split(/\//, $egFileI);
my $egFileTag = $tlist[@tlist-1];
$egFileTag=~s/\.[a-z]+$//g;
print "Formatting the sequences ... ";
$tempi = 0;
open(efO, ">$egDirTemp/$egFileTag.in_fasta") or die "Error when formatting!\n";
my $egSeqObj = new Bio::SeqIO( -file => "$egFileI", -format => "fasta" );
while( my $egSeqO = $egSeqObj->next_seq )
{
	my $egID = $egSeqO->id; my $egSeq = $egSeqO->seq;
	print efO ">$egID\n$egSeq\n";
	$egName[$egNum] = $egID;
	$egSeqLength[$egNum] = length($egSeq);
	$egPrediction[$egNum] = "";
	$egNum++;
	$tempi++;
}
close(efO);
print " [FormattedSeq:$tempi] [done]\n";

my $egCmdLine = "";
my $egKMer = 5;
print "Generating the K-mer profiles ... ";
$egCmdLine = "\"$egDirNow/$egScriptFASTA\" \"$egDirTemp/$egFileTag.in_fasta\" $egKMer \"$egDirTemp/$egFileTag.kmp.txt\"";
if( system("$egCmdLine > /dev/null") )
{
	die "-Error occured!\n";
}
print " [KMer]";

$egCmdLine = "\"$egDirNow/$egScriptARFF\" \"$egDirTemp/$egFileTag.kmp.txt\" /dev/null \"$egDirTemp/$egFileTag.arff\"";
if( system("$egCmdLine > /dev/null") )
{
	die "-Error occurred!\n";
}
print " [Profile]";
print " [done]\n";

my $egFileModel = "cBar.SMO.model";
print "Predicting ... ";
$egCmdLine = "java -classpath \"$egDirNow/weka.jar\" weka.classifiers.functions.SMO -l \"$egDirNow/$egFileModel\" -p 0 -T \"$egDirTemp/$egFileTag.arff\" > \"$egDirTemp/$egFileTag.dat\"";
#if( system("$egCmdLine 2>/dev/null") )
if( system("$egCmdLine") )
{
	die "-Error when predicting!\n";
}
print " [Predicted]";

open(efI, "$egDirTemp/$egFileTag.dat") or die "Error when loading the results!\n";
open(efO, ">$egFileO") or die "Error when saving!\n";
while(<efI>)
{
	if( $_=~/^\s+inst#\s+actual/ )
	{
		last;
	}
}
my $tNumChr = 0; my $tNumPld = 0;
while(<efI>)
{
	$tline = $_; $tline=~s/[\r\n]//g;
	$tline=~s/^\s+//g; $tline=~s/\s+$//g;
	@tlist = split(/\s+/, $tline);
	if( $tline=~/^\d+\s+/ )
	{
		$tid = $tlist[0]-1;
		$pline = $tlist[2]; $pline=~s/^\d+://g;
		if( (defined $tid) and ($tid>=0) and ($tid<$egNum) )
		{
			if( $pline eq "Positive" )
			{
				$egPrediction[$tid] = "Chromosome";
				$tNumChr++;
			}
			else
			{
				$egPrediction[$tid] = "Plasmid";
				$tNumPld++;
			}
		}
	}
}
close(efI);

print efO "#SeqID	Length	Prediction\n";
for($tempi=0;$tempi<$egNum;$tempi++)
{
	print efO "$egName[$tempi]	$egSeqLength[$tempi]	$egPrediction[$tempi]\n";
}
close(efO);
print " [Saved]";

if( system("rm -rf \"$egDirTemp\"/$egFileTag.*") )
{
	print " [TempNotRemoved]";
}
else
{
	print " [Cleared]";
}
print " [done]\n";
print "Prediction statistics: [TotalNum:".($tNumChr+$tNumPld)."] [Chr:$tNumChr] [Pld:$tNumPld] [done]\n";


