#!/usr/bin/perl -w

use strict;


my $tline = "";
my @tlist = ();
my $tempi = 0;
my $tempj = 0;
my $tempk = 0;

my $pline = "";
my @plist = ();

my $tid = 0;
my %tidx = ();

my $egDate = `date \"+%F %T\"`; $egDate=~s/[\r\n]//g;

@ARGV==3 or die "".
"Error in syntax!\n".
"    ./genKMP.arff.pl <input.Positive.txt> <input.Negative.txt> <output.arff>\n";

my $egFileP = $ARGV[0];
my $egFileN = $ARGV[1];
my $egFileARFF = $ARGV[2];

my $egFlagHeader = 0;
my $egDataSize = 0;
my $egTag = "";

print "Scanning the data set and generating the ARFF file ... ";
open(efO, ">$egFileARFF") or die "Error when reading!\n";
$egTag = "Positive"; $tempi = 0;
open(efI, "$egFileP") or die "Error when reading P!\n";
while(<efI>)
{
	$tline=$_; $tline=~s/[\r\n]//g; @tlist = split(/\t/, $tline);
	if( $tline=~/^#/ )
	{
	}
	elsif( @tlist >= 3 )
	{
		my($qID, $qRatio, @qList) = @tlist;

		if( $egFlagHeader==0 )
		{
			&efGenHeader(@qList);
		}

		if( (scalar @qList) == $egDataSize )
		{
			print efO "".join(",", @qList).",$egTag\n";
			$tempi++;
		}
	}
}
close(efI);
print " [P:$tempi]";

$egTag = "Negative"; $tempi = 0;
open(efI, "$egFileN") or die "Error when reading N!\n";
while(<efI>)
{
	$tline=$_; $tline=~s/[\r\n]//g; @tlist = split(/\t/, $tline);
	if( $tline=~/^#/ )
	{
	}
	elsif( @tlist >= 3 )
	{
		my($qID, $qRatio, @qList) = @tlist;

		if( $egFlagHeader==0 )
		{
			&efGenHeader(@qList);
		}

		if( (scalar @qList) == $egDataSize )
		{
			print efO "".join(",", @qList).",$egTag\n";
			$tempi++;
		}
	}
}
close(efI);
print " [N:$tempi]";

close(efO);
print " [done]\n";

sub efGenHeader
{
	my(@ttList) = @_;

	my $egName = $egFileARFF;
	$egName=~s/\s/_/g;
	print efO "\% This data set was generated using the script \"genKMP.arff.pl\" written by Dr. Fengfeng Zhou (FengfengZhou\@gmail.com).\n";
	print efO "\% Generating date: $egDate\n";
	print efO "\% Author: Fengfeng Zhou (FengfengZhou\@gmail.com)\n";
	print efO "\% Positive data file: $egFileP\n";
	print efO "\% Negative data file: $egFileN\n";
	print efO "\% \n";
	print efO "\@RELATION $egName\n";
	$egDataSize = scalar @ttList;
	for(my $tti=0;$tti<$egDataSize;$tti++)
	{
		print efO "\@ATTRIBUTE KMP_$tti real\n";
	}
	print efO "\@ATTRIBUTE VALUE { Positive, Negative }\n";
	print efO "\@DATA\n";

	$egFlagHeader = 1;

}

