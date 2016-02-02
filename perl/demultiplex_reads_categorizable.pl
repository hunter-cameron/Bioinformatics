#!/bin/env perl

#modified version of Sur's demultiplex_reads.pl to demultiplex all_categorizable_reads.fastq

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::SeqIO;

my $usage = "demultiplex_reads_categorizable.pl --input <fastq> --tab ..?\n";

my ($map,$input,$help);

my $opts = GetOptions("input=s" => \$input,
			"tab=s" => \$map,
			"help|map" => \$help);

#print "==$map===\n";
#print "==$input===\n";
die $usage if (! defined $input);

#consensus_demultiplex($input);
categorizable_demultiplex($input);


sub consensus_demultiplex {
    my ($input) = @_;

    my $in = Bio::SeqIO->new(-format => 'fastq', -file => "<$input");
    my @out;
    $out[0] = Bio::SeqIO->new(-format=> 'fastq', -file => ">$input\_TGA", -quality_header => 1); 
    $out[1] = Bio::SeqIO->new(-format=> 'fastq', -file => ">$input\_ACT", -quality_header => 1);

    my $count = 0;
    my $nomatch = 0;
    my $badbarcode = 0;
    print "Processing reads...\n";
    while(my $Seq = $in->next_seq){
	    my $seq = $Seq->seq;
	    my $id = $Seq->id;
	
        my ($tag) = split(";", $id);
	    my ($pnum,$mt) = split(/_/,$id);
	    my ($ft,$rt) = split(/-/,$mt);

	    #print "$id=>$ft-$rt\n";
	    if($ft =~ /([ACGT]{4,9})(TGA|ACT)([ACGT]{4})/){
		    #print "($1)($2)($3)\n";
		    my $demultiplex = $2;
		    if($demultiplex eq 'TGA'){
			    $out[0]->write_fastq($Seq);
		    }elsif($demultiplex eq 'ACT'){
			    $out[1]->write_fastq($Seq);
		    }else{
			    $badbarcode++;
		    }
	    }else{
		    $nomatch++;
	    }

	    #print "$id=>$ft-$rt\n";
        #last if ++$count >= 100;
	    $count++;
	    print "\tProcessed $count reads\n" if (($count % 50000) == 0);
    }

    close OUT;
    print "$nomatch reads out of $count did not match\n";
    print "$badbarcode reads out of $count had the wrong 3 letter barcode\n";
}

sub categorizable_demultiplex {
    my ($input) = @_;

    my $in = Bio::SeqIO->new(-format => 'fastq', -file => "<$input");
    my @out;
    $out[0] = Bio::SeqIO->new(-format=> 'fastq', -file => ">$input\_TGA", -quality_header => 1); 
    $out[1] = Bio::SeqIO->new(-format=> 'fastq', -file => ">$input\_ACT", -quality_header => 1);

    my $count = 0;
    my $nomatch = 0;
    my $badbarcode = 0;
    print "Processing reads...\n";
    while(my $Seq = $in->next_seq){
	    my $seq = $Seq->seq();
	    my $header = $Seq->desc();
        #my $id = $Seq->id()
	
        #my ($tag) = split(";", $id);
	    my ($mt) = split(" ",$header);
	    my ($ft,$rt) = split(/-/,$mt);
        #print "$header=>$mt\n"; 
        #print "$id=>$ft-$rt\n";
	    if($ft =~ /([ACGT]{4,9})(TGA|ACT)([ACGT]{4})/){
		    #print "($1)($2)($3)\n";
		    my $demultiplex = $2;
		    if($demultiplex eq 'TGA'){
			    $out[0]->write_seq($Seq);
		    }elsif($demultiplex eq 'ACT'){
			    $out[1]->write_seq($Seq);
		    }else{
			    $badbarcode++;
		    }
	    }else{
		    $nomatch++;
	    }

	    #print "$id=>$ft-$rt\n";
        #last if ++$count >= 100;
	    $count++;
	    print "\tProcessed $count reads\n" if (($count % 50000) == 0);
    }

    close OUT;
    print "$nomatch reads out of $count did not match\n";
    print "$badbarcode reads out of $count had the wrong 3 letter barcode\n";
}




__END__

sub orig_demultiplex {

    my $in = Bio::SeqIO->new(-format => 'fastq', -file => $input);
    open(OUT,'>all_demultiple_reads.fa') or die $!;
    my $count = 0;
    my $nomatch = 0;
    my $badbarcode = 0;
    print "Processing reads...\n";
    while(my $Seq = $in->next_seq){
	    my $seq = $Seq->seq;
	    my $id = $Seq->id;
	
        my ($tag) = split(";", $id);
	    my ($pnum,$mt) = split(/_/,$id);
	    my ($ft,$rt) = split(/-/,$mt);

	    #print "$id=>$ft-$rt\n";
	    if($ft =~ /([ACGT]{4,9})(TGA|ACT)([ACGT]{4})/){
		    #print "($1)($2)($3)\n";
		    my $demultiplex = $2;
		    if($demultiplex eq 'TGA'){
			    print OUT ">$Map{$pnum}->[0]_$mt\n$seq\n";
		    }elsif($demultiplex eq 'ACT'){
			    print OUT ">$Map{$pnum}->[1]_$mt\n$seq\n";
		    }else{
			    $badbarcode++;
		    }
	    }else{
		    $nomatch++;
	    }

	    #print "$id=>$ft-$rt\n";
	    #last if ++$count >= 100;
	    $count++;
	    print "\tProcessed $count reads\n" if (($count % 50000) == 0);
    }

    close OUT;
    print "$nomatch reads out of $count did not match\n";
    print "$badbarcode reads out of $count had the wrong 3 letter barcode\n";
}

sub process_map {
    my ($map) = @_;
    
    my (%Map);
    open(MAP,$map) or die "Can't open $map ($!)";
    while(<MAP>){
	    chomp;
	    my @line = split(/\t/,$_);
	    $Map{$line[0]} = [$line[3], $line[4]];
	    #print "$Map{$line[0]}->[1]\n";
    close MAP;
    print "Processed map...\n";
    return \%map;
   }




