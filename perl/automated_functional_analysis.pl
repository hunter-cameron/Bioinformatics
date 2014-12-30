#!/usr/bin.perl

###Accepts a csv file with genbank hits labeled gi|987987982| => start_pos..end_pos
##reads the first one from each line and returns information



use strict;
use warnings;

my ($file) = @ARGV;

open my $IN, "<", $file or die "Could not open $file \n";
open my $OUT, ">", "$file\_functional" or die "Could not write file\n";

while ( my $line = readline $IN ) {
    
    print "Processing line $. \n";
    last if ( $. == 1001 );
    
    chomp $line;
    
    #remove windows characters if any
    $line =~ s///g;

    if ( (my $copy = $line) =~ m {gi\|([\d]+)     #match and capture gi| and then one or more digits
                    [^=]*           #match anything after the UID but before the =>
                    =>\s*           #match => and an optional space
                    ([\d]+)..       #capture 1 or more digits (start_pos) followed by ..
                    ([\d]+)         #capture the second group of digits (end pos)
                    }xms ) {

     
        
        my ($UID, $start, $end) = ($1, $2, $3);

        #print "$UID, $start, $end\n";
        #would be nice to include the gene id for easy lookup
        #my $product = `esearch -db nucleotide -query "$UID [UID]" | efetch -format "gb" -mode "text" -seq_start $start -seq_stop $end | grep -m 1 "/product";

        #removed esearch, should be faster
        my $product = `efetch -db nucleotide -id $UID -format "gb" -mode "text" -seq_start $start -seq_stop $end | grep "/product"`;
        
        #remove leading whitespace, tag, and quotes
        chomp $product;
        $product =~ s/^.*\/product=//xms;
        $product =~ s/\n/&&/g;
        $product =~ s/"//g;

        #switch comma with - to use csv
        $product =~ s/,/-/g;

        $product = "No Results" if ( $product !~ m/\w/ );

        chomp $product;

        $line .= ",$product";
        #print "$line \n";

        }

        print $OUT $line, "\n";
}

close $IN;
close $OUT;



        

