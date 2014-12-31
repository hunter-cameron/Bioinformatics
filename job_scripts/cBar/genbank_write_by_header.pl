#!/usr/bin/perl

use strict;
use warnings;

use MyLibs::GenbankTools;
use Directory_obj;


my $usage = "Usage: genbank_write_by_header.pl <folder of genbanks> <folder of cBar results> \n";

@ARGV == 2 or die "Wrong \n";

my $directory_gbk = Directory_obj->new({'directory' => $ARGV[0]});

while ( my $gbk = $directory_gbk->next_file() ) {
    my $filename = $directory_gbk->get_filename();
    print "gbk = $filename\n";
    my $directory_cbar = Directory_obj->new({'directory' => $ARGV[1]});
    while ( my $cbar = $directory_cbar->next_file() ) {
        my $cbar_file = $directory_cbar->get_filename();
        print "cbar = $cbar_file\n";
        if ( $cbar_file =~ m/$filename/ ) {
            print "match = $cbar_file\n";
            main($gbk, $cbar);
            last;


        }

    

    }
    #last;
}






sub main {
    my ($genbank, $cbar) = @_;
    print "$genbank, \t $cbar\n"; 
    `grep "Plasmid" $cbar | awk '{print \$1}' > pheader.temp`;




    #call GenbankTools
    MyLibs::GenbankTools::write_by_header($genbank, "pheader.temp");
}

