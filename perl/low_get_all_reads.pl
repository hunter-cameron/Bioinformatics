#!/usr/bin/perl

###Retrieves raw reads from JGI for genbank files based on genone name 


use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);

use Directory_obj;

use Data::Dumper;

my $USAGE = "Usage: retrieve_raw_reads.pl <folder/file of genbank>\n";


#Set up options

my $genbank;
my $XML;

my $cookie;
GetOptions( 
                'genbank=s' => \$genbank,
                'xml=s' => \$XML,
                'c=s' => \$cookie,
          );

if ( !defined $genbank and !defined $XML  ) {
    die "$USAGE";
}

if ( defined $XML ) {
    use_XML($XML);
}

else {
    #directory object will loop through a directory of files or execute on a lone file; depending on input
    my $directory_obj = Directory_obj->new({directory => $genbank, ext => '.gbk'});


    while ( my $file = $directory_obj->next_file() ) {
        my $filename = $directory_obj->get_filename();
        use_genbank($file);

    }
}



sub use_XML {
    my ($file) = @_;

    open my $XML, "<", $file or die "Couldn't open file\n";

    my %downloads;

    while ( my $line = readline $XML ) {
        chomp $line;
        if ( $line =~ m/<folder name="Raw Data">/ ) {
            
            while ( (my $sub_line = readline $XML) !~ m/<\/folder>/ ) {
                
                my ($key, $url) = "Blank";
                #parse out the name
                $key = $1 if ( $sub_line =~  m{label="       #match label="
                                                ([^"]+)            #capture everything but "...
                                                "               #...until the closing quotes
                                                }xms );
                #parse out the url
                $url = $1 if ($sub_line =~ m{    url="       #match url="
                                                    ([^"]+)        #capture everything but "...
                                                    "           #...until the closing quotes
                                                    }xms );

                #add the information to the hash
                $downloads{$key} = $url;
            }
        }
    }

    #if no links are found, die
    die "Error reading XML file! \n" if ( ! %downloads );

     
    #print Dumper %downloads;
    download(\%downloads);
}

sub download {
    my ($data) = @_;

    foreach my $key ( keys %{$data} ) {

        my $name = $key;
        my $url = $data->{$key};

        open my $OUT, ">", "$name.fastq.gz";
        print "curl \"http://genome.jgi.doe.gov$url\" -c $cookie -b $cookie >test.fast\n";

        `curl "http://genome.jgi.doe.gov$url" -c $cookie -b $cookie > test.fast`;


        #print $OUT system("curl", "http://genome.jgi.doe.gov$url", "-c $cookie", "-b $cookie");
        last; 
    }
}





    
