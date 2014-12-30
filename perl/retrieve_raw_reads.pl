#!/usr/bin/perl

###Retrieves raw reads from JGI for genbank files based on genone name 


use strict;
use warnings;

use Getopt::Long;
use List::MoreUtils qw(any);

use Directory_obj;
#use threads;
#use Thread::Queue;
#use LWP::UserAgent;

use Data::Dumper;

my $USAGE = "Usage: retrieve_raw_reads.pl <folder/file of genbank>\n";


#Set up options

my $genbank;
my $XML;
my $CPU = 16;   #use default 16 threads

my $cookie;
GetOptions( 
                'genbank=s' => \$genbank,
                'xml=s' => \$XML,
                'cpu=s' => \$CPU,
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
        if ( $line eq '<folder name="Raw Data">' ) {
            
            while ( (my $sub_line = readline $XML) ne '</folder>' ) {

                #parse out the name
                my $key = $1 if ( $sub_line =~  m{label="       #match label="
                                                (.*)            #capture everything...
                                                "               #...until the closing quotes
                                                }xms );
                #parse out the url
                my $url = $1 if ($sub_line =~ m{    url="       #match url="
                                                    (.*)        #capture everything...
                                                    "           #...until the closing quotes
                                                    }xms );

                #add the information to the hash
                $downloads{$key} = $url;
            }
        }
    }

    #if no links are found, die
    die "Error reading XML file! \n" if ( ! %downloads );

     
    print Dumper %downloads; 
    #run_threaded(\%downloads);
}

__END__

sub run_threaded {
    my ($data) = @_;
    
    #create a queue
    my $Q = Thread::Queue->new();

    #load the que with the download info
    foreach my $key ( keys %{$data} ) {
        $Q->enqueue({'name' => $key, 'url' => $data->{$key}});
    }

    #initialize threads
    my @threads = init_threads();

    #start each thread
    foreach (@threads) {
        #make a thread and pass the download subroutine
        $_ = $threads->create(\&download);

    }
    
    #signal the que is complete   
    $Q->end();

    #join when the threads are done
    foreach (@threads) {
        $_->join();
    }
}

#initialize an array to hold the threads
sub init_threads {
    my @init_threads;

    for ( my $i = 0; $i <= $CPU; $i++ ) {
        push @init_threads, $i;
    }

    return @init_threads;
}

}

sub download {

    while ( defined my $href = $Q->dequeue()) {
        my ($name) = keys $href;
        my $url = $href->{$name};

        my $id = $threads->tid();

        my $ua = LWP::UserAgent->new();
        
        #add the new cookie
        $ua->cookie_jar({ file => $cookie });
        
        my $response = $ua->get("genome.jgi.doe.gov$url");


        if ( $response->is_success ) {
            my $content = $response->decoded_content;
        }
        else {
            my $content = "URL: genome.jgi.doe.goc$url was not successfully downloaded.";
        }


        open my $OUT, ">", "$name.fastaq.gz";
        print $OUT $content;

    }

        #my $status = getstore("genome.jgi.doe.gov$url", "$name.fastq.gz");

        #if is_success($status) {
        #    print "File $name.fastq.gz successfully downloaded\n";
        #}
        #else {
        #    print "URL:  genome.jgi.doe.gov$url was not downloaded successfully.\n";
        
        #}
    #}
}




    
sub use_genbank {
    my ($genbank) = @_;

    open my $IN, "<", $genbank;

    my $genome_name;
    #read until the definition line
    while ( my $line = readline $IN ) {
        
        #returns 0 if this isn't the right line
        if ( $genome_name = extract_genome_name($line)) {
            last;
        }
    }

    if ( ! $genome_name ) {
        print "File X could not be used\n";
        return;
    }

    print "$genome_name\n";

}







sub extract_genome_name {
    my ($line) = @_;

    my $genome_name;

    #regex to look for the DEFINITION line and capture the genome name
    if ( $line =~ m{\A              #from the beginning of the line
                    DEFINITION\s\s  #the word DEFINITION and 2 whitespace characters
                    (.+)            #capture all the characters...
                    \s:             #until a space followed by a :
                    }xms ) {

        $genome_name = $1;
    }

    else {
        return 0;
    }



    ####this step will be impossible because JGI doesn't use consistent conventions
    #modify the genome name to match JGI XML
    #$genome_name =~ s{\s|\.}{}g;      #replace . and space with nothing

    return $genome_name;

}









     
