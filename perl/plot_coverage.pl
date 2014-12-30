#!/usr/bin/perl

####Graphs contigs with coverage that is at least a standard deviation above the norm to check for repetition or actual coverage


use strict;
use warnings;

use Cwd 'getcwd';
use Getopt::Long;
use Statistics::Basic qw(mean stddev);
use GD::Graph::lines;

use BioUtils::FastaIO;
use BioUtils::FastaSeq;

use Directory_obj;


use Data::Dumper;

my $USAGE = "plot_coverage.pl -file <my_file.txt>\n";

my $in_path;
my $out = getcwd();

GetOptions (    'in=s' =>       \$in_path,
                'out=s' =>      \$out,
                );


if ( !defined $in_path ) {
    die $USAGE;
}


my $directory = Directory_obj->new({'directory' => $in_path, 'ext' => ".txt"});


while ( my $file = $directory->next_file() ) {
    my $prefix = $out;
    my $filename = $directory->get_filename();
    if ( $directory->get_number_of_files() > 1 ) {
        -d "$prefix/$filename" or mkdir "$prefix/$filename";
        $prefix .= "/$filename";
    }    
    main($file, $filename, $prefix);
}




sub main {
    
    my ($file, $filename, $out) = @_;

    #get the coverage
    open my $IN, "<", "$file" or die "Could not open coverage file\n";


    my %data;
    my $data = \%data;

    #cycle through the file and read in all the depth notes
    while ( my $line = readline $IN ) {
        chomp $line; 
        my ($contig, $base, $cov) = split "\t", $line;

        #add all the bases to the appropriate contig
        $data->{$contig}{bases}[@{$data{$contig}{bases}}] = $cov;
    }

    
    my @mean_covs;

    #cycle through the contigs and calculate mean and stdev
    foreach my $contig ( sort keys %{$data} ) {
                
        my $mean = mean( $data->{$contig}{bases} );
        my $stddev = stddev( $data->{$contig}{bases} );
                
        $data->{$contig}{mean} = $mean;
        
        push @mean_covs, $mean;
            
    }
    

    #get the stddev for all the contigs in the file
    my $file_stddev = stddev( @mean_covs );
    my $file_mean = mean ( @mean_covs );
    #cycle through the contigs again to calculate which ones are a stdev away
    foreach my $contig ( sort keys %{$data} ) {

        #if at least 1 stddev away
        if ( $data->{$contig}{mean} > $file_mean + $file_stddev ) {
            
            #make a directory to store the results
            #-d "$out/$filename" or mkdir "$out/$filename";
            
            
            #build the dataset to graph
            my @plot_data = ( [], [] );
            my $y_max = 0;
            for ( my $i = 1; $i <= @{$data->{$contig}{bases}}; $i++ ) {
                my (@base, @coverage);
                
                push @{$plot_data[0]}, $i;
                push @{$plot_data[1]}, $data->{$contig}{bases}[$i-1];

                $y_max = $data->{$contig}{bases}[$i-1] if ( $y_max < $data->{$contig}{bases}[$i-1] );
                  
            }
            
            my $x_max = scalar @{$plot_data[0]};
            my $x_tick = $x_max / 5;
            my $graph = GD::Graph::lines->new(1000, 600);
            
            $graph->set(    transparent => '0',
                            x_label     => "Base Position",
                            y_label     => "Coverage",
                            title       => "Coverage plot for $contig",
                            x_max_value => $x_max,
                            y_max_value => $y_max,
                            bgclr       => 'white',
                            fgclr       => 'black',
                            x_tick_number => 'auto',
                            t_margin    => 5,
                            b_margin    => 5,
                            l_margin    => 5,
                            r_margin    => 40,
                            ) or die $graph->error();

            open my $OUT, ">", "$out/$contig.png" or die "Could not open $out/$filename/$contig.png\n";

            my $graphic = $graph->plot(\@plot_data) or die $graph->error();
            
            print $OUT $graphic->png();
            


        
            #last;
        }

    }
}
