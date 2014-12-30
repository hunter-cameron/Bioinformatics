#!/usr/bin/perl

####Does bwa alignment on two folders; matches filenames

### Reccommended to run with multiple processors 


use strict;
use warnings;

use Cwd 'getcwd';
use Getopt::Long;
use Statistics::Basic qw(mean stddev);

use Directory_obj;
use PlasmidFinder::Genome;


use Data::Dumper;

my $USAGE = "finder.pl -gbk </gbk/> -depth </depth/> \n";

my $gbk_dir;
my $depth_dir;
my $prefix = getcwd();
my $append;



GetOptions (    'gbk=s' =>    \$gbk_dir,
                'depth=s' =>    \$depth_dir,
                'prefix=s' =>   \$prefix,
                'append' =>     \$append,
                );


if ( !defined $gbk_dir  or !defined $depth_dir ) {
    die $USAGE;
}

#print "$fasta_dir\n$fastq_dir\n";

my $GBK = Directory_obj->new({'directory' => $gbk_dir, 'ext' => ".gbk",});
my $DEPTH = Directory_obj->new({'directory' => $depth_dir, 'ext' => ".txt",});

while ( my $gbk = $GBK->next_file() ) {

    my $name = $GBK->get_filename();
    #print "$name\n";
    foreach my $depth ( $DEPTH->get_files() ) {

        if ( $depth =~ m/$name/ ) {

            #mkdir "$prefix/$name\_temp";
            
            my $genome = PlasmidFinder::Genome->new({ 'name' => $name, 'genbank' => $gbk, 'depth' => $depth });
           
            print "Genomic Coverage = ", $genome->get_coverage_depth(), "\n";
            print "Genomic Coding Density = ", $genome->get_coding_density(), "\n";


            foreach my $contig ( $genome->get_contig_refs() ) {

                print join "\t", (  $contig->get_name(), 
                                    $contig->get_length(), 
                                    $contig->get_coverage_depth(), 
                                    $contig->get_coding_density(),
                                    "\n",
                                    );
               
                foreach my $gene ( $contig->get_gene_refs() ) {

                    print join "\t", (  "\t",
                                        $gene->get_name(),
                                        $gene->get_coverage_depth(),
                                        $gene->get_start(),
                                        $gene->get_end(),
                                        $gene->get_product(),
                                        "\n"
                                        );

                    last;
                }
            }
        }
    }

    last;
}

            

            




