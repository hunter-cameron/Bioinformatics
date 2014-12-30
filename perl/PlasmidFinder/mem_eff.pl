#!/usr/bin/perl

####Looks for features within a genome



use strict;
use warnings;

use Cwd 'getcwd';
use Getopt::Long;
use Statistics::Basic qw(mean stddev);

use Directory_obj;
use PlasmidFinder::Genome;
use PlasmidFinder::Contig;
use PlasmidFinder::Gene;


use Data::Dumper;
use Devel::Size qw(size total_size);


my $USAGE = "dump_contigs.pl -gbk </gbk/> -depth </depth/> \n";


my $gbk_dir;
my $depth_dir;
my $prefix = getcwd();



GetOptions (    'gbk=s' =>    \$gbk_dir,
                'depth=s' =>    \$depth_dir,
                'prefix=s' =>   \$prefix,
                );


if ( !defined $gbk_dir  or !defined $depth_dir ) {
    die $USAGE;
}

#print "$fasta_dir\n$fastq_dir\n";

my $GBK = Directory_obj->new({'directory' => $gbk_dir, 'ext' => ".gbk",});
my $DEPTH = Directory_obj->new({'directory' => $depth_dir, 'ext' => ".txt",});


#open the output file and write header lines
open my $OUT, ">", "$prefix/contigs_data.txt" or die "Could not open $prefix/contigs_data.txt\n";

print $OUT join "\t" , (    "Name",
                            "Length",
                            "Coding Density",
                            "Difference from Mean Coding",
                            "Coverage",
                            "Coverage Stddevs",
                            "Number of Genes",
                            "Number of High Coverage Genes",
                            ), "\n";



while ( my $gbk = $GBK->next_file() ) {

    my $name = $GBK->get_filename();
    #print "$name\n";
    foreach my $depth ( $DEPTH->get_files() ) {

        if ( $depth =~ m/$name/ ) {

            #mkdir "$prefix/$name\_temp";
            
            my $genome = PlasmidFinder::Genome->new({ 'name' => $name, 'genbank' => $gbk, 'depth' => $depth });
            calc_stats($genome, $name);

            #forced recycling
            #$genome->DESTROY(); 
        }
    }
    #last;
}

close $OUT;




sub calc_stats {

    my ($genome, $genome_name) = @_;

    my $median = $genome->get_median_coverage_depth();
    my $mean_coding_density = $genome->get_coding_density();
    my ($mean_cov, $stddev) = $genome->get_mean_contig_depth();

    foreach my $contig ( $genome->get_contigs() ) {

        my $name = $genome_name . " | " . $contig->get_name();
        my $length = $contig->get_length();
        my $coding_density = $contig->get_coding_density();
        my $coding_diff = $coding_density - $mean_coding_density;
        my $high_cov = 0;
        my $num_genes = 0;

        my $coverage = $contig->get_coverage_depth();
        my $stddevs = ($coverage - $mean_cov) / $stddev;
         

        foreach my $gene ( $contig->get_genes() ) {

            $high_cov++ if ($gene->get_coverage_depth() / $median >= 1.8);
            
            $num_genes++;            

           # last;
        }

        print $OUT join "\t", ( $name,
                                $length,
                                $coding_density,
                                $coding_diff,
                                $coverage,
                                $stddevs,
                                $num_genes,
                                $high_cov
                                ), "\n";

    }



}




