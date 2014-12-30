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

my $contigs = {};
while ( my $gbk = $GBK->next_file() ) {

    my $name = $GBK->get_filename();
    #print "$name\n";
    foreach my $depth ( $DEPTH->get_files() ) {

        if ( $depth =~ m/$name/ ) {

            #mkdir "$prefix/$name\_temp";
            
            my $genome = PlasmidFinder::Genome->new({ 'name' => $name, 'genbank' => $gbk, 'depth' => $depth });
            print "Genome making done.\n";
            print "Total genome: " . total_size($genome). "\n";
            calc_stats($genome, $contigs, $name);

            
        }
    }
    last;
}

print_results($contigs, "$prefix/contigs_data.txt");





sub calc_stats {

    my ($genome, $contigs, $genome_name) = @_;

    my $median = $genome->get_median_coverage_depth();
    my $coding_density = $genome->get_coding_density();
    my ($mean_cov, $stddev) = $genome->get_mean_contig_depth();

    foreach my $contig ( $genome->get_contigs() ) {

        $contigs->{$contig}{name} = $genome_name . " | " . $contig->get_name();
        $contigs->{$contig}{length} = $contig->get_length();
        $contigs->{$contig}{coding_density} = $contig->get_coding_density();
        $contigs->{$contig}{coding_diff} = $contigs->{$contig}{coding_density} - $coding_density;
        $contigs->{$contig}{high_cov} = 0;
        $contigs->{$contigs}{num_genes} = 0;

        $contigs->{$contig}{coverage} = $contig->get_coverage_depth();
        $contigs->{$contig}{stddevs} = ($contigs->{$contig}{coverage} - $mean_cov) / $stddev;
        

        foreach my $gene ( $contig->get_genes() ) {

            $contigs->{$contig}{high_cov}++ if ($gene->get_coverage_depth() / $median >= 1.8);
            
            $contigs->{$contig}{num_genes}++;            

           # last;
        }
    }



}




sub print_results {

    my ($contigs, $out_file) = @_;

    print "Total results: " . total_size($contigs). "\n";

    open my $OUT, ">", $out_file or die "Could not open $out_file\n";

    print $OUT join "\t" , (    "Name",
                                "Length",
                                "Coding Density",
                                "Difference from Mean Coding",
                                "Coverage",
                                "Coverage Stddevs",
                                "Number of Genes",
                                "Number of High Coverage Genes",
                                ), "\n";
    
            
    foreach my $contig ( keys %{$contigs} ) {


      
        print $OUT join "\t", (     $contigs->{$contig}{name},
                                    $contigs->{$contig}{length},
                                    $contigs->{$contig}{coding_density},
                                    $contigs->{$contig}{coding_diff},
                                    $contigs->{$contig}{coverage},
                                    $contigs->{$contig}{stddevs},
                                    $contigs->{$contig}{num_genes},
                                    $contigs->{$contig}{high_cov},
                                    ), "\n";
    }

    close $OUT;
}











