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

my $USAGE = "dump_gene_pool.pl -gbk </gbk/> -depth </depth/> \n";


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

my $genes = {};
while ( my $gbk = $GBK->next_file() ) {

    my $name = $GBK->get_filename();
    #print "$name\n";
    foreach my $depth ( $DEPTH->get_files() ) {

        if ( $depth =~ m/$name/ ) {

            #mkdir "$prefix/$name\_temp";
            
            my $genome = PlasmidFinder::Genome->new({ 'name' => $name, 'genbank' => $gbk, 'depth' => $depth });
            print "Genome making done.\n"; 
            calc_stats($genome, $genes, $name);

            
        }
    }
    last;
}

print_results($genes, "$prefix/genes_data.txt");





sub calc_stats {

    my ($genome, $genes, $genome_name) = @_;

    my $median = $genome->get_median_coverage_depth();

    foreach my $contig ( $genome->get_contigs() ) {

        foreach my $gene ( $contig->get_genes() ) {

            $genes->{$gene}{name} = $genome_name . " | " . $gene->get_name();
            $genes->{$gene}{location} = join "", (     $gene->get_contig()->get_name(), 
                                                    " => ",
                                                    $gene->get_start(),
                                                    "..",
                                                    $gene->get_end()
                                                    );
            $genes->{$gene}{coverage} = $gene->get_coverage_depth();
            
            $genes->{$gene}{copies} = $genes->{$gene}{coverage} / $median;
            $genes->{$gene}{product} = $gene->get_product();
           # last;
        }
    }

}




sub print_results {

    my ($genes, $out_file) = @_;

    open my $OUT, ">", $out_file or die "Could not open $out_file\n";

    print $OUT join "\t" , (    "Name",
                                "Location", #contig name => start..end
                                "Coverage",
                                "Est. Copies",
                                "Product"
                                ), "\n";
    
            
    foreach my $gene ( keys %{$genes} ) {


      
        print $OUT join "\t", (    $genes->{$gene}{name},
                                    $genes->{$gene}{location},
                                    $genes->{$gene}{coverage},
                                    $genes->{$gene}{copies},
                                    $genes->{$gene}{product}
                                    ), "\n";
    }

    close $OUT;
}











