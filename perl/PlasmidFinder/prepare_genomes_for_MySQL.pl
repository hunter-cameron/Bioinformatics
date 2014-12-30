#!/usr/bin/perl

####Processes genomes for MySQL



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


my $USAGE = "prepare_genomes_for_MySQL.pl -gbk </gbk/> -depth </depth/> \n";


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


my $genomes = {};
my $contigs = {};
my $genes = {};

#write headers
open my $OUT, ">", "$prefix/genomes_stats.txt";
print $OUT "Name\tSize\tMed_Cov\tAvg_Cov\tCoding_Density\n";
close $OUT;

open $OUT, ">", "$prefix/contigs_stats.txt";
print $OUT "Name\tSize\tMed_Cov\tAvg_Cov\tCoding_Density\tRef_Genome\n";
close $OUT;

open $OUT, ">", "$prefix/genes_stats.txt";
print $OUT "Name\tSize\tStat_Pos\tEnd_Pos\tMed_Cov\tAvg_Cov\tProduct\tRef_Genome\tRef_Contig\n";
close $OUT;


while ( my $gbk = $GBK->next_file() ) {

    my $name = $GBK->get_filename();
    #print "$name\n";
    foreach my $depth ( $DEPTH->get_files() ) {

        if ( $depth =~ m/$name/ ) {

            #mkdir "$prefix/$name\_temp";
            
            my $genome = PlasmidFinder::Genome->new({ 'name' => $name, 'genbank' => $gbk, 'depth' => $depth });
            print "Genome making done.\n";
            print "Total genome: " . total_size($genome). "\n";

            ($genomes, $contigs, $genes) = calc_genome($genome, $genomes);
            append_results($genomes, $contigs, $genes, $prefix);
        }
    }
    #last;
}



sub calc_genome {

    my ($genome, $genomes, $contigs, $genes) = @_;
    
    $genomes->{name} = $genome->get_name();
    $genomes->{length} = $genome->get_length();
    $genomes->{med_cov} = $genome->get_median_coverage_depth();
    $genomes->{avg_cov} = $genome->get_coverage_depth();
    $genomes->{coding_density} = $genome->get_coding_density();
    $genomes->{num_genes} = 0;


    #set undefined values to MySQL null
    foreach (keys %{$genomes}) {
        $genomes->{$_} = '\N' if ( ! defined $genomes->{$_} );
    }

    foreach my $contig ( $genome->get_contigs() ) {
        
        $genomes->{num_contigs}++;

        $contigs->{$contig}{name} = $contig->get_name();
        $contigs->{$contig}{length} = $contig->get_length();
        $contigs->{$contig}{num_genes} = 0;
        $contigs->{$contig}{med_cov} = $contig->get_median_coverage_depth();
        $contigs->{$contig}{avg_cov} = $contig->get_coverage_depth();
        $contigs->{$contig}{coding_density} = $contig->get_coding_density();
        $contigs->{$contig}{ref_genome} = $genomes->{name};

        #set undefined values to MySQL null
        foreach (keys %{$contigs->{$contig}}) {
            $contigs->{$contig}{$_} = '\N' if ( ! defined $contigs->{$contig}{$_} );
        }
    

        foreach my $gene ( $contig->get_genes() ) {
            $genomes->{num_genes}++;
            $contigs->{$contig}{num_genes}++;

            $genes->{$gene}{name} = $gene->get_name();
            $genes->{$gene}{length} = $gene->get_length();
            $genes->{$gene}{start_pos} = $gene->get_start();
            $genes->{$gene}{end_pos} = $gene->get_end();
            
            ($genes->{$gene}{avg_cov}, $genes->{$gene}{med_cov}) = $gene->get_coverage_depth();
            $genes->{$gene}{product} = $gene->get_product();
            #replace any , in the product with _
            $genes->{$gene}{product} =~ tr/,/_/;

            $genes->{$gene}{ref_genome} = $genomes->{name};
            $genes->{$gene}{ref_contig} = $contigs->{$contig}{name};

            #set undefined values to MySQL null
            foreach (keys %{$genes->{$gene}}) {
                $genes->{$gene}{$_} = '\N' if ( ! defined $genes->{$gene}{$_} );
            }

        }

    }
    return($genomes, $contigs, $genes);
}


sub append_results {

    my ($genomes, $contigs, $genes, $prefix) = @_;

    open my $G_OUT, ">>", "$prefix/genomes_stats.txt" or die "Could not write\n";

    print $G_OUT join("\t", ( $genomes->{name},
                            $genomes->{length},
                            #$genomes->{num_contigs},
                            #$genomes->{num_genes},
                            $genomes->{med_cov},
                            $genomes->{avg_cov},
                            $genomes->{coding_density},
                            )), "\n";


    close $G_OUT;



    open my $C_OUT, ">>", "$prefix/contigs_stats.txt" or die "Could not write\n";

    foreach my $contig ( keys %{$contigs} ) {
        print $C_OUT join("\t", ( $contigs->{$contig}{name},
                                $contigs->{$contig}{length},
                                #$contigs->{$contig}{num_genes},
                                scalar $contigs->{$contig}{med_cov},
                                $contigs->{$contig}{avg_cov},
                                $contigs->{$contig}{coding_density},
                                $contigs->{$contig}{ref_genome}
                                )), "\n";
    }

    close $C_OUT;



    open my $Ge_OUT, ">>", "$prefix/genes_stats.txt" or die "Could not write\n";
    
    foreach my $gene ( keys %{$genes} ) {
        print $Ge_OUT join("\t", ( $genes->{$gene}{name},
                                $genes->{$gene}{length},
                                $genes->{$gene}{start_pos},
                                $genes->{$gene}{end_pos},
                                $genes->{$gene}{med_cov},
                                $genes->{$gene}{avg_cov},
                                $genes->{$gene}{product},
                                $genes->{$gene}{ref_genome},
                                $genes->{$gene}{ref_contig}
                                )), "\n";
    }

    close $Ge_OUT;


}





