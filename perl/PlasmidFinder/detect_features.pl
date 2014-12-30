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

my $USAGE = "finder.pl -gbk </gbk/> -depth </depth/> \n";

#deviation to consider significant
my $DEV = 1;





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
            my $gene_coverage = calc_gene_coverage($genome);
            
            open my $OUT, ">>", "$prefix/$name\_out.txt";

            foreach my $contig ( sort $genome->get_contigs() ) {

                print $OUT $contig->get_name(), "\n";

                foreach ( sort keys %{$gene_coverage->{contigs}{$contig}} ) {

                    print $OUT "\t$_ = \t", $gene_coverage->{contigs}{$contig}{$_}, "\n";
                }

                if ( $gene_coverage->{contigs}{$contig}{frac_high} >= .90 ) {

                print $OUT "\n";

                foreach my $gene ( sort $contig->get_genes() ) {

                    print $OUT "\t", $gene->get_name(), "\n" if ( defined $gene_coverage->{genes}{$gene} );

                    foreach (sort keys %{$gene_coverage->{genes}{$gene}} ) {

                        print $OUT "\t\t$_ = \t", $gene_coverage->{genes}{$gene}{$_}, "\n";
                    }

                }
                }

                print $OUT "\n", "-" x 40, "\n";
            }
            
            print $OUT "x" x 40, "\nEND OF GENOME $name\n", "x" x 40, "\n";
            close $OUT;

        }
    }

    #last;
}


sub calc_gene_coverage {

    my ($genome) = @_;



    #calculate mean coverage and stdev
    my @coverages;
    foreach my $contig ( $genome->get_contigs() ) {

        foreach my $gene ( $contig->get_genes() ) {

            push @coverages, $gene->get_coverage_depth();
        }
    }

    my $mean = mean @coverages;
    my $stddev = stddev @coverages;


    #calculate genes that deviate by the specified value
    my $data = {};
    foreach my $contig ( $genome->get_contigs() ) {

        #how would I call scalar on the get_genes() to avoid the temp array?
        my @temp = $contig->get_genes();
        $data->{contigs}{$contig}{num_genes} = scalar @temp;
        $data->{contigs}{$contig}{high_cov} = 0; 
        foreach my $gene ( $contig->get_genes() ) {

            my $coverage = $gene->get_coverage_depth;

            my $dev = ($coverage - $mean) / $stddev;

            if ( $dev >= $DEV ) {

                $data->{genes}{$gene}{depth} = $dev;

                $data->{contigs}{$contig}{high_cov}++;
            }

        }
        
        #calculate the fraction of high coverage genes per contig
        if ( $data->{contigs}{$contig}{num_genes} != 0 and defined $data->{contigs}{$contig}{num_genes} ) {
            $data->{contigs}{$contig}{frac_high} = $data->{contigs}{$contig}{high_cov} / $data->{contigs}{$contig}{num_genes};
        }

        else {
            $data->{contigs}{$contig}{frac_high} = 0;
            $data->{contigs}{$contig}{num_genes} = 0;
        }
        

    }
    
    return $data;
}


    

            









