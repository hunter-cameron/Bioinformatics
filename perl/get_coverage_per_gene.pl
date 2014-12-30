#!/usr/bin/perl

####Does bwa alignment on two folders; matches filenames

### Reccommended to run with multiple processors 


use strict;
use warnings;

use Cwd 'getcwd';
use Getopt::Long;
use Statistics::Basic qw(mean stddev);


use BioUtils::FastaIO;
use BioUtils::FastaSeq;

use Directory_obj;


use Data::Dumper;

my $USAGE = "get_coverage_per_gene.pl -gbk </gbk/> -depth </depth/> \n";

#number of threads to run
my $CPU = 10;

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

            my $path = "$prefix/$name\_temp/";

            my $genes = get_genes($gbk);

            my $output = "$prefix/$name\_genecov.txt";
            get_coverage($depth, $genes, $output);
        }
    }

    #last;
}

            

            



sub get_genes {
    my ($gbk) = @_;

    my %genes;

    use Bio::SeqIO;

    my $seqio_object = Bio::SeqIO->new( -file => $gbk );       
    while (my $seq_object = $seqio_object->next_seq) {

        #get the contig
        my $id = $seq_object->display_id();
        $genes{$id} = {};
        #print "$id\n";

        #get the tag objects
        for my $feat_object ($seq_object->get_SeqFeatures) {
        
            #only get info from the CDS tag
            if ($feat_object->primary_tag() eq "CDS" ) {
        
                my ($locus_tag) = $feat_object->get_tag_values('locus_tag') if ($feat_object->has_tag('locus_tag'));
                $genes{$id}->{$locus_tag}{start} = $feat_object->location()->start();
                $genes{$id}->{$locus_tag}{end} = $feat_object->location()->end();
        
         
                ($genes{$id}->{$locus_tag}{mol_type}) = $feat_object->get_tag_values("mol_type") if ($feat_object->has_tag("mol_type"));
                ($genes{$id}->{$locus_tag}{product}) = $feat_object->get_tag_values("product") if ($feat_object->has_tag("product"));

                

                #replace undefined tags with "null"
                $genes{$id}->{$locus_tag}{mol_type} = "null" if ( !defined $genes{$id}->{$locus_tag}{mol_type} );
                $genes{$id}->{$locus_tag}{product} = "null" if ( !defined $genes{$id}->{$locus_tag}{product} );

                #last;
            }
        }
    }

        return \%genes;
        #print Dumper %genes;
}




sub get_coverage {

    my ($file, $genes, $output) = @_;

    open my $IN, "<", $file or die "Could not open $file for reading\n";

    while ( my $line = readline $IN ) {
        chomp $line; 
        my ($contig, $base, $cov) = split "\t", $line;

        #add all the bases to the appropriate contig
        $genes->{$contig}{bases}[@{$genes->{$contig}{bases}}] = $cov;
    }

    close $IN; 
    open my $OUT, ">", $output or die "Could not write $output \n";
    print $OUT join "\t", ( "ID",
                            "Contig",
                            "Mean Coverage",
                            "Contig Mean Coverage",
                            "Genome Mean Coverage",
                            "Contig Coding Density",
                            "Genome Coding Density",
                            "Product",
                            "Molecule Type",
                            "Start",
                            "End"
                            ), "\n";

    #get the coverage per contig and genomic
    my @coverages;


    ($genes, my $genome_coverage, my $genome_density) = group_calculations($genes);
    



    my $contig_ID = 1;
    my $gene_ID = 1;
    foreach my $contig ( sort keys %{$genes} ) {


        #sort each contig in ascending order of genes
        foreach my $gene ( sort keys %{$genes->{$contig}} ) {
           

            

            #make a reference to avoid typing out this hash path again
            my $ref = $genes->{$contig}{$gene};

            #skip if the gene is not a hash ref (skips the mean and the bases array)
            next if ( ref $ref ne ref {} );

            #get the start and end; subtract 1 to turn into base 0 array notation
            my $start = $ref->{start} - 1;
            my $end = $ref->{end} - 1;

            #take an array slice with only the region of the gene
            my @gene_slice = @{$genes->{$contig}{bases}}[$start..$end];
                
            my $mean = mean( @gene_slice );
                
            $ref->{mean} = $mean;
                
            print $OUT join "\t", ( "C$contig_ID-G$gene_ID",
                                    $contig,
                                    $mean,
                                    $genes->{$contig}{mean},
                                    $genome_coverage,
                                    $genes->{$contig}{density},
                                    $genome_density,
                                    $ref->{product},
                                    $ref->{mol_type},
                                    $start,
                                    $end
                                    ), "\n";

            #print Dumper $data->{$contig};
            #last;

            $gene_ID++;
        }
         
        $contig_ID++;
        $gene_ID = 1;
    }

    close $OUT;

}


sub group_calculations {

    my ($genes) = @_;

    my @coverages;
    my @densities;

    foreach my $contig (keys %{$genes}) {
        my $mean = mean( $genes->{$contig}{bases} );
        $genes->{$contig}{mean} = $mean;
        

        #get the total number of coding bp
        my $coding_bp = 0;
        foreach my $gene ( keys %{$genes->{$contig}} ) {

            my $ref = $genes->{$contig}{$gene};
            next if ( ref $ref ne ref {} );
            
            my $length = $ref->{end} - $ref->{start};

            $coding_bp += $length;
        }

        #divide by the number of bases for the density
        my $density = $coding_bp / scalar @{$genes->{$contig}{bases}};

        $genes->{$contig}{density} = $density;
        
        push @densities, $density;
        push @coverages, $mean;
    }

    my $genome_coverage = mean ( @coverages );
    my $genome_density = mean ( @densities );

    return ( $genes, $genome_coverage, $genome_density );
}




