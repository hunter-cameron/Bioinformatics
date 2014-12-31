
use strict;
use warnings;
use Bio::SeqIO;

use Data::Dumper;

scalar @ARGV == 2 or die("Usage: name.pl <sort_list_aln> <gbk_dir>\n");

my ($list, $dir) = @ARGV;
my $list_w_gene = $list . "plus_gene";

open my $LIST, "<", $list;
open my $NEW, ">", $list_w_gene;




my $previous = {"contig" => ''};
my $gbk_file = '';
my $bio_gbk;
my $indx = 0;
my @genes;      # use a real copy to avoid making a local copy for each alignment
my $gene_data;

# reporting statistics
my $total_alns = 0;
my $not_found = 0;


print $NEW "Sample\tGenome\tContig\tRead\tMapQ\tLocus tag\tGene start\tGene end\tGene product\n";

while ( my $aln = get_next_alignment($LIST) ) {
    $total_alns++;
    
    # continue to use same set of genes if contig is the same
    #print $aln->{contig}, "\t", $previous->{contig}, "\n";
    
    if ( $aln->{contig} ne $previous->{contig} ) {
        
        # check if the contig is in the same gbk
        my $file = get_gbk_from_contig($aln->{contig}, $dir);
        print "$gbk_file\t$file\n";
        if ( $file ne $gbk_file ) {       # open a new one if not
            print "Opening $file\n";
            $gbk_file = $file;
            $bio_gbk = Bio::SeqIO->new( -file => $gbk_file );
        }

        # get genes ref and unpack it
        print "Getting genes from $aln->{contig}\n";
        (my $genes, $gene_data) = get_contig_genes_from_gbk($bio_gbk, $aln->{contig});
        
        # if the contig is not found in order, open again and check from beginning
        # die if not found again
        if ($genes == 0) {
            warn "List or gbk probably poorly sorted...searching again.\n";
            $bio_gbk = Bio::SeqIO->new( -file=> $gbk_file );
            ($genes, $gene_data) = get_contig_genes_from_gbk($bio_gbk, $aln->{contig});
            die "Cound not find GenBank info for contig $aln->{contig} in $bio_gbk\n" if $genes == 0;
        }

        @genes = @{$genes};
        $indx = 0;
    }

    $previous = $aln;

    #
    ## Begin processing genes
    #
    
    # need to ensure each alignment has a fair chance to be categorized: must be greater than at least one of the genes
    #

    my $flag = 0;
    while ( $indx < scalar @genes ) {
   
        my $gene = gene2hash($genes[$indx], $gene_data);
        #print Dumper $gene;

        if ( $aln->{start} > $gene->{end_pos} ) {
            # alignment is after this gene
            # keep moving- possibly remove gene
            $flag = 1;
            $indx++;
            if ( $indx == scalar @genes ) {     # test is the current gene is the last gene in the contig
                # didn't find in that contig
                warn "Aln $aln->{read} at position $aln->{start} was not included in a gene.\n";
                $not_found++;
                $gene->{start_pos} = $aln->{start};
                $gene->{end_pos} = "None";
                $gene->{product} = "None";
                $gene->{name} = "Nongene";
                print_new($aln, $gene);
                last;
            }

        }

        elsif ( $aln->{start} < $gene->{start_pos} ) {
            # gene is past aln, mapping place was skipped    
            if ( $flag ) {      # aln given fair chance - just didn't work
                warn "Aln $aln->{read} at position $aln->{start} was not included in a gene.\n";
                $not_found++;
                $gene->{start_pos} = $aln->{start};
                $gene->{end_pos} = "None";
                $gene->{product} = "None";
                $gene->{name} = "Nongene";
                print_new($aln, $gene);
                last; 
            }
            else {      # aln not given fair chance repeat going back one
                # go down an index
                if ( $indx > 0 ) {
                    $indx -= 1;
                }
                else { 
                    warn "Aln $aln->{read} at position $aln->{start} was not included in a gene.\n";
                    $not_found++;
                    $gene->{start_pos} = $aln->{start};
                    $gene->{end_pos} = "None";
                    $gene->{product} = "None";
                    $gene->{name} = "Nongene";
                    print_new($aln, $gene);
                    last; 
                }
            }
        }

        elsif ( $aln->{start} >= $gene->{start_pos} and $aln->{start} <= $gene->{end_pos} ) {
            # alignment is to this gene !!! need something to skip the warning message
            print_new($aln, $gene);
            last;
        }

        else {  # should never get here
            die "An error occurred in the loop that detects matches. This line should have never been reached.\n";
        }

    }
    
    # last;
}

print "$not_found of $total_alns total alignments were not found in the CDS, tRNA, or rRNA of the Genbanks.\n";

sub get_gbk_from_contig {
    # returns a filepath from a contig and the directory of gbk
    # will mess up if the contig name is found in two files
    
    my ( $contig, $dir ) = @_;

    my $result = `grep -m 1 "$contig" $dir/*`;

    my $file = '';
    if ( $result =~ m/^$dir\/(\d+\.gbk):/ ) {
        $file = $1;
        $file = "$dir/$file";
        print "Found in gbk = $file\n" ;
    }
    else {
        die "Contig $contig not found!\n";
    }

    return $file;
}

sub get_contig_genes_from_gbk {
    # parses gbk files for a given contig ( given in order they appear in gbk )

    my ($gbk, $contig) = @_;
 
    while ( my $seq_object = $gbk->next_seq() ) {
        my %gene_data = ();
        #get the contig
        $gene_data{contig} = $seq_object->display_id();
        print "$gene_data{contig}\n";
        print $contig, "\n";
        if ( $gene_data{contig} =~ m/$contig/ ) {       # the matching is somewhat of a hack b/c I forgot formats work as scaffold_0001.1 instead of scaffold_0001
            
            #remove the contig data from the DEFINITION line
            my $def = ($seq_object->desc());
            #print "$description\n";
            ($gene_data{genome}) = split ":", $def;

            my @genes = ();

            for my $candidate ( $seq_object->get_SeqFeatures ) {
                if ( $candidate->primary_tag() eq "CDS" or $candidate->primary_tag() eq "rRNA" or $candidate->primary_tag() eq "tRNA" ) {
                    push @genes, $candidate;
                }
            }
        return ( \@genes, \%gene_data );
        }

    }
    # if we get here, something has gone wrong (probably LIST is out of order)
    # return 0 and let parent figure out what to do
    return 0;    
}

sub get_next_alignment {
    # parses the next line of the sorted file and returns an alignment "object"
    my ($LIST) = @_;
    
    my $line = readline $LIST;
    return 0 if ! defined $line;    # check for end of file

    chomp $line;
    my %aln; 

    ($aln{sample}, $aln{read}, $aln{contig}, $aln{start}, $aln{map_q}) = split "\t", $line;

    return \%aln;
}

sub gene2hash {
    # converts a Bio Feature object (or something like that) to a hash with the params I care about
   
    my ( $gene, $gene_data ) = @_;

    ($gene_data->{name}) = $gene->get_tag_values('locus_tag') if ($gene->has_tag('locus_tag' ));
    $gene_data->{start_pos} = $gene->location()->start();
    $gene_data->{end_pos} = $gene->location()->end();

    #set optional parameter
    ($gene_data->{product}) = $gene->get_tag_values("product")   if ($gene->has_tag("product"));

    return $gene_data;
}

sub print_new {
    my ( $aln, $gene ) = @_;

    #print Dumper $aln, $gene;
    print $NEW join("\t", ( $aln->{sample},
                            $gene->{genome},
                            $gene->{contig},        # to get the full contig name
                            $aln->{read},
                            $aln->{map_q},
                            $gene->{name},
                            $gene->{start_pos},
                            $gene->{end_pos},
                            $gene->{product}
                        )), "\n";
}
