#!/nas02/apps/perl-5.18.2/bin/perl


####NEEDS TO CHECK FOR FASTAS W/O MATCHING FSQ AND NEEDS BETTER REPORTING ON UNDEFINED VALUES




#Does bwa alignment on two folders; matches filenames

### Reccommended to run with multiple processors 
### run from the same directory to skip already processed files. 

use strict;
use warnings;

use Cwd 'getcwd';
use Getopt::Long;
use Statistics::Basic qw(mean stddev);


use BioUtils::FastaIO;
use BioUtils::FastaSeq;

use Directory_obj;


use Data::Dumper;

my $USAGE = "get_coverage_by_directory.pl -fasta </fasta/> -fastq </fastq/> [-prefix mypath] -gbk </gbk/> \n";

#number of threads to run
my $CPU = 10;

my $fasta_dir;
my $fastq_dir;
my $gbk_dir;
my $PAIRED_END;
my $prefix = getcwd();



GetOptions (    'fasta=s' =>    \$fasta_dir,
                'fastq=s' =>    \$fastq_dir,
                'gbk=s'   =>    \$gbk_dir,
                'pe|PE' =>      \$PAIRED_END,
                'prefix=s' =>   \$prefix,
                );


if ( !defined $fasta_dir  or !defined $fastq_dir ) {
    die $USAGE;
}

#print "$fasta_dir\n$fastq_dir\n";

#open directories corresponding to the 3 files needed
my $FSA = Directory_obj->new({'directory' => $fasta_dir, 'ext' => ".fasta",});
my $FSQ = Directory_obj->new({'directory' => $fastq_dir, 'ext' => ".fastq",});
my $GBK = Directory_obj->new({'directory' => $gbk_dir,   'ext' => ".gbk",});

while ( my $fasta = $FSA->next_file() ) {

    my $name = $FSA->get_filename();

    #get corresponding gbk file
    my $gbk_file;
    foreach my $gbk ( $GBK->get_files() ) {
       print "$name \t $gbk\n"; 
        if ( $gbk =~ m/\/$name\.gbk/ ) {
            $gbk_file = $gbk;
            #last;
        }
    }
    if ( !defined $gbk_file ) {
        warn "No Genbank found for fasta: $fasta. Skipping file.\n";
        next;
    }
    

    #get corresponding fastq file
    my @fastq;
    foreach my $fastq ( $FSQ->get_files() ) {

        if ( $fastq =~ m/\/$name\_?R?\d?\.fastq/ ) {  #match /,  the name, and optional _Rx for read number
            push @fastq, $fastq;
        }

    }

    #if not the right number of fastq files, skip this one
    if ( ($PAIRED_END and @fastq != 2) or (!$PAIRED_END and @fastq != 1) ) {

        warn "Incorrect number of fastq files (", scalar @fastq, ") for fasta: $fasta. Skipping file.\n";
        next;
    }

            
    print "Current genome: $name\n";
    my $path = "$prefix/$name\_data/";
    my $depth = process_fastq(\@fastq, $fasta, $path, $name);

    if (-e "$path/$name\_genome_stats.txt" and -e "$path/$name\_contig_stats.txt" and -e "$path/$name\_gene_stats.txt") {
        print "Genome $name already completed\n";
    }

    else {
        my ($contigs, $genome) = calc_stats($gbk_file, $fasta, $depth, $name);

        #print Dumper $contigs;
        print_results($path, $name, $genome, $contigs);
    }

}



sub print_results {
    my ($path, $name, $genome, $contigs) = @_;

    open my $G_OUT, ">", "$path/$name\_genome_stats.txt";
    #print headers
    print $G_OUT "name\tdescription\tsize\tmed_cov\tavg_cov\tcode_density\n";

    print $G_OUT join( "\t", (  $name,
                                $genome->{description},
                                $genome->{length},
                                $genome->{med_cov},
                                $genome->{avg_cov},
                                $genome->{code_density}
                                )), "\n";

    close $G_OUT;

    
    open my $C_OUT, ">", "$path/$name\_contig_stats.txt";
    
    #print headers
    print $C_OUT "name\tsize\tmed_cov\tavg_cov\tcode_density\tgc_content\tref_Genome\n";
    
    open my $Ge_OUT, ">", "$path/$name\_gene_stats.txt";

    #print headers
    print $Ge_OUT "name\tsize\tstart_pos\tend_pos\tmed_cov\tavg_cov\tproduct\tref_Genome\tref_Contig\n";

    foreach my $contig ( keys %{$contigs} ) {

        print $C_OUT join ( "\t", ( $contig,
                                    $contigs->{$contig}{length},
                                    $contigs->{$contig}{med_cov},
                                    $contigs->{$contig}{avg_cov},
                                    $contigs->{$contig}{code_density},
                                    $contigs->{$contig}{gc_content},
                                    $name
                                    )), "\n";

        foreach my $gene ( keys %{$contigs->{$contig}{genes}} ) {

            my $ref = $contigs->{$contig}{genes}{$gene};
            print $Ge_OUT join ( "\t", (    $gene,
                                            $ref->{length},
                                            $ref->{start_pos},
                                            $ref->{end_pos},
                                            $ref->{med_cov},
                                            $ref->{avg_cov},
                                            $ref->{product},
                                            $name,
                                            $contig
                                            )), "\n";
        }
    }
    close $C_OUT;
    close $Ge_OUT;
}




sub calc_stats {
    my ($gbk, $fasta, $depth, $name) = @_;

    my $data = {};
    my $genome = {};
 
 
    #foreach ( keys %{$data} ) {print "Calc_stats = $_\n";}
    #parse the fasta for GC content, size
    parse_fasta($data, $fasta);
    
    #get genes and set up the contig names for the fasta (in case fasta has appended names; eg length & cov
    parse_genbank($data, $genome, $gbk);


    
    #parse depth file and calculate coverage
    parse_depth($data, $genome, $depth);

       
    #calculate the genomic and contig coding density
    calculate_density($data, $genome);
    

    #check data, mostly for depth (if no reads mapped back, values are null)
    foreach my $contig (keys %{$data}) {
        if ( !defined $data->{$contig}{med_cov} ) {
            $data->{$contig}{med_cov} = 0;
            $data->{$contig}{avg_cov} = 0;
        }
    }

    return $data, $genome;
}




sub calculate_density {

    my ($data, $genome) = @_;

    #print Dumper $data;

    my $genomic_coding = 0;

    #foreach (keys %{$data} ) {
    #    print "Sample: $_\n";
    #}

    foreach my $contig ( keys %{$data} ) {
     
        my $contig_coding = 0;

        #add bases from each gene
        foreach my $gene ( keys %{$data->{$contig}{genes}} ) {
            
            $contig_coding += $data->{$contig}{genes}{$gene}{length};
        }
        
        #print "Contig = $contig\n";
        $data->{$contig}{code_density} = $contig_coding / $data->{$contig}{length};
        #add coding bases to genome
        $genomic_coding += $contig_coding;
    }

    #calculate genomic density
    $genome->{code_density} = $genomic_coding / $genome->{length};
}

   

#review algorithm for getting base coverage and calculating; since I have to do genomic coverage, there may be no advantage to reading only a contig worth of bases at a time (rather than a hash for the whole genome)
sub parse_depth {
    my ($data, $genome, $depth_file) = @_;


    #foreach ( keys %{$data} ) {print "Depth = $_\n";}
    open my $IN, "<", $depth_file or die "Could not open $depth_file\n";

    my $current_con;
    my @bases;          #array to hold contig bases
    my @genome;         #array to hold genomic bases

    my $contig;
    while ( my $line = readline $IN ) {
            
        chomp $line;
            
        ($contig, my $base_pos, my $coverage) = split "\t", $line;

        #set the initial value of current_con
        $current_con = $contig if ( !defined $current_con );

            
        #if the contig doesn't match then end the current contig and send the array of bases and reset
        if ( ($contig ne $current_con) and @bases ) {
                
            #$depth{$current_con} = \@bases;

            #add in the uncovered bases - this should be unnecessary for bedtools
            #print scalar @bases, "\t", $data->{$current_con}{length}, "\n";
            
            #for (my $length = scalar @bases; $length <= $data->{$current_con}{length}; $length++ ) {
            #    push @bases, 0;
            #    warn "Bases not same length for $depth_file\n";
            #}
            my $contig_alias = _contig_match($current_con, $data);
            calc_cov_stats($data, \@bases, $contig_alias);
            
            #reset bases
            @bases = ();

            #set the new current contig
            $current_con = $contig;
         }
            
           
        push @bases, $coverage;
        push @genome, $coverage;
    }

    #add the bases to the last contig after the file ends
    #$depth{$contig} = \@bases;

    my $contig_alias = _contig_match($contig, $data);
    calc_cov_stats($data, \@bases, $contig_alias);
     
    ($genome->{avg_cov}, $genome->{med_cov}) = mean_median(\@genome);
    $genome->{length} = scalar @genome; 
    close $IN;

    #foreach ( keys %{$data} ) {print "Depth-after = $_\n";}
}


sub calc_cov_stats {
    my ($data, $bases, $contig) = @_;
    
    die "Duplicate of contig $contig. Aborting.\n" if (defined $data->{$contig}{avg_cov});

    #calc contig cov
    ($data->{$contig}{avg_cov}, $data->{$contig}{med_cov}) = mean_median($bases);

    #calc cov for each gene on the contig
    foreach my $gene ( keys %{$data->{$contig}{genes}} ) {

        my $ref = $data->{$contig}{genes}{$gene};

        my $start = $ref->{start_pos};
        my $end = $ref->{end_pos};

        #get an array slice (-1 for base 0)
        my @slice = @{$bases}[$start - 1..$end - 1];

        ($ref->{avg_cov}, $ref->{med_cov}) = mean_median(\@slice);

    }
}



sub mean_median {
    my ($data) = @_;
    my @data = @{$data};

    #calculate median
    @data = sort {$a <=> $b} @data;
    
    my $median;
    my $mid = int @data / 2;
    if ( @data % 2 ) {
        #if data set is odd, use middle
        $median = $data[$mid];
    }
    else {
        #if even, average 2 middle; use -1 because scalar @ is base 1 and @ index is base 0 
        $median = ($data[$mid - 1] + $data[$mid]) / 2;
    }

    #calculate mean
    my $sum = 0;
    foreach ( @data ) {
        $sum += $_;
    }

    my $mean = $sum / scalar @data;

    return $mean, $median;
}








#get the length (could be done with the depth file) and GC count
sub parse_fasta {

    my ($data, $fasta) = @_;

    
    #foreach ( keys %{$data} ) {print "Fasta = $_\n";}
    my $in = BioUtils::FastaIO->new({stream_type => '<', file => $fasta});

    while ( my $seq_obj = $in->get_next_seq() ) {
        my $header = $seq_obj->get_header();
        my $seq = $seq_obj->get_seq();

        #get the length
        my $length = length $seq;

        #figure out which contig to store it as
        #my $contig = _contig_match($header, $data);
        #die "Duplicated match to $contig for $header.\n" if (defined $data->{$contig}{length});
       
        my $contig = $header;
        $data->{$contig}{length} = $length;

        #get GC content
        my $GC = ($seq =~ tr/GgCc//);   #transliterate with nothing = count

        $data->{$contig}{gc_content} = $GC / $length;
    }

}

#finds the best match for a given contig based on 3mers
sub _contig_match {
    my ($header, $data) = @_;

    #try the easy way
    my $contig;
    foreach my $key ( keys %{$data} ) {
        #   print "$header\t$key\n";
        if ( $header =~ m/$key/ ) {        #this way assumes that the fasta header might have extra info
        #   print "MATCH!\n";
            $contig = $key;             
            last;
        }
    }

    #if there is no contig, then the naming conventions are off, try to salvage
    if ( defined $contig ) {
        return $contig;
    }
    else {

        warn "No match for $header attempting 10mer contig match \n";
        #die "Duplicated match to $contig for $header.\n" if (defined $data->{$contig}{length});
        
        my $max_score = 0;
        my $max_key;

        foreach my $key ( keys %{$data} ) {
            my $count = 0;
            for ( my $i = 1; $i <= length($header) - 10; $i++ ) {
                my $sub = substr($header, $i, 10);

                $count++ if ($key =~ m/$sub/);
            }

            if ( $count > $max_score ) {
                $max_score = $count;
                $max_key = $key;
            }
        }

        die "No match for $header\n" if (!defined $max_key);
        warn "Matched to $max_key\t\tScore=$max_score\n"; 
        return $max_key;
    }
}



#parse the genbank to get stats about the genes and contigs
sub parse_genbank {
    my ($data, $genome, $gbk_file) = @_;

    use Bio::SeqIO;
        
    #foreach ( keys %{$data} ) {print "Genbank = $_\n";}
    my $seqio_object = Bio::SeqIO->new( -file => $gbk_file ); 
    while ( my $seq_object = $seqio_object->next_seq() ) {
        
        #get genome description (Family and lab name if available)
        if (!defined $genome->{description} ) {

            #remove the contig data from the DEFINITION line
            my $description = ($seq_object->desc());
            #print "$description\n";
            ($description) = split ":", $description;
            #print "$description\n";
            chomp($description);
            #print "$description\n";
            $genome->{description} = $description;
        }

        #get the contig
        my $header = $seq_object->display_id();

        #figure out which contig to store it as
        my $contig = _contig_match($header, $data);
        die "Duplicated match to $contig for $header.\n" if (defined $data->{$contig}{genes});
        



        #initialize contig in case there are no genes
        $data->{$contig}{genes} = {};
        #print "$contig\n"; 
        #get the tag objects
        for my $feat_object ($seq_object->get_SeqFeatures) {
        
            #only get info from the CDS tag; include the rRNA tag incase the cds is blank.. will one entry ever have both?
                    
            if ($feat_object->primary_tag() eq "CDS" or $feat_object->primary_tag() eq "rRNA" ) {
            #if ($feat_object->primary_tag() eq "CDS" ) {
                my ($gene_tag) = $feat_object->get_tag_values('locus_tag') if ($feat_object->has_tag('locus_tag' ));
                 
                #set a reference to avoid typing out the path 
                my $gene = $data->{$contig}{genes}{$gene_tag};

                $data->{$contig}{genes}{$gene_tag}{start_pos} = $feat_object->location()->start();
                $data->{$contig}{genes}{$gene_tag}{end_pos} = $feat_object->location()->end();

                #set optional parameters
                #$gene->{mol_type} = $feat_object->get_tag_values("mol_type") ) if ($feat_object->has_tag("mol_type"));
                ($data->{$contig}{genes}{$gene_tag}{product}) = $feat_object->get_tag_values("product")   if ($feat_object->has_tag("product"));
                $data->{$contig}{genes}{$gene_tag}{length} = $data->{$contig}{genes}{$gene_tag}{end_pos} - $data->{$contig}{genes}{$gene_tag}{start_pos};
            }

            
        }
    }
    #print Dumper $data;

}

#does many system calls to eventually return a depth file
sub process_fastq {

    my ($fastq, $fasta, $path, $name) = @_;


    #check for the appropriate programs
    #system("module add bwa");
    #system("module add samtools");
    #system("module add bedtools");

    #make a directory to store info
    mkdir $path or die "Could not make directory $path\n" if ( ! -d $path );


    #index the fasta - fast enough to redo if necessary
    system ( "bwa index -p $path/$name\_database $fasta") == 0 or die "Could not index fasta $fasta\n";
            
    #align the fasta and the fastq
    if ( ! -e "$path/$name\_align.sam" ) {
        if ( $PAIRED_END ) {
            system ( "bwa aln $path/$name\_database @{$fastq}[0] > $path/$name\_1.sai") == 0 or die "Could not align 1 \n";
            system ( "bwa aln $path/$name\_database @{$fastq}[1] > $path/$name\_2.sai") == 0 or die "Could not align 2 \n";
            system ( "bwa sampe $path/$name\_database $path/$name\_1.sai $path/$name\_2.sai @{$fastq}[0] @{$fastq}[1] > $path/$name\_align.sam" ) == 0 or die "Could not align fasta $fasta with fastqs\n";
        }
        else {  #SINGLE END

            system ( "bwa mem $path/$name\_database @{$fastq}[0] > $path/$name\_align.sam" ) == 0 or die "Could not align fasta $fasta with fastq @{$fastq}[0]\n";
        }
    }
            
    #convert to BAM
    if ( ! -e "$path/$name\_align.bam" ) {
        system ( "samtools view -S -b -o $path/$name\_align.bam $path/$name\_align.sam" ) == 0 or die "Could not convert SAM to BAM: $name\n";
    }

    #sort the BAM
    if ( ! -e "$path/$name\_align_sorted.bam" ) {
        system ( "samtools sort $path/$name\_align.bam $path/$name\_align_sorted" ) == 0 or die "Could not sort BAM: $name\n";
    }
    
    #copy the fasta
    if ( ! -e "$path/$name\_copy.fasta" ) {
        system ( "cp $fasta $path/$name\_copy.fasta" ) == 0 or die "Could not copy fasta\n";
    }

    #index the fasta and make the genome file
    if ( ! -e "$path/$name\_copy.fasta.fai" ) {
        system ( "samtools faidx $path/$name\_copy.fasta" ) == 0 or die "Could not call samtools faidx $name\n";

        system ( "cut -f 1-2 $path/$name\_copy.fasta.fai > $path/$name.genome" ) == 0 or die "Couldn't make the genome file for $name\n";
    }


    #make the depth file
    if ( ! -e "$path/$name\_depth.txt" ) {
        system ( "bedtools genomecov -d -ibam $path/$name\_align_sorted.bam -g $path/$name.genome > $path/$name\_depth.txt" ) == 0 or die "Couldn't make bedtools depth file for $name\n";
    }

    my $depth = "$path/$name\_depth.txt";
    return $depth;

    ####generate a samtools depth file
    ####the samtools suite doesn't work because it doesn't report 0 coverage bases

    #this takes 100s (-A include anomalous; -B disable BAQ (save time); -Q base quality; -d max depth; -D report depth
    #system "samtools mpileup -A -BQ0 -d10000000 -D -f $fasta $path/align_sorted.bam > $path/mpileup.txt" == 0 or die "Could not run mpileup\n";

    #this takes 40s and seems to give the same depth output
    #if ( ! -e "$path/$name\_depth.txt" ) {
    #    system ( "samtools depth $path/$name\_align_sorted.bam > $path/$name\_depth.txt" ) == 0 or die "Could not run depth: $name\n";
    #}

}
