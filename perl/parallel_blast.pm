package parallel_blast;

#used to run a blast search on multiple forks

use strict;
use warnings;

use List::MoreUtils qw(any);

use Data::Dumper;


sub parallel_blast {

    my ( $args_href ) = @_;

    #would it be easier to just send blast command and number of forks and parse everything out of command?
    
    
    my $command = $args_href->{command};
    my $num_forks = $args_href->{num_forks};
    my $file = $1 if ( $command =~ m/-query\s([^\s]+)/ );   #match everything after -query until the next space
    my $out = $1 if ( $command =~ m/-out\s([^\s]+)/ );     #match everything after -out until the next space
    
    print($file);


    if ( any {! defined} $num_forks, $command, $file, $out ) {
        die "Could not parallelize blast.\n";
    
    }

    my $prefix = $1 if ($out =~ m{(\/.*\/)});   #match anything between 2 /... will break if trailing / is left off;


    my $contigs = get_contigs($file);


    # make the bins
    my $bin_size = (length(keys $contigs) // $num_forks) + 1

    my $counter = 1
    foreach my $contig ( keys $contig ) {
        if $counter <= $bin_size

        $counter++

    my $num_files = bsub_blast($contigs, $command, $prefix);

    wait_all_jobs();


    #write the results to an output file; it it becomes important to have the results in exactly the same order, the contigs could be read in as an array
    open my $OUT, ">", $out;
    for ( my $i = 1; $i <= $num_files; $i++ ) {

        open my $IN, "<", "$prefix/blst_temp/$i.txt";

        while (readline $IN) {print $OUT $_}

        close $IN;
    }


    #cleanup
    #print "system(rm -rf $prefix/blast_temp/ \n";
    system("rm -rf $prefix/blst_temp/");
}


#minimalistic sub to turn a fasta into a hash
sub get_contigs {
    my ($file) = @_;

    open my $IN, "<", $file;

    my %contigs; 
    my $header;
    while (my $line = readline $IN) {

        chomp $line;

        if ( $line =~ m/\A>/ ) {

            $header = $line;
        
        }
        else {
            $contigs{$header} .= $line;
        }
    }

    return \%contigs;
}



sub bsub_blast {

    my ($contigs, $command, $prefix) = @_;


    my $temp_dir = "$prefix/blst_temp";

    (-d $temp_dir) or mkdir $temp_dir;
    
    my $id = 1;
    foreach my $contig ( keys %{$contigs} ) {

        my $temp_fas = "$temp_dir/$id.fasta";
        my $temp_result = "$temp_dir/$id.txt";
        open my $OUT, ">", $temp_fas;

        print $OUT $contig . "\n" . $contigs->{$contig} . "\n";

        close $OUT;

        my $sub_command = substitute_command($command, $temp_fas, $temp_result);
        
        my $bsub = join " ", ( "bsub",
                                "-q week",
                                "-o $prefix/multi-out.bjobs",
                                "-e $prefix/multi-err.bjobs",
                                "-n 1",
                                "-J parallelBLAST",
                                "$sub_command"
                                );

        #system($bsub) == 0 or die "Could not use bsub: \n $bsub\n";


        $id++;
    }
    #return the number of files
    return $id - 1;
}


sub wait_all_jobs {

    my $is_running = 1;


    #checks every 30s to see if the any job is still running based on the file size of the bjobs output
    #might fail w/o write permissions
    while ( $is_running ) {
        sleep(30);

        system("bjobs -J parallelBLAST > ACTIVE_JOBS?");
        $is_running = -s "ACTIVE_JOBS?";
    }
    system("rm ACTIVE_JOBS?");

}




sub substitute_command {

    my ($command, $sub_file, $sub_out) = @_;

    my $query = $1 if ( $command =~ m/-query\s([^\s]+)/ );  #get query
    my $out = $1 if ( $command =~ m/-out\s([^\s]+)/ );      #get out file
    
    #print "$command\n\n$query\t$sub_file\n\n$out\t$sub_out\n";

    #replace the original query and outs with a subfile
    $command =~ s/$out/$sub_out/;
    $command =~ s/$query/$sub_file/;
   

    return $command;
}

















1;

