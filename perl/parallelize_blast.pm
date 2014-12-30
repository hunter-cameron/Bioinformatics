package parallelize_blast;

#used to run a blast search on multiple forks

use strict;
use warnings;

use List::MoreUtils qw(any);

use Data::Dumper;


sub parallelize_blast {

    my ( $args_href ) = @_;

    #would it be easier to just send blast command and number of forks and parse everything out of command?
    
    
    my $command = $args_href->{command};
    my $num_forks = $args_href->{num_forks};
    my $file = $1 if ( $command =~ m/-query\s([^\s]+)/ );   #match everything after -query until the next space
    my $out = $1 if ( $command =~ m/-out\s([^\s]+)/ );     #match everything after -out until the next space
    


    if ( any {! defined} $num_forks, $command, $file, $out ) {
        die "Could not parallelize blast.\n";
    
    }
    #open the file once to get a line count
    open my $IN, "<", $file;

    while (readline $IN) {}

    my $num_lines = $.;

    close $IN;

    #each sequence has 2 lines (header + seq)
    my $sequences = $num_lines / 2;

    #get sequences per fork, add 1 because in scales down (last fork may have fewer rather than last fork holding remainder.)
    my $sequences_per_fork = int( $sequences / $num_forks ) + 1;


    #open the file again to rewrite
    open $IN, "<", $file;
    my $n = 0;
    my $sub_file = 0;
    my $OUT;
    my @sub_files;
    while ( my $line = readline $IN ) {
        chomp $line;
        #if no lines have been written, open the output file
        if ( $n == 0 ) {
            open $OUT, ">", "$out\_s$sub_file.tfas" or die "Could not open $out\_s$sub_file.tfas for writing\n";
            push @sub_files, "$out\_s$sub_file.tfas";
            $sub_file++;
        }

        print $OUT "$line\n";

        $n++;

        #reset n if the appropriate number of lines have been written
        if ( $n == $sequences_per_fork * 2 ) {
            $n = 0;
            close $OUT;
        }
    }
    close $IN;



    #print "Number of files: " . scalar @sub_files . "\n";

    my @children = _init_forks(\@sub_files, $command);

    #print "Number of children: " . scalar @children . "\n";
    #print Dumper @children; 
    #wait for the children to finish
    foreach (@children) {
        my $tmp = waitpid($_, 0);
        print "Done with pid $_\n";
    }


    #combine the subfiles into a large one and clean up
    open $OUT, ">", $out or die "Could not print results\n";
    foreach my $file ( @sub_files ) {

        #sub the extension to the blast results
        unlink $file;       #delete the temp fasta
        $file =~ s/.tfas/.ttxt/;

        open my $IN, "<", $file or die "Could not open $file for reading\n";
        while ( readline $IN ) {
            print $OUT $_;
        }
        close $IN;
        unlink $file;       #delete the temp output
        
    }

    close $OUT;

}








sub _init_forks {

    my ($files_ref, $command) = @_;

    my @files = @{$files_ref};


    my @children;
    for ( my $i = 0; $i < @files; $i++ ) {

        my $pid = fork();

        if ($pid) {     #add the child pid to the parent
            push @children, $pid;
        }
        elsif ($pid == 0) {     #start the child process
            
            _blast($command, $files[$i]);
            #childproc();

        }

        else {
            die "Couldn't Fork.\n";
        }
    }

    return @children;
}


sub childproc {
    sleep 1;
    exit();
}

sub _blast {

    my ($command, $sub_file) = @_;

    #replace the temp fasta ext with the temp text
    (my $sub_out = $sub_file) =~ s/.tfas/.ttxt/;
    
    

    my $query = $1 if ( $command =~ m/-query\s([^\s]+)/ );  #get query
    my $out = $1 if ( $command =~ m/-out\s([^\s]+)/ );      #get out file
    
    #print "$command\n\n$query\t$sub_file\n\n$out\t$sub_out\n";

    #replace the filler query with a subfile
    #order is important here; can this substitution still go wrong?
    $command =~ s/$out/$sub_out/;
    $command =~ s/$query/$sub_file/;
   

    #print "\n\n$command\n";
    system($command) == 0 or die "System command: $command failed. \n";
    
    #exit the thread
    exit;
}


        
















1;

