#!/nas02/apps/perl-5.18.2/bin/perl


####Does bwa alignment on two folders; matches filenames

### Reccommended to run with multiple processors 
### run from the same directory to skip already processed files. 

use strict;
use warnings;

use Cwd 'getcwd';
use Getopt::Long;


use Directory_obj;

my $USAGE = "driver_coverage_pipeline.pl -fasta </fasta/> -fastq </fastq/> -gbk </gbk> [-prefix mypath] \n";

#number of threads to run
my $CPU = 10;

my $fasta_dir;
my $fastq_dir;
my $gbk_dir;
my $prefix = getcwd();
my $PAIRED_END;



GetOptions (    'fasta=s' =>    \$fasta_dir,
                'fastq=s' =>    \$fastq_dir,
                'gbk=s'   =>    \$gbk_dir,
                'prefix=s' =>   \$prefix,
                'pe|PE' =>      \$PAIRED_END,
                );


if ( !defined $fasta_dir  or !defined $fastq_dir or !defined $gbk_dir ) {
    die $USAGE;
}

my $FSA = Directory_obj->new({'directory' => $fasta_dir, 'ext' => ".fasta",});
my $FSQ = Directory_obj->new({'directory' => $fastq_dir, 'ext' => ".fastq",});
my $GBK = Directory_obj->new({'directory' => $gbk_dir,   'ext' => ".gbk",});

while ( my $fasta = $FSA->next_file() ) {

    my $name = $FSA->get_filename();

            bsub($fasta, $fastq_dir, $gbk_dir, $prefix);
            
            #exit the loop if a match was found.
            #last;    
        
    

    #print the name of the fasta if no match was found
    #print "No match for fasta: $name\n";
    #last;
}


#wait for all jobs to finish
#wait_all_jobs();


sub bsub {
    my ($fasta, $fastq, $gbk, $prefix) = @_;


    my $command = join( " ", (  "coverage_pipeline.pl",
                                "-fasta $fasta",
                                "-fastq $fastq",
                                "-gbk $gbk",
                                "-prefix $prefix"
                                ));

    $command = join( " ", (  $command, "-pe" )) if ($PAIRED_END);

    my $bsub = join( " ", ( "bsub",
                            "-q week",
                            "-o $prefix/driver-out.bjobs",
                            "-e $prefix/driver-err_\%J.bjobs",
                            "-n 1",
                            "-J parallelCov",
                            "$command"
                            ));

    system( "$bsub" ) == 0 or die "Could not use bsub\n";

    return;
}



sub wait_all_jobs {

    my $is_running = 1;

    while($is_running) {

        sleep(30);

        system("bjobs -J parallelCov > ACTIVE_JOBS?");
        $is_running = -s "ACTIVE_JOBS?";
    }

    system("rm ACTIVE_JOBS?");
}


