#!usr/bin/perl

use strict;
use warnings; 


####
####Setup and error handling
####


#Check and intake arguments for source dir and output dir
if (@ARGV != 2) {
    die("Usage: plasgraph.pl <source directory> <output directory> \n");
}
my ($src_directory, $out_directory) = @ARGV;


#try to open the source directory
opendir(my $DIR, $src_directory) or die "Directory $src_directory not found.\n";

#store all file names to an array and close directory
my @files = readdir $DIR;
my @not_fasta;
closedir $DIR;

#check if output directory exists, if not make it
#could use to put this in positive conditional
unless(-d $out_directory) {
    mkdir($out_directory) or die "Output directory $out_directory could not be created. \n";
    print("Directory created \n");
}


####
####Begin loop to process each file
####

my $file_n = 1;
foreach my $file (@files) {
    #update user on progress also counts non-fasta files...
    print("Processing file ",  $file_n++, " of ", scalar @files, "... \n");

    ###check if fasta format by looking for ">" in header
	open(my $FASFILE, "<", "$src_directory/$file");
	my $header = readline($FASFILE);

    #if the file is a fasta file proceed, consider moving this to a method
    if (defined($header) and substr($header, 0, 1) eq ">") {   
	    #format output file by removing extension
        (my $noext = $file) =~ s{\.[^.]+$}{};

        #check if log has already been created, if so skip the file. 
    	unless(-e "$out_directory/$noext.cBar.txt") {
          
          
          #call cBar and pass source file and output directory. System call outputs cBar stuff    
    	    #system("perl", "cBar.1.2/cBar.pl", "$src_directory/$file", "$out_directory/$noext.cBar.txt");
    	
            unless(-d "$out_directory/cBar_logs") {
                mkdir("$out_directory/cBar_logs") or die("Directory cBar_logs could not be created");
            }
            #use this call rather than system because it stops output from cBar
           my $errs =  `perl cBar.1.2/cBar.pl $src_directory/$file $out_directory/cBar_logs/$noext.cBar.txt`;
           print($errs);

            #set up directory for graphs
            unless(-d "$out_directory/Graphs") {
                mkdir("$out_directory/Graphs") or die "Directory Graphs could now be made";
            }

            #call R with the input file, file name, and output file
            `Rscript pgraph.r $out_directory/cBar_logs/$noext.cBar.txt $noext $out_directory/Graphs/`;
            }

    }
	else {
	    push(@not_fasta, "$src_directory/$file");
	}
	close($FASFILE);
}	


#Print out the files within the directory that were determined not to be fasta files. 
print("\n\nThe following files were not processed through cBar because they were determined not to be fasta files:\n");
foreach my $norun (@not_fasta) {
    print("$norun \n");
}
