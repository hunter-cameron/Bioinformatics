#!/usr/bin/perl


#Accepts fasta files or folders for refs and qrys and crosses them against each other to create a dotplot. 
#Begins making the output structure in the location specified by --prefix
#Mimics output structure of the input directorys
#
#Wildcard * treats all files in the directory as a query
#
#Example: to cross one fasta against a folder of others
#mummer_dotplot_by_fasta.pl --prefix /my_output/ --ref ref_fasta qry_directory/*


use warnings;
use strict;

use Getopt::Long;

my $usage = "\n\nmummer_dotplot_by_fasta.pl --ref <fasta> ... --qry <fasta>...
[--qry <fastafile>, --prefix <path>] [--Help]\n\n ";

#set GetOpt defaults
my @refs = ();
my @qrys = ();
my $PREFIX = '';

GetOptions ('ref=s' => \@refs,
            'qry=s' => \@qrys,
            'prefix=s' => \$PREFIX,

                        );


#ensure there is a reference sequence
if (! @refs) {
    die("$usage");
}



#assume any additional argument to be a query
for(my $i = 0; $i < @ARGV; $i++) {
        push(@qrys, $ARGV[$i]);
}

#make the output directory if necessary
if (! -d $PREFIX) {
    mkdir($PREFIX) or die ("Couldn't make output directory $PREFIX \n");
}

#make sure paths entered exist
&check_paths(\@refs, \@qrys);

#set the output path to the prefix the user specified
my $path = $PREFIX;


####
####Begin main loop to create dotplots
####

#loop through references
foreach my $ref (@refs) {
    #reset path to the prefix
    $path = $PREFIX;
   
   #if the path is to a directory, begin building the output structure
    my @ref_files = ();
    if (-d $ref) {

        #pass the reference directory and path and return the files (as a reference) and the modified path
        (my $ref_files_r, $path) = &is_directory($ref, $path);

        #dereference the files
        @ref_files = @$ref_files_r;
    }
    #if not a directory, store the single file in the array
    else {
        $ref_files[0] = $ref;
    }
    
    #make sure the file(s) is a fasta
    my @ref_fastas = &check_fasta(\@ref_files);
    

    #loop through all query files using the same procedure as references
    my @qry_files = ();
    foreach my $qry (@qrys) {

       if (-d $qry) {
            (my $qry_files_r, $path) = &is_directory($qry, $path);

            @qry_files = @$qry_files_r;
        }
        else {
            $qry_files[0] = $qry;
        }

        #make sure the file(s) is a fasta
        my @qry_fastas = &check_fasta(\@qry_files);

        &plot(\@ref_fastas, \@qry_fastas, $path);
    }
}





########
# SUBS #
########


####checks if paths in an array of files exist, dies if not

sub check_paths {
    #check all inputted arrays
    for (my $arg = 0; $arg < @_; $arg++) {

        #check if each path in the array exists
        foreach my $path (@{$_[$arg]}) {
            if (! -e $path) {
                die("Aborting - Path: $path does not exist\n");
            }
        }
    }
}     



####returns an array of fasta files from an array of files

sub check_fasta {
    my @files = @{$_[0]};
    my @fasta_files = ();

    #loop through all files
    foreach my $file (@files) {
        open(my $IN, "<", "$file") or die("Could not open file $file\n");
        if(readline($IN) =~ m{^>}) {     #match the opening > of a fasta file
            #push file to array if it was a fasta
            push(@fasta_files, $file);
        }
        close($IN);
    }
    return(@fasta_files);
}


####trims the prefix and extension from a file name
sub trim {
    my $original = $_[0];

    #remove prefix
    (my $trimmed = $original) =~ s{/     #begin with a /
                                 .*    #any chars
                                 /     #end with a /
                                 }
                                 {}xms;  #replace with nothing

    #remove extension
    $trimmed =~ s{\.     #begin with a period
                     [^.]+  #match one or more of not . (new line?)??
                     $      #end of line
                    }
                    {}xms;   #replace with nothing
    
    
    return($trimmed);
}


####builds output directory structure for directories; returns an array of files and a new path
sub is_directory {
    my $dir = $_[0];
    my $path = $_[1];
    
    opendir(my $DIR, $dir); 
    my @files = readdir($DIR) or die("Could not read files from $dir\n");
    closedir($DIR);
    
    #add the directory to the path for the files
    for (my $i = 0; $i < @files; $i++) {
        $files[$i] = "$dir/" . $files[$i];
    }

    #trim to the top level directory to add to output path
    (my $trim_dir = $dir) =~ m{\/         #begins with a /
                              ([^\/]+    #capture one or more non / characters
                              \/)$        #ends with a slash
                              }xms;

    $trim_dir = $1;

    #add top level directory to output path
    my $out_dir = "$path/$trim_dir/";

    #make output directory
    (-d $out_dir) or mkdir($out_dir) or die("Couldn't make output directory: $out_dir");

    &check_paths(\@files);   
    
    return(\@files, $out_dir);
}


####call nucmer and mummerplot with arguments (depends on having then in user's PATH
sub plot {
    my @refs = @{$_[0]};
    my @qrys = @{$_[1]};
    my $path = $_[2];

    foreach my $ref (@refs) {
        my $trim_ref = &trim($ref);

        foreach my $qry (@qrys) {
            my $trim_qry = &trim($qry);
            my $out_path = "$path/$trim_ref=X=$trim_qry";

            (-d $out_path) or mkdir($out_path) or die("Couldn't make output directory $out_path \n");
            
            system("nucmer -p $out_path/out $ref $qry") == 0 or die("Couldn't call nucmer. Make sure it is included in \$PATH \n");
            
            #Call mummerplot, --layout = arrange contigs to make the best line
                            # --small = output a small plot (--medium, --large)
                            # --png = output in png format (--postscript)
                            # -R = fasta file to use as ref for --layout
                            # -Q = fasta file to use as qry for --layout
                            # -p = prefix for output files
            system("mummerplot -p '$out_path/$trim_ref=X=$trim_qry' --layout --small --png -R $ref -Q $qry $out_path/out.delta") == 0 or die("Couldn't call mummerplot. Make sure it is included in \$PATH \n");
        }
    }
}



