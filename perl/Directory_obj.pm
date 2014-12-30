package Directory_obj;

use strict;
use warnings;

use List::MoreUtils qw(any);
use Class::Std::Utils;

{   


    #Class Attributes
    my %files_of;
    my %file_index;
    my %file_ext;
    
    sub new {
        my ($class, $args_href) = @_;
        
        #check if essential parameters defined, could use to give better error report
        if ( any {! defined} $args_href->{directory} ) {
            die "Cannot create new folder object. Essential parameter undefined.\n";
        }


        #bless takes two args, a reference to the variable and a string containing th ename of the class
        #/do{my $anon_scalar} = reference to a scalar that only has scope within the do statement. So the object doesn't "exist" after the end of the do but it still kept alive because of the blessed reference to it
        #could be done with the anon_scalar() method from class::std::utils
        my $new_obj = bless \do{my $anon_scalar}, $class;

        #set parameters for object; ident() returns a unique id for the object
        $file_ext{ident $new_obj} = $args_href->{ext} if ( defined $args_href->{ext} );
        $files_of{ident $new_obj} = _get_files($args_href->{directory}, $args_href->{ext});

        #set the file index to -1 so when it is incremented in the first call of next_file if will be 0
        #$file_index{ident $new_obj} = -1;
        

        return $new_obj;
    
    }

    sub get_files {
        my ($self) = @_;
        
        return @{$files_of{ident $self}};
    }

    sub next_file {
        my ($self) = @_;
        if ( ! defined $file_index{ident $self} ) {
            $file_index{ident $self} = 0;
        }
        else {
        $file_index{ident $self}++;
        }

        if ( $file_index{ident $self} >= @{$files_of{ident $self}} ) {
            return 0;
        }
        
        else {
            return $files_of{ident $self}->[$file_index{ident $self}];
        }
    }


    sub get_filename {
        my ($self, $index) = @_;

        #if the index wasn't sent, assume the index on file
        $index = $file_index{ident $self} if ( ! defined $index );

        my $fullpath = $files_of{ident $self}->[$index];
        chomp($fullpath);
        my $filename;
        if ( $fullpath =~ m{/*           #start at last / in folder (optional)
                            ([\w\.\-]+)     #match and capture one or more word characters, periods, or hyphens (filename)
                            \.[\w]+\z    #match a . and one or more word characters and then the end of line (extension)
                            }xms ) {
            $filename = $1;
        }
        
        return $filename;
    }

    sub get_number_of_files {
        my ($self) = @_;
        return scalar @{$files_of{ident $self}};
    }

    sub _get_files {
        my ($folder, $file_ext) = @_;

        #check if path is a folder; if not, return the lone file 
        if ( ! -d $folder ) {
            my @lone_file = ($folder);
            return \@lone_file;
        }

        #read all the files
        opendir my $DIR, $folder or die "Could not open $folder. Check permissions\n";
        my @files = readdir $DIR;

        #parse the list of files
        for  ( my $i = 0; $i < @files; $i++ ) {
            
            #remove if it begins with a .
            if ( $files[$i] =~ m/\A[.]/ ) {
                splice(@files, $i, 1);
                redo;       #redo because the index went down by one for the splice
            }
            
            #remove if doesn't match specified extension
            if ( defined $file_ext and $files[$i] !~ m/$file_ext\z/ ) {
                splice(@files, $i, 1);
                redo;
            }

            #remove if it is another directory
            if ( -d $files[$i] ) {
                splice(@files, $i, 1);
                redo;
            }

            #otherwise, append the folder name to the beginning for a full path
            $files[$i] = $folder . $files[$i];

        }

        return \@files;



    }



}

1;
