#!/usr/bin/evn perl

# gets fastX sequences by their id

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use File::Temp qw(tempfile tempdir);
use Carp;
use Readonly;
use File::Basename;

use BioUtils::FastaIO;
use BioUtils::FastqIO;

# Subroutines #
sub _guess_type;
sub _store_names;

# Variables #
Readonly::Hash my %FA_TYPES => ('fasta' => 1,
                             'FASTA' => 1,
                             'fas' => 1,
                             'FAS' => 1,
                             'fa' => 1,
                             'FA' => 1,);
Readonly::Hash my %FQ_TYPES => ('fastq' => 1,
                                'FASTQ' => 1,
                                'fq' => 1,
                                'FQ' => 1);
my $in;
my $out;
my %names;
my %found;

# parameters #
my ($fastx_file, $type, $names_file, $out_prefix, $invert, $help);

my $options_okay = GetOptions (
    "fastx_file:s" => \$fastx_file,
    "type:s" => \$type,
    "names_file:s" => \$names_file,
    "out_prefix:s" => \$out_prefix,
    "invert" => \$invert,              # flag
    "help" => \$help,                  # flag
);

# check for input errors
if ( $help ) { pod2usage(2) }
if ( ! defined $fastx_file ) { pod2usage(2) }
if ( defined $fastx_file and ! -e $fastx_file ) { pod2usage(2) }
if ( ! defined $names_file ) { pod2usage(2) }
if ( ! defined $out_prefix ) {
    $out_prefix = "seq_subset";
}
if ( ! defined $type ) {
    $type = _guess_type($fastx_file, \$type);
}
if ( ! defined $FA_TYPES{$type} and
     ! defined $FQ_TYPES{$type} ) { pod2usage(2) }



# store all the names
_store_names($names_file, \%names);

# build the in and out objects based on $type
if ( defined $FA_TYPES{$type} ) {
    # this is a fasta file
    my $file_name = $out_prefix . "." . $type;
    $in = BioUtils::FastaIO->new({stream_type => '<', file => $fastx_file});
    $out = BioUtils::FastaIO->new({stream_type => '>', file => $file_name});
}
elsif ( defined $FQ_TYPES{$type} ) {
    # this is a fastq file
    my $file_name = $out_prefix . "." . $type;
    $in = BioUtils::FastqIO->new({stream_type => '<', file => $fastx_file});
    $out = BioUtils::FastqIO->new({stream_type => '>', file => $file_name});
}
else {
    croak "Unrecognized file type: $type";
}

while ( my $seq = $in->get_next_seq() ) {
    if ( ! $invert ) {
        if ( defined $names{$seq->get_id()} ) {
            $out->write_seq($seq);
            $found{$seq->get_id()} = 1;
        }
    }
    else {
        if ( ! defined $names{$seq->get_id()} ) {
            $out->write_seq($seq);
            $found{$seq->get_id()} = 1;
        }
    }
}

# some summary stuff about names that were not found in the fastx file
print "\n";
if ( keys %found == keys %names ) {
	print "All ids stored successfully\n";
}
else {
	print keys %found, " IDs stored successfully.\n";
	
	# Print out the ones that were not found
	print "Missing IDs: ";
	foreach ( keys %names ) {
		if ( !defined $found{$_} ) {
			print $_ . ", ";
		}
	}
	print "\n";
}



########
# Subs #
########
sub _guess_type {
    my ($fastx_file) = @_;
    
    my ($filename, $directories, $suffix) = fileparse($fastx_file);
    
    if ( ! defined $FA_TYPES{$suffix} and
         ! defined $FQ_TYPES{$suffix} ) {
        croak "Unrecognized file type: $suffix";
    }
    
    return $suffix;
}

sub _store_names {
    my ($file, $names_href) = @_;
    
    open my $NAMES, "<", "$file" or croak "Cannot open file: $file";
    
    foreach my $line ( <$NAMES> ) {
        chomp $line;
        my @vals = split /\t/, $line;
        $names_href->{$vals[0]} = 1;
    }
    
    close($NAMES);
}


__END__

# POD

=head1 NAME

fastX_get_seqs_by_id.pl - Automates version number updates


=head1 VERSION

This documentation refers to fastX_get_seqs_by_id.pl version 0.0.1


=head1 SYNOPSIS

    fastX_get_seqs_by_id.pl
        --fastx_file my_file.fasta
        [--type fasta | fas | fa | fastq | fq]
        --names_file get_these_seqs.txt
        --out_prefix my_subset.fasta
        [--invert]
        [--help]

    --fastx_file  = Path to a PERL input file
    --type  = Either a fasta or fastq file
    --names_file = a file with sequence names to keep [or exclude with --invert]
    --out_prefix = output file prefix
    --invert = use this flag to get all sequences NOT in names_file
    --help  = Prints USAGE statement


=head1 ARGUMENTS
    
=head2 --fastx_file

Path to PERL input fasta or fastq formated file
    
=head2 --type

Designates the type of fastx file.  Must be either fasta, fas, fa, fastq, or fq.
    
=head2 --names_file

File path with names to keep [or exclude with --invert].  One name on each line.

=head2 --out_prefix

Name of the output file.  The suffix will be the same as the input or --type

=head2 --invert

Pass this flag to get all the sequences NOT in --names_file

=head2 [--help]
    
An optional parameter to print a usage statement.
    

=head1 DESCRIPTION

This Perl script gets a subset of sequences from a sequence file.  It keesp the
sequences whose names are passed in via --names_file.  It can also get all the
sequences that are not in --names_file by passing the --invert flag.


=head1 CONFIGURATION AND ENVIRONMENT
    
No special configurations or environment variables needed
    
    
=head1 DEPENDANCIES

version
Getopt::Long
Pod::Usage


=head1 AUTHOR

Scott Yourstone     scott.yourstone81@gmail.com
    
    
=head1 LICENCE AND COPYRIGHT

Copyright (c) 2013, Scott Yourstone
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.


=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.


=cut
