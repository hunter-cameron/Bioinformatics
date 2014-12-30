#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long;

use DBI;
use Cwd qw(getcwd);

my $USAGE = "USAGE = insert_new_genome.pl <file>\n" .
            "file =\tusername\n\tpassword\n";

            

#get username and password
my ($file) = @ARGV;

open my $I, "<", $file or die "Could not open file:$file.\n\n$USAGE\n"; 
my @creds = <$I>;

my ($user, $password) = @creds;
chomp $user;
chomp $password;


my $prefix = getcwd();



# Connect to the database.

my $dbh = DBI->connect("DBI:mysql:database=dangl_lb;host=localhost","$user", "$password", {'RaiseError' => 1});


#read all the files to delete into an array
open my $IN, "<", $delete_file;
#my @to_delete = <$IN>;

my $G_select = $dbh->prepare("SELECT * FROM Genomes WHERE id = ?;");
my $C_select = $dbh->prepare("SELECT * FROM Contigs WHERE ref_Genome = ?;);
my $Ge_select = $dbh->prepare("SELECT * FROM Genes WHERE ref_Genome = ?;);

#need to test this one
my $delete = $dbh->prepare("DELETE FROM Genomes, Contigs, Genes WHERE Genomes.id = ? AND Contigs.ref_Genome = Genomes.id AND Genes.ref_Genome = Genomes.id;");

#print headers to the back up files
open my $G_OUT, ">", "";
open my $C_OUT, ">", "";
open my $Ge_OUT, ">", "";


###check these fields to make sure they are right
print $G_OUT "name\tdescription\tsize\tmed_cov\tavg_cov\tcode_density\n";
print $C_OUT "name\tsize\tmed_cov\tavg_cov\tgc_content\tcode_density\tref_Genome\n";
print $Ge_OUT "name\tsize\tstart_pos\tend_pos\tmed_cov\tavg_cov\tproduct\tref_Genome\tref_Contig\n";

while ( my $line = readline $IN ) {

    chomp $line;
    
    my ($id) = split "\t", $line;       #grab the id

    #execute the search
    $G_select->execute($id);
    $C_select->execute($id);
    $Ge_select->execute($id);

    #store the results
    my @genomes = $G_select->fetch_row_array();
    my @contigs = $C_select->fetch_row_array();
    my @genes = $Ge_select->fetch_row_array();


    #alter the results to get the appropriate data
    

    #print the results to the backup
    foreach (@genomes) { print $G_OUT $_, "\n" }
    foreach (@contigs) { print $C_OUT $_, "\n" }
    foreach (@genes)  { print $Ge_OUT $_, "\n" }

    #$delete->execute($id);
    
}

