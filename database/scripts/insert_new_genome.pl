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


#insert genomes

#prepare the SQL statement to only compile once
my $G_insert = $dbh->prepare("INSERT IGNORE INTO Genomes (name, size, med_cov, avg_cov, code_density) VALUES (?, ?, ?, ?, ?);");
open my $IN, "<", "$prefix/genome_stats.txt" or die "Couldn't open genome file $prefix/genome_stats.txt\n";
while ( my $line = readline $IN ) {

    #skip if header line
    next if ($. == 1);
    
    chomp $line;

    $G_insert->execute((split "\t", $line))  or die "SQL ERR: $DBI::errstr\n";


}
close $IN;





#prepare a query with placeholders for genome and contig name
my $C_select_G = $dbh->prepare("SELECT id FROM Genomes WHERE name = ?;");

#insert contigs
my $C_insert = $dbh->prepare("INSERT IGNORE INTO Contigs (name, size, med_cov, avg_cov, code_density, gc_content, ref_genome) VALUES (?, ?, ?, ?, ?, ?, ?);");
open $IN, "<", "$prefix/contig_stats.txt" or die "Couldn't open contigs file $prefix/contig_stats.txt\n";

while ( my $line = readline $IN ) {
    #skip if header line
    next if ($. == 1);
    
    chomp $line;
    my @fields = split "\t", $line;

    my $genome = pop @fields;


    $C_select_G->execute($genome);
    my $n = 1;
    my $id;
    while ( my @row = $C_select_G->fetchrow_array() ) {
        die "Multiple genomes by the same name! : $genome\n" if $n == 2;
        $id = $row[0];
        $n++;
    }

    #don't allow the insertion of contigs without a proper reference genome
    if ( ! defined $id ) {
        warn "Genome $genome not found for entry $line\n";
        next;
    }
    push @fields, $id;
    
    $C_insert->execute(@fields);

}
close $IN;


#insert genes
open $IN, "<", "$prefix/gene_stats.txt" or die "Couldn't open genes file $prefix/gene_stats.txt\n";

#prepare a query with placeholders for genome and contig name
my $Ge_select_C = $dbh->prepare("SELECT id FROM Contigs WHERE ref_genome = ? AND name = ?;");
my $Ge_select_G = $dbh->prepare("SELECT id FROM Genomes WHERE name = ?;");
my $Ge_insert = $dbh->prepare("INSERT IGNORE INTO Genes (name, size, start_pos, end_pos, med_cov, avg_cov, product, ref_Genome, ref_Contig) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?);");
while ( my $line = readline $IN ) {
    #skip header
    next if ($. == 1);
    
    #get the genome and contig of the gene
    my @fields = split "\t", $line;

    my $contig = pop @fields;
    my $genome = pop @fields;

    #get the genome ID associated with the name. 
    $Ge_select_G->execute($genome);
    my $n = 1;
    my $id;
    while (my @row = $Ge_select_G->fetchrow_array() ) {
        die "Multiple genomes by the same name! : $genome\n" if ($n == 2);
        $id = $row[0];
        $n++;
    }

    push @fields, $id;

    #get the contig ID associated with the genome and contig
    $Ge_select_C->execute($id, $contig);

    my @row = $Ge_select_C->fetchrow_array();
    
    #grab the first field in the row (which is the id)
    $id =  $row[0];

    #add the id as the last field
    push @fields, $id;

    $Ge_insert->execute(@fields);
}
close $IN;
