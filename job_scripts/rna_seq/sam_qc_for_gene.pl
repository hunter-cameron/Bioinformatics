
### Used to parse the sam files and print tables of quality and individual OTU tables for each sample

use strict;
use warnings;

use List::MoreUtils qw(any);

sam_qc($ARGV[0]);

sub sam_qc { 
    my ($sam) = @_;
    open(my $IN, "<", $sam);
    open(my $QUAL, ">", "$sam-quality.txt");

    my $flags = {
                    # unambig mapping - in proper pair
                    '99' => ["ppaired, 2map, fwd, r1"],
                    '147' => ["ppaired, 2map, rev, r2"],

                    '83' => ["ppaired, 2map, rev, r1"],
                    '163' => ["ppaired, 2map, fwd, r2"],
                    
                    # unambig unmapped - no proper pair
                    '77' => ["npaired, 2unmap, NA, r1"],
                    '141' => ["npaired, 2unmap, NA, r2"],

                    # one unmap, other map; different strand
                    '89' => ["npaired, 1map, rev, r1"],
                    '165' => ["npaired, 1unmap, fwd, r2"],

                    '101' => ["npaired, 1unmap, fwd, r1"],   # may not be in BWA?
                    '153' => ["npaired, 1map, rev, r2"],

                    # one unmap, other map; same strand
                    '73' => ["npaired, 1map, fwd, r1"],
                    '133' => ["npaired, 1unmap, fwd, r2"],

                    '69' => ["npaired, 1unmap, fwd, r1"],
                    '137' => ["npaired, 1map, fwd, r2"],
                    
                    '181' => ["npaired, 1map, rev, r1"],    # may not be in BWA?
                    '121' => ["npaired, 1unmap, rev, r2"],
                    
                    '117' => ["npaired, 1unmap, rev, r1"],
                    '185' => ["npaired, 1map, rev, r2"],


                    # both map - not proper pair; different strands (wrong insert size, may be in diff contig)
                    '97' => ["npaired, 2map, fwd, r1"],
                    '145' => ["npaired, 2map, rev, r2"],

                    '81' => ["npaired, 2map, rev, r1"],
                    '161' => ["npaired, 2map, fwd, r2"],

                    # both map - not proper pair; same strand (wrong insert size, may be diff contig)
                    '65' => ["npaired, 2map, fwd, r1"],
                    '129' => ["npaired, 2map, fwd, r2"],

                    '113' => ["npaired, 2map, rev, r1"],
                    '177' => ["npaired, 2map, rev, r2"],

                    # both map - in proper pair; same strand   (within insert size, wrong orientation)
                    '67' => ["paired, 2map, fwd, r1"],
                    '131' => ["paired, 2map, fwd, r2"],

                    '113' => ["paired, 2map, rev, r1"],
                    '177' => ["paired, 2map, rev, r2"],
                };



    my $counts = {};

    while ( my $line = readline $IN ) {
        next if $line =~ m/^@/;     # skip header lines

        my @fields = split("\t", $line);

        my $query = $fields[0];
        my $flag = $fields[1];
        my $ref = $fields[2];
        my $pos = $fields[3];
        my $mapq = $fields[4];
       

        my $length = length($fields[9]);
        my $mismatches = $1 if ( $line =~ m/NM:i:(\d+)/ );

        $flags->{$flag}[1]++;
      
        # only count reads where both mapped to have unambiguous results
        if ( $flag == 83 or $flag == 99 ) {             # read and mate mapped; only count first read
                                                        # 99 = 1pe-fwd; 147 = 2pe-rev
                                                        # 83 = 1pe-rev; 163 = 2pe-fwd

            next if ( $length < 145 );                  # skip alignments that aren't long enough
            next if ( $mismatches  > 0 );               # skip reads with too many mismatches 
            next if ( $mapq < 30 );                     # skip reads with bad quality

            # $previous = $query;
            print $QUAL "$query\t$ref\_p$pos\t$mapq\n";
            $counts->{$ref}++;
        }
    }
    print "Flags:\n";
    foreach (sort {$flags->{$b}[1] <=> $flags->{$a}[1]} keys %{$flags}) {
        if ( ! defined (my $flag = $flags->{$_}[0]) ) {
            if ( $_ >= 2048 ) {
                $flags->{$_}[0] = "supplemental";
            }
            else {
                $flags->{$_}[0] = '';
            }
        }
        print "\t$_\t$flags->{$_}[1]\t$flags->{$_}[0]\n";
    }

    $counts->{'unmapped'} = $flags->{'77'}[1];
    print_counts($counts, "$sam-counts.txt");

}


sub print_counts {
    my ($counts, $out) = @_;

    open my $OUT, ">", $out;

    print $OUT "RefId\tCount\n";
    foreach my $ref ( keys %{$counts} ) {
        print $OUT "$ref\t$counts->{$ref}\n";
    }
}
