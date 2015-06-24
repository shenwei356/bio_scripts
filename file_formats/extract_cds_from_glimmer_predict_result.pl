#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;
use File::Basename;
use BioUtil::Seq;
use BioUtil::Util;

$0 = basename($0);
my $usage = qq(
    usage: $0 <glimmer .predict file> <genome file> [gff]

);
die $usage unless @ARGV == 2 or @ARGV == 3;
my $prfile  = shift @ARGV;
my $seqfile = shift @ARGV;
my $gtf     = shift @ARGV;

my $genome = ( values %{ read_sequence_from_fasta_file($seqfile) } )[0];

my @data = ();
my ( $genome, $name, $a, $b, $frame, $seq );
open my $fh, '<', $prfile or die "fail to open file: $prfile\n";
while (<$fh>) {
    $genome = $1 if /^>(.+)/;
    @data = split /\s+/, $_;
    next unless scalar(@data) == 5;
    ( $name, $a, $b, $frame ) = @data;
    next unless $a =~ /^\d+$/;

    if ( $gtf eq 'gff' ) {
        my $strand = $frame > 0 ? '+' : '-';
        printf "%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n",
             $name, 'glimmer', 'CDS', $a, $b, '.', $strand, '.', $genome;
    }
    else {
        if ( $frame > 0 ) {
            $seq = substr( $genome, $a - 1, ( $b - $a + 1 ) );
        }
        else {
            $seq = revcom( substr( $genome, $b - 1, ( $a - $b + 1 ) ) );
        }
        print ">${name}_${a}..${b}..$frame\n", format_seq($seq);
    }

}
close $fh;
