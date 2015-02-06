#!/usr/bin/env perl
use strict;

my $usage = <<USAGE;

Function: Reset start position for circular genome.
   Usage: reset_start_position_for_circular_genome <fasta file> <new start>
 Example: 
    1. Set the 100th base as the new start position
        reset_start_position_for_circular_genome seq.fa 100

Author: Wei Shen <shenwei356#gmail.com> <http://shenwei.me>
Change history:
    - 2014-04-30 rewrite.
    - 2011 first edition.

USAGE

die $usage unless @ARGV == 2;

my ( $infile, $newstart, $head, $seq, $newseq, $buffer, $outfile );

$infile   = shift;
$newstart = shift;

die "newstart should be integer greater than 0, you input $newstart.\n"
    unless $newstart =~ /^\d+$/ and $newstart > 0;

$buffer = '';
open IN, $infile or die "fail to open sequence file $infile!\n";
local $/ = '>';
<IN>;

while (<IN>) {
    s/>$//;
    ( $head, $seq ) = split "\r?\n", $_, 2;
    $seq =~ s/\s+//g;

    $newseq = substr( $seq, $newstart - 1  ) . substr( $seq, 0, $newstart - 1 );
    
    $buffer .= ">$head (start position move to $newstart)\n"
        . format_seq( $newseq, 70 ) . "\n";
}
close IN;
$/ = "\n";

$outfile = "$infile.newstart$newstart.fa";
if ( $infile =~ /(.+)\.(.+?)$/ ) {
    $outfile = "$1.newstart$newstart.$2";
}
open OUT, ">", $outfile or die "failed to open file $outfile\n";
print OUT $buffer;
close OUT;

sub format_seq($$) {
    my ( $s, $n ) = @_;
    my $s2 = '';
    my ( $j, $int );
    $int = int( ( length $s ) / $n );
    for ( $j = 0; $j <= $int - 1; $j++ ) {
        $s2 .= substr( $s, $j * $n, $n ) . "\n";
    }
    $s2 .= substr( $s, $int * $n );
    return $s2;
}
