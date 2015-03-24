#!/usr/bin/env perl

use strict;
use File::Basename;
use BioUtil::Seq;
use BioUtil::Util;

$0 = basename $0;
die "\nusage: $0 <seq_file> <w.size start> <w.size end> <w.size step> <slid step>\n\n"
    unless @ARGV == 5;

my ( $file_query, $win_start, $win_end, $win_step, $step ) = @ARGV;
check_positive_integer($win_start);
check_positive_integer($win_end);
check_positive_integer($win_step);
check_positive_integer($step);

die "win_start should not be larger han win_end\n"
    unless $win_end >= $win_start;

my $next_seq = FastaReader($file_query);
while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;
    my $len_seq = length $seq;

    for ( my $win = $win_start; $win <= $win_end; $win += $win_step ) {
        my $end = $len_seq - $win < 0 ? 0 : $len_seq - $win;
        for ( my $i = 0; $i <= $end; $i += $step ) {
            my $s = substr( $seq, $i, $win );
            printf ">%s_window(%d,%d)\n%s\n", $header, $i+1, $win, $s;
        }
    }
}
