#!/usr/bin/env perl

# Copyright 2014 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
# https://github.com/shenwei356/bio_scripts/

use strict;

use BioUtil::Seq;
use BioUtil::Util;

my $usage = <<USAGE;

usage: extract_sequence_randomly.pl <probability> [fasta file ...]

USAGE
die $usage unless @ARGV >= 2;

my $p = shift @ARGV;
die "Probability should between 0 and 1\n"
    unless $p =~ /^[\d\.]+$/
    and $p > 0
    and $p <= 1;

my $n = 0;
for my $file (@ARGV) {
    srand();

    my $next_seq = FastaReader( $file, 1 );
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;

        next unless rand() < $p;
        $n++;
        print ">$header\t$seq";
    }
}

print STDERR "sum: $n\n";
