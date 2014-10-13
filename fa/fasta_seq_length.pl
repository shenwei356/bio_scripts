#!/usr/bin/env perl

# Copyright 2014 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
# https://github.com/shenwei356

use strict;

# try to use BioUtil::Seq
if ( eval { require BioUtil::Seq; 1; } ne 1 ) {
    die "\nPlease install BioUtil::Seq by CPAN:\n"
        . "  cpan install BioUtil\n\n";
}
else {
    BioUtil::Seq->import();
}

my $usage = <<USAGE;

usage: fasta_seq_length <infile> [infile ...]

USAGE

die $usage unless @ARGV > 0;

for my $file (@ARGV) {
    my $outfile = "$file.len";
    open OUT, ">", $outfile
        or die "failed to open file: $outfile\n";

    my $next_seq = FastaReader($file);
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;

        my $len = length($seq);
        print OUT "$header\t$len\n";
    }

    close OUT;
}
