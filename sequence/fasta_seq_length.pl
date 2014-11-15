#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;

use File::Basename;
use BioUtil::Seq;

$0 = basename($0);
my $usage = <<USAGE;

Usage: $0 <infile> [infile ...]

https://github.com/shenwei356/bio_scripts

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
