#!/usr/bin/env perl

use strict;
use File::Basename;
use BioUtil::Util;

$0 = basename($0);
my $usage = <<USAGE;

Usage: $0 <threads number> <fastafile> [fastafile...]

https://github.com/shenwei356/bio_scripts

USAGE

die $usage unless @ARGV >= 2;

my $threads = shift @ARGV;

for my $file (@ARGV) {
    my $fileout = "$file.align.fa";
    my $cmd     = "clustalo -i $file -o $fileout --force --outfmt fasta --threads=$threads";
    my $fail = run($cmd);
    die "failed to run:$cmd\n" if $fail;
}
