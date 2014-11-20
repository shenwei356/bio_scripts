#!/usr/bin/env perl

# Function: Delete illegal charactors of head line in fasta file before blast.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    : 2014-08-14

use strict;
use BioUtil::Seq;

die "\nUsage: $0 fasta_file [fasta_file ...]\n\n"
    unless @ARGV > 0;

while (@ARGV) {
    my $file = shift @ARGV;
    my $n    = rename_fasta_header( '[^a-z\d\s\-\_\(\)\[\]\|]', '_', $file,
        "$file.rename.fa" );
    print "$file: $n records renamed\n";
}
