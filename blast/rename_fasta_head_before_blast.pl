#!/usr/bin/env perl

# Function: Delete illegal charactors of head line in fasta file before blast.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    : 2014-08-14

use strict;

# try to use BioUtil::Seq
if ( eval { require BioUtil::Seq; 1; } ne 1 ) {
    die "\nPlease install BioUtil::Seq by CPAN:\n"
        . "  cpan install BioUtil::Seq\n\n";
}
else {
    BioUtil::Seq->import();
}

die "\nUsage: $0 fasta_file [fasta_file ...]\n\n"
    unless @ARGV > 0;

while (@ARGV) {
    my $file = shift @ARGV;
    my $n    = rename_fasta_header( '[^a-z\d\s\-\_\(\)\[\]\|]', '_', $file,
        "$file.rename.fa" );
    print "$file: $n records renamed\n";
}
