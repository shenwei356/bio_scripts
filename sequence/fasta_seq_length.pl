#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;

use BioUtil::Seq;

# get the file list
my @files = ();
for my $file (@ARGV) {
    for my $f ( glob $file ) {
        push @files, $f;
    }
}
if ( @files == 0 ) {
    push @files, 'STDIN';
}

for my $file (@files) {
    my $next_seq = FastaReader($file);
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;

        my $len = length($seq);
        print "$header\t$len\n";
    }
}
