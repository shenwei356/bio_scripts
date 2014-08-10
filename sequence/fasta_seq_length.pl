#!/usr/bin/env perl

# Copyright 2014 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
# https://github.com/shenwei356

use strict;

my $usage = <<USAGE;

usage: fasta_seq_length <infile> [infile ...]

USAGE

die $usage unless @ARGV > 0;

for my $file (@ARGV) {
    my $outfile = "$file.len";
    open OUT, ">", $outfile
        or die "failed to open file: $outfile\n";

    my $next_seq = FastaReader($file);

    my ( $head, $seq ) = ( '', '' );
    my $len = 0;

    while (1) {
        ( $head, $seq ) = &$next_seq();
        last if $head eq "" and $seq eq "";

        $len = length($seq);
        print OUT "$head\t$len\n";
    }

    close OUT;
}

# FastaReader is a fasta file parser, which returns a function that
# returns a pair of head and sequence when it was called
#
# Example:
#    use lib './FastaReader.pm';
#
#    my $next_seq = FastaReader("test.fa");
#
#    my ( $head, $seq );
#    while (1) {
#        ( $head, $seq ) = &$next_seq();
#        last
#          if $head eq "" and $seq eq "";
#        print ">$head\n$seq\n";
#    }
sub FastaReader {
    my ($file) = @_;
    open IN, "<", $file
        or die "Fail to open file: $file!\n";
    local $/ = '>';
    <IN>;
    $/ = '\n';

    my ( $line, $head, $seq );
    return sub {
        local $/ = '>';
        while ( $line = <IN> ) {
            $line =~ s/\r?\n>?$//;
            ( $head, $seq ) = split /\r?\n/, $line, 2;
            $seq =~ s/\s+//g;
            return ( $head, $seq );
        }
        close IN;
        $/ = "\n";
        return ( "", "" );
    };
}

1;
