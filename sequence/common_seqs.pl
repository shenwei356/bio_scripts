#!/usr/bin/env perl
# Copyright 2013 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
use Getopt::Long;
use strict;

my $usage = <<"USAGE";
===============================================================================
Name   : Find common sequences in fasta files.
		 Features:
	         1) Comparing by name or sequence are both supported.
			 2) No files number limit.
	         3) Low RAM usage.
         Note that sequence names and case of sequence letters may be different
Contact: Wei Shen <shenwei356#gmail.com>

Usage  : $0 [-s] [-i] fastafile fastafile2 [fastafile3 ...]
Options:
	-s Comparing by sequence.           [false]
	-i Ignore case.                     [false]
===============================================================================

USAGE

my $by_seq      = 0;
my $ignore_case = 0;
GetOptions(
    "s" => \$by_seq,
    "i" => \$ignore_case,
) or die $usage;
die $usage unless @ARGV >= 2;    # at least two files;

if ($by_seq) {
    use Digest::Perl::MD5 'md5_hex';
}

my $file_num = scalar @ARGV;

my $counts = {};
my $names  = {};

my ( $file, $next_seq, $head, $head0, $seq, $seq_md5 );

for $file (@ARGV) {
    $next_seq = FastaReader($file);
    while (1) {
        ( $head, $seq ) = &$next_seq();
        last
          if $head eq "" and $seq eq "";

        $head0 = $head;                     # orgin sequence name
        $head = lc $head if $ignore_case;

        if ($by_seq) {
            $seq =~ tr/A-Z/a-z/ if $ignore_case;
            $seq_md5 = md5_hex($seq);
            $$counts{$seq_md5}++;
            $$names{$seq_md5}{$file} = $head0;
        }
        else {
            $$counts{$head}++;
            $$names{$head}{$file} = $head0;
        }
    }
}

# output common sequences
$file = $ARGV[0];    # export from the first file.

my $names_ok = {};   # get the sequence names in the first file.
for my $key ( keys %$counts ) {
    next unless $$counts{$key} == $file_num;    # every file has a same record
    $$names_ok{ $$names{$key}{$file} } = 1;     # save to a hash.
}

$next_seq = FastaReader($file);
while (1) {
    ( $head, $seq ) = &$next_seq();
    last
      if $head eq "" and $seq eq "";

    if ( exists $$names_ok{$head} ) {
        print ">$head\n$seq\n";
    }
}

# FastaReader is a fasta file parser, which returns a function that
# returns a pair of head and sequence when it was called
#
# Example:
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
        while (1) {
            $line = <IN>;
            last
              if $line eq "";

            $line =~ s/\r?\n>?$//;
            ( $head, $seq ) = split /\r?\n/, $line, 2;
            $seq =~ s/\s+//g;
            return ( $head, $seq );
        }
        close IN;
        $/ = "\n";
        return ( "", "" );
      }
}
