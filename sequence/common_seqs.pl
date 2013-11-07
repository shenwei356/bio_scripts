#!/usr/bin/env perl
# Copyright 2013 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
use Getopt::Long;
use Digest::MD5 'md5_hex';
use strict;

my $usage = <<"USAGE";
===============================================================================
Fnction: Find common sequences in fasta files.
         Features:
             1) Comparing by name or sequence are both supported.
             2) No files number limit.
             3) Low RAM usage.
         Note that:
             1) Records with different names may have same sequences.
             2) Case of sequence letters or name may be different.
             3) Duplicated records may exist in a fasta file.
Contact: Wei Shen <shenwei356#gmail.com>
Date   : 2013-11-07
Site   : https://github.com/shenwei356/bio_scripts

Usage  : $0 [-s] [-i] fastafile fastafile2 [fastafile3 ...]
Options:
	-s Comparing by sequence.
	-i Ignore case.
===============================================================================

USAGE

my $by_seq      = 0;
my $ignore_case = 0;
GetOptions(
    "s" => \$by_seq,
    "i" => \$ignore_case,
) or die $usage;

# at least two files;
die "$usage\n>2 sequence file needed.\n" unless @ARGV >= 2;

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

            # count sequences with md5 $seq_md5 in $file
            $$counts{$seq_md5}{$file}++;    #
                                            # record the origin sequence name.
            $$names{$seq_md5}{$file} = $head0;
        }
        else {
            # count sequences with name $head in $file
            $$counts{$head}{$file}++;
            $$names{$head}{$file} = $head0;
        }
    }
}

# output common sequences
my $file_num = scalar @ARGV;
$file = $ARGV[0];    # extract sequences from the first file.
my $names_ok = {};
for my $key ( keys %$counts ) {

    # all files have a same record
    next unless ( scalar keys %{ $$counts{$key} } ) == $file_num;

    $$names_ok{ $$names{$key}{$file} } =
      $$counts{$key}{$file};    # save to a hash.
}

$next_seq = FastaReader($file);
while (1) {
    ( $head, $seq ) = &$next_seq();
    last
      if $head eq "" and $seq eq "";

    if ( exists $$names_ok{$head} and $$names_ok{$head} > 0 ) {
        print ">$head\n$seq\n";

        # just export one record for duplicated records.
        $$names_ok{$head} = 0;
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
