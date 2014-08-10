#!/usr/bin/env perl
# Copyright 2013 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
use Getopt::Long;
use Digest::MD5;
use strict;

my $usage = <<"USAGE";
===============================================================================
Function: Find common sequences in fasta files. Version 2.
          Features:
              1) Comparing by name or sequence are both supported.
              2) No files number limit.
              3) Very low RAM usage. (Lower than Version 1).
          Note that:
              1) Records with different names may have same sequences.
              2) Case of sequence letters or name may be different.
              3) Duplicated records may exist in a fasta file.
Contact : Wei Shen <shenwei356#gmail.com>
Date    : 2013-12-05
Update  : 2014-08-10
Site    : https://github.com/shenwei356/bio_scripts

Usage   : $0 [-s] [-i] fastafile fastafile2 [fastafile3 ...]
Options :
   -s Comparing by sequence.
   -i Ignore case.
===============================================================================

USAGE

our $by_seq      = 0;
our $ignore_case = 0;
GetOptions(
    "s" => \$by_seq,
    "i" => \$ignore_case,
) or die $usage;

# at least two files;
die "$usage\n>2 sequence file needed.\n" unless @ARGV >= 2;

our $counts = {};
our $names  = {};

our ( $file, $has_head, $last_head, $head, $head0, $seq_len, $seq_md5 );
our $md5;
$md5 = Digest::MD5->new if $by_seq;

# check files
for $file (@ARGV) {
    die "File ($file) does not exists.\n" unless -e $file;
}

for $file (@ARGV) {
    open IN, "<", $file
        or die "Fail to open file: $file!\n";

    $has_head = 0;
    $seq_len  = 0;
    $md5->reset if $by_seq;

    while (<IN>) {
        s/\r?\n//;
        if (/^\s*>/) {    # fasta head
            s/>\s*//;
            s/\s+$//;

            recording();

            $seq_len  = 0;
            $has_head = 1;
        }
        elsif ( $has_head == 1 ) {    # sequence          
            next if $_ eq "";

            $seq_len += length $_;

            next unless $by_seq;
            tr/A-Z/a-z/ if $ignore_case;
            $md5->add($_);
        }
    }
    close IN;

    # do not forget the last record
    recording() if $seq_len > 0;
}

sub recording {
    $head0     = $last_head;    # orgin sequence name
    $last_head = $_;            # store this head for next turn;

    $head = $head0;
    $head = lc $head if $ignore_case;
    if ($by_seq) {
        $seq_md5 = $md5->hexdigest;
        $md5->reset;

        # ingore sequence records without sequence.
        return if $seq_len == 0;

        # count sequences with md5 $seq_md5 in $file
        $$counts{$seq_md5}{$file}++;

        # record the origin sequence name.
        $$names{$seq_md5}{$file} = $head0;
    }
    else {
        # ingore sequence records without head
        return if $head eq '';

        # count sequences with name $head in $file
        $$counts{$head}{$file}++;
        $$names{$head}{$file} = $head0;
    }
}

# find common sequences
my $file_num = scalar @ARGV;

# extract sequences from the first file.
$file = $ARGV[0];
my $names_ok = {};
for my $key ( keys %$counts ) {

    # all files have a same record
    next unless ( scalar keys %{ $$counts{$key} } ) == $file_num;

    # save into a hash.
    $$names_ok{ $$names{$key}{$file} }
        = $$counts{$key}{$file};
}

# print common sequences
my $is_target = 0;
open IN, "<", $file
    or die "Fail to open file: $file!\n";
while (<IN>) {
    if (/^\s*>/) {
        s/>\s*//;
        s/\s+$//;
        next if $_ eq '';

        $head      = $_;
        $is_target = 0;
        if ( exists $$names_ok{$head} and $$names_ok{$head} > 0 ) {
            print ">$head\n";
            $is_target = 1;

            # just export one record for duplicated records.
            $$names_ok{$head} = 0;
        }
    }
    elsif ( $is_target == 1 ) {
        print $_;
    }
}
close IN;
