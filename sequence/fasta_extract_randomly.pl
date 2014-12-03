#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;

use File::Basename;
use BioUtil::Seq;
use BioUtil::Util;

$0 = basename($0);
my $usage = <<USAGE;

Randomly extract fasta sequences by a given proportion

Examples: 

    1) $0 0.1 seq.fa
    2) $0 0.1 seq.fa seq2.fa   # multi files supported
    3) $0 0.1 seq*.fa          # glob expression
    4) cat seq.fa | $0 0.1     # read from STDIN

https://github.com/shenwei356/bio_scripts

USAGE
die $usage unless @ARGV >= 1;

my $p = shift @ARGV;
die "Probability should between 0 and 1\n"
    unless $p =~ /^[\d\.]+$/
    and $p > 0
    and $p <= 1;

srand();

my @files = file_list_from_argv(@ARGV);

my $n = 0;
for my $file (@files) {
    my $next_seq = FastaReader( $file, 1 );
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;

        next unless rand() < $p;
        $n++;
        print ">$header\n$seq";
    }
}

print STDERR "sum: $n\n";
