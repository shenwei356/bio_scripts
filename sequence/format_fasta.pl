#!/usr/bin/env perl
use strict;
use BioUtil::Seq;

die "$0 <seq.fa> \n" unless @ARGV == 1;
my $seqfile = shift @ARGV;

my $next_seq = FastaReader($seqfile);
while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;
    $seq =~ s/\-+//g;
    printf ">%s\n%s", $header, format_seq($seq);
}
