#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts
use strict;
use BioUtil::Seq;

die "format_fasta.pl <seq.fa> \n" unless @ARGV == 1;
my $seqfile = shift @ARGV;

my $next_seq = FastaReader($seqfile);
while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;
    $seq =~ s/\-+//g;
    printf ">%s\n%s", $header, format_seq($seq);
}
