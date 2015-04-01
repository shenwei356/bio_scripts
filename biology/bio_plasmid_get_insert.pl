#!/usr/bin/perl

use strict;
use File::Basename;
use BioUtil::Seq;

# M13
my $prefix = "AGCGGCCGCGAATTGCCCTT";
my $suffix = "AAGGGCAATTCGTTTAAACCT";

$0 = basename $0;
my $usage = <<USAGE;

usage: $0 <forward seq> <reverse seq>

USAGE

die $usage unless @ARGV == 2;

my $seqf = get_the_one_seq( shift @ARGV );
my $seqr = revcom (get_the_one_seq( shift @ARGV ) );

my $sf =  extract_insert( $prefix, $suffix, $seqf );
my $sr =  extract_insert( $prefix, $suffix, $seqr );

if ( $sf ne $sr ) {
    print "forward: $sf\nreverse: $sr\n";
    die "forward and reverse sequences are not equal!";
}

print $sf, "\n";



sub extract_insert {
    my ( $prefix, $suffix, $seq ) = @_;
    die "prefix and suffix do not match sequence!\n"
        unless $seq =~ /$prefix(.+)$suffix/;
    return $1;
}

sub get_the_one_seq {
    my ($file) = @_;
    my $seqs = read_sequence_from_fasta_file($file);
    die "only one sequence should be in $file. Please check it.\n"
        unless keys %$seqs == 1;
    return ( values %$seqs )[0];
}
