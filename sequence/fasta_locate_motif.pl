#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;
use File::Basename;
use BioUtil::Seq;

$0 = basename($0);
die "\nusage: $0 <motif fasta file> <subject fasta file>\n\n"
    unless @ARGV == 2;

my $queries = read_sequence_from_fasta_file( shift @ARGV );

my $next_seq = FastaReader( shift @ARGV );

print "subject\tquery\tstart\tend\tstrand\tmatched\n";
while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;

    for my $qname ( sort keys %$queries ) {
        my $qseq = $$queries{$qname};

        my $qseq_r = degenerate_seq_to_regexp($qseq);
        my $sites = degenerate_seq_match_sites($qseq_r, $seq);
        for my $site (@$sites) { 
            my ($start, $end, $matched) = @$site;           
            print "$header\t$qname\t$start\t$end\t+\t$matched\n";
        }

        my $qseq_r = degenerate_seq_to_regexp(revcom($qseq));
        my $sites = degenerate_seq_match_sites($qseq_r, $seq);
        for my $site (@$sites) { 
            my ($start, $end, $matched) = @$site;           
            print "$header\t$qname\t$start\t$end\t-\t".revcom($matched)."\n";
        }
    }
}
