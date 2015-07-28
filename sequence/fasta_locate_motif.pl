#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;
use Getopt::Long;
use File::Basename;
use BioUtil::Seq;

$0 = basename($0);
my $usage = <<USAGE;

$0 - locating motif in genomes. 

Motifs could be EITHER plain sequence containing "ACTGN" OR regular
expression like "A[TU]G(?:.{3})+?[TU](?:AG|AA|GA)" for ORFs. 
Degenerate bases like "RYMM.." are also supported by option -d/

Usage: $0 <motif fasta file> <subject fasta file>
Options:

    -d, --degenerate   Motif contains egenerate base
    -h, --help         Show this help information

Attention: In default, motifs are treated as regular expression.
           When option -d given, regular expression may be wrong. 
           For example: "\\w" -> "\\[AT]". In this case you can use "\\.+?"

Example

USAGE

my $args = {};
GetOptions(
    'help|h'       => \$$args{help},
    'degenerate|d' => \$$args{degenerate},
) or die $usage;
die $usage if $$args{help};
die $usage unless @ARGV == 2;

my $queries = read_sequence_from_fasta_file( shift @ARGV );

my $next_seq = FastaReader( shift @ARGV );

print "subject\tquery\tstart\tend\tstrand\tmatched\n";
while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;

    for my $qname ( sort keys %$queries ) {
        my $qseq = $$queries{$qname};

        my $qseq_r = $qseq;
        $qseq_r = degenerate_seq_to_regexp($qseq_r) if $$args{degenerate};

        my $matches = match_regexp( $qseq_r, $seq );
        for my $match (@$matches) {
            my ( $start, $end, $matched ) = @$match;
            $start += 1;
            $end   += 1;
            print "$header\t$qname\t$start\t$end\t+\t$matched\n";
        }

        my $qseq_r = revcom($qseq);
        $qseq_r = degenerate_seq_to_regexp($qseq_r) if $$args{degenerate};
        my $matches = match_regexp( $qseq_r, $seq );
        for my $match (@$matches) {
            my ( $start, $end, $matched ) = @$match;
            $start += 1;
            $end   += 1;
            print "$header\t$qname\t$start\t$end\t-\t"
                . revcom($matched) . "\n";
        }
    }
}

=head2 degenerate_seq_to_regexp

Translate degenerate sequence to regular expression.

=cut

sub degenerate_seq_to_regexp {
    my ($seq) = @_;
    my %bases = (
        'A' => 'A',
        'T' => 'T',
        'U' => 'U',
        'C' => 'C',
        'G' => 'G',
        'R' => '[AG]',
        'Y' => '[CT]',
        'M' => '[AC]',
        'K' => '[GT]',
        'S' => '[CG]',
        'W' => '[AT]',
        'H' => '[ACT]',
        'B' => '[CGT]',
        'V' => '[ACG]',
        'D' => '[AGT]',
        'N' => '[ACGT]',
    );
    return join '', map { exists $bases{$_} ? $bases{$_} : $_ }
        split //, uc $seq;
}

=head2 match_regexp

Find all sites matching the regular expression.

See https://github.com/shenwei356/bio_scripts/blob/master/sequence/fasta_locate_motif.pl

=cut

sub match_regexp {
    my ( $r, $s ) = @_;
    my @matched = ();
    my $pos     = -1;
    while ( $s =~ /($r)/ig ) {
        $pos = pos $s;

        # return start, end, matched string
        # start and end are 0-based
        push @matched, [ $pos - length($1), $pos - 1, $1 ];
        pos $s = $pos - length($1) + 1;
    }
    return \@matched;
}
