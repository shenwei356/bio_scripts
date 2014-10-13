#!/usr/bin/env perl
# http://doc.bioperl.org/bioperl-live/Bio/SeqIO/fastq.html
# http://doc.bioperl.org/bioperl-live/Bio/Seq/Quality.html

# make sure the reads in the two fastq files has same order!
use strict;

use Bio::SeqIO;
use Parallel::Runner;

die "usage: $0 <read.1.fq> <read.2.fq>\n"
    unless @ARGV == 2;

my $fqfile1 = shift @ARGV;
my $fqfile2 = shift @ARGV;

# ===========================================================

print "read $fqfile1\n";
my $headers1 = get_headers($fqfile1);

print "read $fqfile2\n";
my $headers2 = get_headers($fqfile2);

# ===========================================================

print "find common IDs: ";
my $headers = {};
for my $header ( keys %$headers1 ) {
    next unless exists $$headers2{$header};
    $$headers{$header} = 1;
}
my $n = keys %$headers;
print "$n\n";

die "sadly, no paired reads found\n" if $n == 0;

# ===========================================================

my $runner = Parallel::Runner->new(2);

print "extract $fqfile1\n";
$runner->run( sub { extract( $headers, $fqfile1 ); } );

print "extract $fqfile2\n";
$runner->run( sub { extract( $headers, $fqfile2 ); } );

$runner->finish;

# ===========================================================

sub extract {
    my ( $headers, $fqfile ) = @_;

    my $in = Bio::SeqIO->new(
        -format => 'fastq',
        -file   => $fqfile
    );

    my $fqfileout = $fqfile;
    $fqfileout =~ s/\.(fq|fastq)$//i;
    $fqfileout .= ".pe.fq";
    my $out = Bio::SeqIO->new(
        -format => 'fastq',
        -file   => ">$fqfileout"
    );
    while ( my $seq = $in->next_seq ) {
        next unless exists $$headers{ $seq->id };
        $out->write_seq($seq);
    }
    $out->close();
    $in->close();
}

sub get_headers {
    my ($fqfile) = @_;

    my $headers = {};

    my $in = Bio::SeqIO->new(
        -format => 'fastq',
        -file   => $fqfile
    );

    while ( my $seq = $in->next_seq ) {
        $$headers{ $seq->id } = "";
    }

    $in->close();
    return $headers;
}
