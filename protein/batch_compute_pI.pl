#!/usr/bin/env perl

# Function: Batch compute pI (isoelectric point) and Mw (molecular weight)
#           via submiting sequences to Compute pI/Mw tool at ExPASy.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    : 2013-10-16
# Update  : 2014-07-29

use strict;

my $usage = <<"USAGE";

Function: Batch compute pI (isoelectric point) and Mw (molecular weight)
          via submiting sequences to Compute pI/Mw tool at ExPASy
 Contact: Wei Shen <shenwei356#gmail.com>
   Usage: $0 amino_acid_fasta_file
    
USAGE
die $usage
    unless @ARGV == 1;
my $aa_file = shift @ARGV;

# initialize fasta file parser
my $next_seq = FastaReader($aa_file);

# initialize pI request
my $PI = &compute_pi();

my ( $head, $seq );
my ( $success, $pi, $mw );
my $out_file = "$aa_file.result.txt";
open OUT, ">", $out_file
    or die "fail to write file $out_file\n";

while (1) {
    ( $head, $seq ) = &$next_seq();
    last if $head eq "" and $seq eq "";

    # print ">$head\n$seq\n";

    ( $success, $pi, $mw ) = &$PI($seq);
    unless ($success) {
        print
            "$pi. Please check whether the amino acid sequence contains illegal characters.\r\n"
            ;    # here $pi is the status_line of response
        next;
    }
    print "$head\t$pi\t$mw\r\n";
    print OUT "$head\t$pi\t$mw\r\n";
}

close OUT;

# Compute pI/Mw via submiting sequence to Compute pI/Mw tool at ExPASy.
#
# See more: http://web.expasy.org/compute_pi/
#
# Example:
#
#    my @proteins = qw/AYYAYYAYAYAY ACACAGACG ---/;
#    my $PI = &compute_pi();
#    my ( $success, $pi, $mw );
#    for my $protein (@proteins) {
#        ( $success, $pi, $mw ) = &$PI($protein, "average");
#        # ( $success, $pi, $mw ) = &$PI($protein, "monoisotopic");
#        unless ($success) {
#            print "$pi\n";    # here $pi is the status_line of response
#            next;
#        }
#        print "($pi, $mw)\n";
#    }
sub compute_pi() {
    use LWP::UserAgent;

    my $ua  = LWP::UserAgent->new;
    my $url = "http://web.expasy.org/cgi-bin/compute_pi/pi_tool";
    my ( $res, $formdata, $result );

    return sub($$) {
        my ( $protein, $resolution ) = @_;
        $resolution = "average" unless defined $resolution; # or  monoisotopic
        $formdata = [
            protein    => $protein,
            resolution => $resolution,
            file       => ""
        ];

        $res = $ua->post( $url, $formdata );

        # 0 means failed
        return ( 0, $res->status_line )
            unless $res->is_success;

        $result = $res->content;
        $result =~ /Theoretical pI\/Mw: ([\d\.]+)\s\/\s([\d\.]+)/;

        # 1 means success
        return ( 1, $1, $2 );
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
sub FastaReader($) {
    my ($file) = @_;
    open IN, "<", $file
        or die "Fail to open file: $file!\n";
    local $/ = '>';
    <IN>;
    $/ = '\n';

    my ( $line, $head, $seq );
    return sub() {
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
    };
}
