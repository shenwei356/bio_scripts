#!/usr/bin/env perl

# Function: Batch compute pI (isoelectric point) and Mw (molecular weight)
#           via submiting sequences to Compute pI/Mw tool at ExPASy.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    : 2013-10-16
# Update  : 2014-07-29

use strict;

# try to use BioUtil::Seq
if ( eval { require BioUtil::Seq; 1; } ne 1 ) {
    die "\nPlease install BioUtil::Seq by CPAN:\n"
        . "  cpan install BioUtil\n\n";
}
else {
    BioUtil::Seq->import();
}

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

while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;

    ( $success, $pi, $mw ) = &$PI($seq);
    unless ($success) {
        print
            "$pi. Please check whether the amino acid sequence contains illegal characters.\r\n"
            ;    # here $pi is the status_line of response
        next;
    }
    print "$header\t$pi\t$mw\r\n";
    print OUT "$header\t$pi\t$mw\r\n";
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
