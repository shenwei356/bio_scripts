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

# FastaReader is a fasta file parser using closure.
# FastaReader returns an anonymous subroutine, when called, it
# will return a fasta record which is reference of an array
# containing fasta header and sequence.
#
# A boolean argument is optional. If set as "true", "return" ("\r") and
# "new line" ("\n") symbols in sequence will not be trimed.
#
# Example:
#
#    # my $next_seq = FastaReader("test.fa", 1);
#    my $next_seq = FastaReader("test.fa");
#
#    while ( my $fa = &$next_seq() ) {
#        my ( $header, $seq ) = @$fa;
#
#        print ">$header\n$seq\n";
#    }
#
sub FastaReader {
    my ( $file, $not_trim ) = @_;

    my ( $last_header, $seq_buffer ) = ( '', '' ); # buffer for header and seq
    my ( $header,      $seq )        = ( '', '' ); # current header and seq
    my $finished = 0;

    open FH, "<", $file
        or die "fail to open file: $file!\n";

    return sub {

        if ($finished) {                           # end of file
            return undef;
        }

        while (<FH>) {
            s/^\s+//;    # remove the space at the front of line

            if (/^>(.*)/) {    # header line
                ( $header, $last_header ) = ( $last_header, $1 );
                ( $seq,    $seq_buffer )  = ( $seq_buffer,  '' );

                # only output fasta records with non-blank header
                if ( $header ne '' ) {
                    $seq =~ s/\s+//g unless $not_trim;
                    return [ $header, $seq ];
                }
            }
            else {
                $seq_buffer .= $_;    # append seq
            }
        }
        close FH;
        $finished = 1;

        # last record
        # only output fasta records with non-blank header
        if ( $last_header ne '' ) {
            $seq_buffer =~ s/\s+//g unless $not_trim;
            return [ $last_header, $seq_buffer ];
        }
    };
}
