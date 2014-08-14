#!/usr/bin/perl

usage() unless @ARGV == 2;
$regex_file = shift;
$seq_file   = shift;

open IN, $regex_file or die "fail to open file $regex_file!\n";
while (<IN>) {
    s/\r?\n//g;
    next if /^\s*$/;
    s/^\s+|\s+$//;
    $data{ quotemeta $_ } = 1;
}
close IN;
@keys = keys %data;

my $next_seq = FastaReader($seq_file);
while ( my $fa = &$next_seq() ) {
    my ( $head, $seq ) = @$fa;

    for (@keys) {
        if ( $head =~ /$_/ ) {
            print ">$head\n$seq\n";
            last;
        }
    }

}

sub usage {
    die qq(
Usage: $0 <regular expression file> <sequence_file> 
    
);
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
