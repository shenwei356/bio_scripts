#!/usr/bin/env perl

use strict;

# try to use BioUtil::Seq
if ( eval { require BioUtil::Seq; require BioUtil::Util; 1; } ne 1 ) {
    die "\nPlease install BioUtil::Seq by CPAN:\n"
        . "  cpan install BioUtil\n\n";
}
else {
    BioUtil::Seq->import();
    BioUtil::Util->import();
}

usage() unless @ARGV == 2;
my $regex    = shift;
my $seq_file = shift;

my $next_seq = FastaReader($seq_file);
while ( my $fa = &$next_seq() ) {
    my ( $head, $seq ) = @$fa;

    if ( $seq =~ /$regex/i ){
        print ">$head\n$seq\n";
    }
}

sub usage{
    die qq(
Usage: $0 <regular expression> <sequence_file> 
    
);
}