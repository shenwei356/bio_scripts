#!/usr/bin/perl

usage() unless @ARGV == 2;
$regex    = shift;
$seq_file = shift;

open IN, $seq_file or die "fail to open sequence file $seq_file!\n";
$/ = '>';
<IN>;
while (<IN>) {
    s/>//;
    ( $head, $seq ) = split "\r?\n", $_, 2;
    if ( $head =~ /$regex/ ) {
        print ">$head\n$seq";
    }
}
$/ = "\n";
close IN;

sub usage {
    die qq(
Usage: $0 <regular expression> <sequence_file> 
    
);
}
