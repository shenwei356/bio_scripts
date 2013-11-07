#!/usr/bin/perl

usage() unless @ARGV == 2;
$regex_file = shift;
$seq_file   = shift;

open IN, $regex_file or die "fail to open file $regex_file!\n";
while ( <IN> ) {
    s/\r?\n//g;
    s/\s+$//;
    $data{$_} = 1;
}
close IN;


open IN, $seq_file or die "fail to open sequence file $seq_file!\n";
$/ = '>';<IN>;
while ( <IN> ) {
    s/>//;
    ( $head, $seq ) = split "\r?\n", $_, 2;
    if ( defined $data{$head} ){
        print ">$head\n$seq";
    }

}
close IN;

sub usage{
    die qq(
Usage: $0 <name_file> <sequence_file> 
    
);
}
    