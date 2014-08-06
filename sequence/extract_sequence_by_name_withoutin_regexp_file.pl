#!/usr/bin/perl

usage() unless @ARGV == 2;
$regex_file = shift;
$seq_file   = shift;

open IN, $regex_file or die "fail to open file $regex_file!\n";
while ( <IN> ) {
    s/\r?\n//g;
    next if /^\s*$/;
    s/^\s+|\s+$//;
    $data{quotemeta $_} = 1;
}
close IN;
@keys = keys %data;

open IN, $seq_file or die "fail to open sequence file $seq_file!\n";
$/ = '>';<IN>;
while ( <IN> ) {
    s/>//;
    ( $head, $seq ) = split "\r?\n", $_, 2;
    @in = grep { $_ if $head =~ /$_/i } @keys;
    # print "@in\n" if @in > 0;
    if ( @in == 0){
        print ">$head\n$seq";
    }
}
$/ = "\n";
close IN;

sub usage{
    die qq(
Usage: $0 <regular expression file> <sequence_file> 
    
);
}
