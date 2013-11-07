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
@keys = keys %data;

open IN, $seq_file or die "fail to open sequence file $seq_file!\n";
$/ = '>';<IN>;
while ( <IN> ) {
    s/>//;
    ( $head, $seq ) = split "\r?\n", $_, 2;    
    for (@keys){
        
        if ( $head =~ /$_/ ){
            print ">$head\n$seq";
            last;
        }
    }

}
close IN;

sub usage{
    die qq(
Usage: $0 <regular expression file> <sequence_file> 
    
);
}

sub format_for_regexp($){
    my ($s) = @_;
    my @t = qw(\( \) { } [ ] . * + ? - /);
    my (@t1, @t2);
    push @t1, "\\"."$_" for @t;
    push @t2, "\\\\"."$_" for @t1;
    print "$_\t" for @t1;print "\n";
    print "$_\t" for @t2;print "\n";
    print "$s\n";
    for (0..$#t1) {
        print "$s           $t1[$_], $t2[$_]\n";
        $s =~ s/$t1[$_]/$t2[$_]/g;
    }
}