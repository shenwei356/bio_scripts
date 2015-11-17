#!/usr/bin/env perl
# Author      : Wei Shen
# Email       : shenwei356@gmail.com
# Date        : 2011-07-20
# Last Update : 2011-07-20
use strict;
use File::Basename;

$0 = basename($0);
die "Usage: $0 gb_file\n" unless @ARGV == 1;
my $file = shift;

my ($definition, $version, $gi, $seq);

open IN, $file or die "failed to open file: $file\n";
$/ = "\n//";
while (<IN>) {
    next unless /DEFINITION  (.+)\./;
    $definition = $1;
    #print "$definition\n";
    next unless /VERSION     (.+)  GI\:(.+)\r?\n/;
    $version = $1;
    $gi      = $2;
    #print "$version, $gi\n";
    $seq     = substr($_, index($_, 'ORIGIN') + 6);
    $seq     =~ s/\/\/.*//s;
    $seq     =~ s/\s+//g;
    $seq     =~ s/\d+//g;
    #print "$seq\n";
    #print length($seq),"\n";
    print ">gi|$gi|gb|$version| $definition\n".(format_seq($seq, 60))."\n";
}
$/ = "\n";
close IN;


sub format_seq($$){
    my ($s, $n) = @_;
    my $s2 ='';
    my ($j, $int);
    $int = int ((length $s) / $n);
    for($j = 0 ; $j <= $int - 1; $j ++){
        $s2 .= substr($s, $j * $n, $n)."\n";
    }
    $s2 .= substr($s, $int * $n);
    return $s2;
}
