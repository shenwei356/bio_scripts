#!/usr/bin/env perl
use strict;
use File::Basename;
use Data::Dumper;

$0 = basename($0);
my $usage = <<USAGE;
Usage: $0 <tsv1> <index1> <tsv2> <index2>

USAGE
die($usage) unless scalar(@ARGV) == 4;
my ( $tsv1, $index1, $tsv2, $index2 ) = @ARGV;

sub tsv2map ($$) {
    my ( $file, $index ) = @_;
    $index = 1 unless defined($index);    # index column, defautl: 1

    my $data = {};    # data is a hash reference, I prefer this.
    open( my $fh, "<", $file ) or die("failed to open file: $file\n");
    while (<$fh>) {
        chomp($_);
        my @items = split( /\t/, $_ );
        if ( scalar(@items) < $index ) {    # verify $index
            die "number of column in file ($file) < index ($index).\n";
        }
        my $key = $items[ $index - 1 ];     # get the key
        $$data{$key} = $_;                  # store key => value
    }
    close $fh;

    return $data;
}

my $data_tsv2 = tsv2map( $tsv2, $index2 );

# print Dumper($data_tsv2);
# result:
# $VAR1 = {
#           '123' => '123 onetwothree',
#           'str' => 'str string',
#           '245' => '245 twofourfive'
#         };

# parse tsv1
open( my $fh, "<", $tsv1 ) or die("failed to open file: $tsv1\n");
while (<$fh>) {
    chomp($_);
    my @items = split( /\t/, $_ );
    if ( scalar(@items) < $index1 ) {
        die "number of column in file ($tsv1) < index ($index1).\n";
    }
    my $key = $items[ $index1 - 1 ];    # get the key

    if ( exists $$data_tsv2{$key} ) {   # check if key existed in tsv2
        print "$_\t$$data_tsv2{$key}\n";
    }
    else {
        print "$_\n";
    }
}
close $fh;
