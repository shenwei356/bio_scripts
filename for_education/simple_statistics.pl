#!/usr/bin/env perl

# Copyright 2014 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
# https://github.com/shenwei356

use strict;

my $usage = <<USAGE;

usage: simple_statistics <infile> <column>

<infile> is a plain text file. each column should be seperated by TAB(\\t)
<column> is the column number of the table.

USAGE

die $usage unless @ARGV == 2;

my $file   = shift @ARGV;
my $column = shift @ARGV;

my $data = get_column_data( $file, $column );

printf "#.\t%d\n", scalar @$data;
printf "min.\t%d\n", min($data);
printf "max.\t%d\n", max($data);

my ($mean, $stdev) = mean_and_stdev($data);

printf "mean.\t%.2f\n", $mean;
printf "stdev.\t%.2f\n", $stdev;


sub get_column_data {
    my ( $file, $column ) = @_;
    unless ( $column =~ /^(\d+)$/ and $column > 0 ) {
        warn
            "column number ($column) should be an integer and greater than 0.\n";
        $column = 1;
    }

    open IN, "<", $file or die "failed to open file: $file\n";
    my @linedata = ();
    my @data     = ();
    my $n        = 0;
    while (<IN>) {
        s/\r?\n//;
        @linedata = split /\t/, $_;
        $n = scalar @linedata;
        next unless $n > 0;
        
        if ( $column > $n ) {
            die
                "number of columns of this line ($n) is less than given column number ($column)\n";
        }

        push @data, $linedata[ $column - 1 ];
    }
    close IN;

    return \@data;
}

# you can also modules
# use List::Util qw/max min sum/;

sub max {
    my ($list) = @_;
    my $max = shift @$list;
    for (@$list) {
        $max = $_ if $_ > $max;
    }
    return $max;
}

sub min {
    my ($list) = @_;
    my $min = shift @$list;
    for (@$list) {
        $min = $_ if $_ < $min;
    }
    return $min;
}

sub mean_and_stdev($) {
    my ($list) = @_;
    return ( 0, 0 ) if @$list == 0;
    my $sum = 0;
    $sum += $_ for @$list;
    my $sum_square = 0;
    $sum_square += $_ * $_ for @$list;
    my $mean     = $sum / @$list;
    my $variance = $sum_square / @$list - $mean * $mean;
    my $std      = sqrt $variance;
    return ( $mean, $std );
}
