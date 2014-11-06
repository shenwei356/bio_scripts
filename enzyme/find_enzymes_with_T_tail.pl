#!/usr/bin/env perl
#
# Copyright 2014 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
# https://github.com/shenwei356/bio_scripts/

use strict;
use BioUtil::Misc;

die "usage: $0 embossre.enz\n"
    unless @ARGV == 1;

my $file = shift @ARGV;
my $d = shift @ARGV;

my $enzs = parse_embossre($file);

for my $enz (sort keys %$enzs) {
    my $e = $$enzs{$enz};
    next unless $$e{cuts_number} == 2
        and $$e{c1} - $$e{c2} == 1
        and substr ($$e{pattern}, $$e{c1} - 1, 1) =~ /[aN]/i;
    print "$enz\n";
}

# there's no enzyme meeting this condition
