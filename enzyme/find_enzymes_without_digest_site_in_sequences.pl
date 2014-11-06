#!/usr/bin/env perl
#
# Copyright 2014 Wei Shen (shenwei356#gmail.com). All rights reserved.
# Use of this source code is governed by a MIT-license
# that can be found in the LICENSE file.
# https://github.com/shenwei356/bio_scripts/

# embossre.enz
#   ftp://ftp.neb.com/pub/rebase/

use strict;
use File::Basename;
use BioUtil::Misc;
use BioUtil::Seq;
use BioUtil::Util;

my $usage = sprintf "
Usage: %s <embossre.enz> <fasta file> [enzyme list file]

", basename($0);
die $usage unless @ARGV == 3 or @ARGV == 2;

my $enzymefile = shift @ARGV;
my $seqfile    = shift @ARGV;

my $enzs    = parse_embossre($enzymefile);
my %subenzs = ();

my $listfile = shift @ARGV;
if ( defined $listfile ) {
    my $list = get_list_from_file($listfile);
    my %listhash = map { $_ => 0 } @$list;
    for my $enz ( keys %$enzs ) {
        if ( exists $listhash{$enz} ) {
            $subenzs{$enz} = $$enzs{$enz};
        }
    }
}
else {
    %subenzs = %$enzs;
}

# show process
local $| = 1;
my $n    = 0;
my $sum  = scalar keys %subenzs;
my $left = $sum;

my $next_seq = FastaReader($seqfile);
while ( my $fa = &$next_seq() ) {
    my ( $header, $seq ) = @$fa;
    $seq = uc $seq;
    my $revcom = revcom($seq);

    for my $enz ( keys %subenzs ) {
        my $e       = $subenzs{$enz};
        my $pattern = $$e{pattern_regexp};
        # check enzyme digest site
        if ( $seq =~ /$pattern/ or $revcom =~ /$pattern/ ) {
            delete $subenzs{$enz};
        }
    }

    # show process
    $n++;
    $left = scalar keys %subenzs;
    print STDERR "\rcheck seq $n, candidate: $left / $sum";
}
$| = 0;

print STDERR "\n";
print "$_\n" for sort keys %subenzs;
