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
use Getopt::Long;

use BioUtil::Misc;
use BioUtil::Seq;
use BioUtil::Util;

my $usage = sprintf "
Usage: %s [options]

Options:
    -e FILE    Enzymefile (from Rebase)
    -i FILE    Fasta file
    -l FILE    Enzyme list file
    -t INT     Threshold  [%d]

Example:
    
    %s -e embossre.enz -i test.fasta -t 10

See more: https://github.com/shenwei356/bio_scripts
", basename($0), 1 << 30, basename($0);

my $help       = 0;
my $enzymefile = "";
my $seqfile    = "";
my $threshold  = 1 << 30;

GetOptions(
    'help|h' => \$help,
    'e=s'    => \$enzymefile,
    'i=s'    => \$seqfile,
    't=i'    => \$threshold,
) or die $usage;

die $usage if $help;
die $usage unless $enzymefile ne "" and $seqfile ne "";
die "threshold should > 0\n" unless $threshold > 0;

# ===============================================================

my $enzs     = parse_embossre($enzymefile);
my %subenzs  = ();
my %listhash = ();

my $listfile = shift @ARGV;
if ( defined $listfile ) {
    my $list = get_list_from_file($listfile);
    %listhash = map { $_ => 0 } @$list;
    for my $enz ( keys %$enzs ) {
        if ( exists $listhash{$enz} ) {
            $subenzs{$enz} = $$enzs{$enz};
        }
    }
}
else {
    %subenzs = %$enzs;
}

%listhash = ();
%listhash = map { $_ => 0 } keys %subenzs;

# show process
local $| = 1;
my $n = 0;

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
            $listhash{$enz}++;
            if ($listhash{$enz} >= $threshold) {
                delete $subenzs{$enz};
                delete $listhash{$enz};
            }
        }
    }

    # show process
    $n++;
    print STDERR "\rcheck seq $n";
}
$| = 0;

print STDERR "\n";
for ( sort { $listhash{$b} <=> $listhash{$a} } keys %listhash ) {
    print "$_\t$listhash{$_}\n";
}
