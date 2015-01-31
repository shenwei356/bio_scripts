#!/usr/bin/env perl
# Fuction     : Run EMBOSS restrict, and parse restriction fragments
# Author      : Wei Shen
# Email       : shenwei356@gmail.com, http://shenwei.me
# Date        : 2014-12-18
# Last Update : 2015-01-31

use strict;
use File::Basename;
use Getopt::Long;
use Parallel::Runner;

use BioUtil::Misc;
use BioUtil::Seq;
use BioUtil::Util;

local $| = 1;

$0 = basename($0);
my $usage = <<USAGE;
restrict_fragments.pl

    Parsing restriction fragments for chosing a appropriate 
    enzyme to identify multi genomes

Usage: $0 [options]

Options:
    -e FILE    Enzymefile (from Rebase)
    -i FILE    Fasta file
    -l FILE    Enzyme list file
    -t INT     Thread number
    -linear    Linear genome

Example:

See more: https://github.com/shenwei356/bio_scripts
USAGE

my $help          = 0;
my $enzymefile    = "";
my $seqfile       = "";
my $listfile      = "";
my $threads       = 4;
my $linear_genome = 0;

GetOptions(
    'help|h' => \$help,
    'e=s'    => \$enzymefile,
    'i=s'    => \$seqfile,
    'l=s'    => \$listfile,
    't=i'    => \$threads,
    'linear' => \$linear_genome,
) or die $usage;

die $usage if $help;
die $usage unless $enzymefile ne "" and $seqfile ne "";

#=====================[ run restrict ]====================

my $enzs     = parse_embossre($enzymefile);
my %subenzs  = ();
my %listhash = ();

my $dir = '';
if ( $listfile ne "" ) {
    my $list = get_column_data( $listfile, 1 );
    %listhash = map { $_ => 0 } @$list;
    for my $enz ( keys %$enzs ) {
        if ( exists $listhash{$enz} ) {
            $subenzs{$enz} = $$enzs{$enz};
        }
    }
    $dir = "re.$seqfile.digestedby.$listfile";
}
else {
    %subenzs = %$enzs;
    $dir     = "re.$seqfile.digestedby.$enzymefile";
}

%listhash = ();
%listhash = map { $_ => 0 } keys %subenzs;

# unless ( -e $dir and -d $dir ) {
rm_and_mkdir($dir);
my $runner = Parallel::Runner->new($threads);
for my $enz ( keys %subenzs ) {
    $runner->run(
        sub {
            run_emboss_restrict( $dir, $enz );
        }
    );
}
$runner->finish;

# }

sub run_emboss_restrict {
    my ( $dir, $enzyme ) = @_;
    my $resultfile = "$dir/$seqfile.$enzyme.re";
    return if -e $resultfile;
    print STDERR "$enzyme\n";
    my $cmd = "restrict -auto -solofragment -limit "
        . "-sequence $seqfile -outfile $resultfile -enzymes $enzyme ";
    $cmd .= " -plasmid " unless $linear_genome;

    my $fail = run($cmd);
    die "failed to run:$cmd\n" if $fail;
}

# ===========[ Parsing restriction fragments ]=============

my @files = glob "$dir/*.re";

my $fragments = {};
my $stats     = {};
for my $file (@files) {
    my ( $enzyme, $seq ) = (undef) x 2;

    open my $fh, $file
        or die "fail to read enzyme file $file\n";
    while (<$fh>) {
        if (/^#\s+\-enzymes (.+)/) {    # enzyme name
            $enzyme = $1;
        }
        elsif (/^# Sequence: (.+)\s+from/) {    # sequence name
            $seq = $1;
        }
        elsif (/^# \t([\d\t]+)$/) {             # fragment size
            if ( ref $$fragments{$enzyme}{$seq} ne ref [] ) {
                $$fragments{$enzyme}{$seq} = [];
            }
            push @{ $$fragments{$enzyme}{$seq} }, split( /\t/, $1 );
        }
    }
    close $fh;

    my $n = 0;
    for my $seq ( keys %{ $$fragments{$enzyme} } ) {
        my @frags = sort { $b <=> $a } @{ $$fragments{$enzyme}{$seq} };

        # print "$enzyme\n  $seq\n  @frags\n";
        $$fragments{$enzyme}{$seq} = \@frags;
        $n += scalar @frags;
    }
    $$stats{$enzyme}{nfrags} = $n;
}

# ===========[ Output restriction fragments ]=============

my $outfile = "$seqfile.digestedby.$listfile.frag";

open OUT, ">", $outfile
    or die "fail to write file $outfile\n";

my $frag = {};
for my $enzyme (
    sort { $$stats{$a}{nfrags} <=> $$stats{$b}{nfrags} }
    keys %$fragments
    )
{

    print OUT "-" x 79, "\n", "$enzyme\n";
    for my $seq ( sort keys %{ $$fragments{$enzyme} } ) {
        my @frags = @{ $$fragments{$enzyme}{$seq} };
        print OUT "$seq: @frags\n";
    }
}

close OUT;
