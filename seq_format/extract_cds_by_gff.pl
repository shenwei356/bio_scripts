#!/usr/bin/env perl

use strict;
use File::Basename;
use Getopt::Long;
use BioUtil::Seq;
use BioUtil::Util;
use Data::Dumper;

$0 = basename($0);
my $usage = qq(
Usage: $0 [options] gff_file fasta_file
Options:
    -t,   --type         gene type (CDS or mRNA) [CDS]
    -us,  --up-stream    up stream length [0]
    -ds,  --down-stream  down stream length [0]
    -h,   --help         show this usage

);

my $argv = {};
$$argv{type}        = 'CDS';
$$argv{up_stream}   = 0;
$$argv{down_stream} = 0;

GetOptions(
    'help|h'           => \$$argv{help},
    'type|t=s'         => \$$argv{type},
    'up-stream|us=s'   => \$$argv{up_stream},
    'down-stream|ds=s' => \$$argv{down_stream},
);

die $usage if $$argv{help};
die $usage if scalar(@ARGV) != 2;

check_positive_integer( $$argv{up_stream} + 1 );
check_positive_integer( $$argv{down_stream} + 1 );

my ( $gff_file, $fasta_file ) = @ARGV;

my $genes = read_gff_file($gff_file);

# print Dumper($genes);

my $next_seq = FastaReader($fasta_file);
while ( my $fa = &$next_seq() ) {
    my ( $name, $genome ) = @$fa;
    next if not exists $$genes{$name};

    for my $gene ( @{ $$genes{$name} } ) {
        next if lc $$gene{type} ne lc $$argv{type};    # specific type
        my $seq = '';

        if ( $$gene{strand} eq '+' ) {
            my $s = $$gene{start} - $$argv{up_stream} - 1;
            $s = 0 if $s < 0;
            $seq = substr(
                $genome, $s,
                $$gene{end}
                    - $$gene{start}
                    + $$argv{down_stream} + 1

            );
        }
        else {
            my $s = $$gene{start} - $$argv{down_stream} - 1;
            $s = 0 if $s < 0;
            $seq = revcom(
                substr(
                    $genome, $s,
                    $$gene{end} - $$gene{start} + $$argv{up_stream} + 1
                )
            );
        }
        printf( ">%s_%d..%d..%s\n%s",
            $name, $$gene{start}, $$gene{end}, $$gene{strand},
            format_seq($seq) );
    }

}

sub read_gff_file {
    my ($file) = @_;
    my $genes = {};
    open( my $fh, "<", $file ) or die "fail to open file: $file\n";
    while (<$fh>) {
        my @data = split( /\s+/, $_ );
        next unless scalar(@data) >= 9;
        my $name = $data[0];
        my $gene = {};
        ( $$gene{type}, $$gene{start}, $$gene{end}, $$gene{strand} )
            = ( $data[2], $data[3], $data[4], $data[6] );
        if ( not exists $$genes{$name} ) {
            $$genes{$name} = [];
        }
        push @{ $$genes{$name} }, $gene;

    }
    close($fh);
    return $genes;
}

