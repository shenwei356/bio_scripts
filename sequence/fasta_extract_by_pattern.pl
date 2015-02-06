#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;

use Getopt::Long;
use File::Basename;
use BioUtil::Seq;
use BioUtil::Util;

$0 = basename($0);
my $usage = <<USAGE;

Extract fasta sequences by header (list file) or regular expression (list file)

Version: 2015.02.06
Usage: $0 [options] [fastafiles...]
Options:
    
    -p,  --pattern STRING      Search pattern
    -pf, --patternfile FILE    Pattern list file (use first column)
    -r,  --useregexp           Use regular expression, case ignored
    -d,  --speedup             Delete matched pattern, if you know what it means
    -n,  --not                 Invert match, extract sequences NOT match the pattern
    -s,  --byseq               Match by sequence 
    -h,  --help                Show this help information

Examples:

    1) sequences WITH "bacteria" in header
        $0 -r -p Bacteria *.fa > result.fa
    2) sequences WITHOUT "bacteria" in header
        $0 -r -n -p Bacteria seq1.fa seq2.fa > result.fa
    3) sequences with TTSAA (AgsI digest site) in SEQUENCE. 
       Base S stands for C or G.
        $0 -r -s -p 'TT[C|G]AA' seq.fa > result.fa
    4) sequences (read from STDIN ) with header that matches any patterns
       in list file
        zcat seq.fa.gz | $0 -pf name_list.txt > result.fa

https://github.com/shenwei356/bio_scripts

USAGE

my $para = {};
GetOptions(
    'help|h'           => \$$para{help},
    'useregexp|r'      => \$$para{useregexp},
    'speedup|d'        => \$$para{speedup},
    'not|n'            => \$$para{not},
    'pattern|p=s'      => \$$para{pattern},
    'patternfile|pf=s' => \$$para{patternfile},
    'byseq|s'          => \$$para{byseq},
) or die $usage;
die $usage if $$para{help};

# get patterns
my $patterns = {};
$$patterns{$$para{pattern}} = 1 if $$para{pattern};
if ( $$para{patternfile} ){
    $$patterns{$_} = 1 for @{ get_column_data( $$para{patternfile}, 1 ) };
}
die "no patterns given. Type \"$0 -h\" for help.\n" if keys %$patterns == 0;

# get the file list
my @files = file_list_from_argv(@ARGV);

my $not_trim = 1;
$not_trim = 0 if $$para{byseq};

my ( $sum, $n ) = ( 0, 0 );

for my $file (@files) {

    my $next_seq = FastaReader( $file, $not_trim );
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;
        $sum++;

        # matching object, by header or sequence
        my $object = $header;
        if ( $$para{byseq} ) {
            $object = $seq;
        }

        my $hit = undef;
        if ( $$para{useregexp} ) {    # use regular expression
            for my $p (keys %$patterns) {
                if ( $object =~ /$p/i ) {
                    $hit = 1;
                    delete $$patterns{$p} if $$para{speedup};
                    last;
                }
            }
        }
        else {                        # compare with full header | sequence
            if ( exists $$patterns{$object} ) {
                $hit = 1;
            }
        }

        if ( $$para{not} ) {          # NOT
            next if $hit;
        }
        else {
            next unless $hit;
        }

        $n++;
        if ( $$para{byseq} ) {
            print ">$header\n", format_seq($seq);
        }
        else {
            print ">$header\n$seq";
        }
    }
}

print STDERR "\rHits: $n / $sum\n";
