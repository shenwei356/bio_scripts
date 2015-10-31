#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts

use strict;
use File::Basename;
use Getopt::Long;
use BioUtil::Seq;
use BioUtil::Util;

$0 = basename($0);
my $usage = <<USAGE;

Remove duplicated fasta names

Usage: $0 [options] [fastafiles...]
Options:

   -l   Output line length. should be >= 0, 0 for no formating [70]
   -h   Show this help information.
https://github.com/shenwei356/bio_scripts

USAGE

my $help       = 0;
my $linelength = 70;
GetOptions(
    'help|h' => \$help,
    'l=i'    => \$linelength,
) or die $usage;

die $usage if $help;
if ( $linelength < 0 ) {
    die sprintf "value of -l (%d) should be greatter or equal to 0\n",
      $linelength;
}

# get the file list
my @files = file_list_from_argv(@ARGV);

my $names = {};
for my $file (@files) {
    my $next_seq = FastaReader($file);
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;
        if ( exists $$names{$header} ) {
            $$names{$header}++;
            $header = "$header r$$names{$header}";
        }
        else {
            $$names{$header} = 1;
        }

        if ( $linelength > 0 ) {
            print ">$header\n", format_seq( $seq, $linelength );
        }
        else {
            print ">$header\n", $seq, "\n";
        }
    }
}
