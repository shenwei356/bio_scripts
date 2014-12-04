#!/usr/bin/env perl
# https://github.com/shenwei356/bio_scripts
use strict;
use Getopt::Long;
use File::Temp qw/ tempfile/;
use BioUtil::Seq;
use BioUtil::Util;

local $| = 1;

my @GAPS = ( '-', '.' );
my $tmpfile_prefix = "fasta_trim_aligned_fasta_tmpfile_";

my $usage = <<USAGE;

Remove common gaps in aligned fasta sequences.

Reading from STDIN is supported. But in this case, sequences will be saved to
disk temporarily. Because this script reads sequences twice to reduce memery 
usage.

Usage: $0 [options] [aligned fasta file...]
Options:
    -h,  --help                Show this help information
    -g,  --gaps                Gap symbols [-.]
    -l,  --linelength          Line length

https://github.com/shenwei356/bio_scripts

USAGE

my $para = {};
$$para{linelength} = 70;

GetOptions(
    'help|h'         => \$$para{help},
    'gaps|g=s'       => \$$para{gaps},
    'linelength|l=i' => \$$para{linelength},
) or die $usage;

die $usage if $$para{help};

# gap symbols
my %GAPSMAP = ();
if ( $$para{gaps} ) {
    @GAPS = split //, $$para{gaps};
}
$GAPSMAP{$_} = 1 for @GAPS;

my $use_stdin = 0;
my ( $tmp_file_fh, $tmp_file ) = (undef) x 2;

my @files = ();
for my $file (@ARGV) {
    for my $f ( glob $file ) {
        push @files, $f;
    }
}
if ( @files == 0 ) {
    push @files, 'STDIN';
    ( $tmp_file_fh, $tmp_file )
        = tempfile( $tmpfile_prefix . "XXXXXX", DIR => ".", SUFFIX => '.fa' );

    $use_stdin = 1;
}

print STDERR "sequences from STDIN is saved in $tmp_file\n" if $use_stdin;
print STDERR "check...\n";

my $gaploc  = {};    # store the gap location
my $do_once = 1;
my ( $header, $seq, $len, $i, $base ) = (undef) x 5;
my ( $sum, $n ) = (0) x 2;
for my $file (@files) {
    my $next_seq = FastaReader($file);
    while ( my $fa = &$next_seq() ) {
        ( $header, $seq ) = @$fa;
        $sum++;
        print STDERR "\rsum: $sum";
        if ($do_once) {
            $len = length $seq;
            $$gaploc{$_} = 1 for 0 .. ( $len - 1 );
            $do_once = 0;
        }

        for $i ( 0 .. ( $len - 1 ) ) {
            $base = substr $seq, $i, 1;
            if ( $GAPSMAP{$base} != 1 ) {    # it's not a gap!
                delete $$gaploc{$i};
            }
        }

        if ( scalar keys %$gaploc == 0 ) {
            close $tmp_file_fh if $use_stdin;
            remove_tmpfile()   if $use_stdin;
            die "no gap to trim\n";
        }

        print $tmp_file_fh ">$header\n$seq\n" if $use_stdin;
    }
}

close $tmp_file_fh if $use_stdin;

print STDERR scalar keys %$gaploc, " gaps to trim\n";
print STDERR "\nextract sequences...\n";

@files = ($tmp_file) if $use_stdin;

my @index = keys %$gaploc;
for my $file (@files) {
    my $next_seq = FastaReader($file);
    while ( my $fa = &$next_seq() ) {
        ( $header, $seq ) = @$fa;
        $n++;
        print STDERR "\r$n / $sum";
        print ">$header\n",
            format_seq( delete_string_elements_by_indexes( \$seq, \@index ),
            $$para{linelength} );
    }
}

print STDERR "\n";

remove_tmpfile() if $use_stdin;

sub remove_tmpfile {
    print STDERR "\nremove temporary files\n";
    for ( glob "$tmpfile_prefix*" ) {
        unlink $_ or die "fail to remove $_\n";
    }
}
