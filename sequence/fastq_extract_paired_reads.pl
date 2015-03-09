#!/usr/bin/env perl

# make sure the reads in the two fastq files has same order!
use strict;
use Parallel::Runner;
use File::Basename;

$0 = basename($0);
die "usage: $0 <read.1.fq> <read.2.fq>\n"
    unless @ARGV == 2;

my $fqfile1 = shift @ARGV;
my $fqfile2 = shift @ARGV;

# ===========================================================

print "read $fqfile1\n";
my $headers1 = get_headers($fqfile1);

print "read $fqfile2\n";
my $headers2 = get_headers($fqfile2);

# ===========================================================

print "find common IDs: ";
my $headers = {};
for my $header ( keys %$headers1 ) {
    next unless exists $$headers2{$header};
    $$headers{$header} = 1;
}
my $n = keys %$headers;
print "$n\n";

die "sadly, no paired reads found\n" if $n == 0;

# ===========================================================

my $runner = Parallel::Runner->new(2);

print "extract $fqfile1\n";
$runner->run( sub { extract( $headers, $fqfile1 ); } );

print "extract $fqfile2\n";
$runner->run( sub { extract( $headers, $fqfile2 ); } );

$runner->finish;

# ===========================================================

sub extract {
    my ( $headers, $fqfile ) = @_;

    my $fqfileout = $fqfile;
    $fqfileout =~ s/\.(fq|fastq)$//i;
    $fqfileout .= ".pe.fq";
    open my $fh, ">", $fqfileout or die "fail to wrtie file: $fqfileout\n";

    my $next_seq = FastqReader($fqfile);
    my $id = '';
    while ( my $fq = &$next_seq() ) {
        my ( $head, $seq, $qual ) = @$fq;
        $id = (split / /, $head )[0];
        if ($id =~ /(.+)\/\d$/){
            $id = $1;
        }
        next unless exists $$headers{ $id };
        print $fh "\@$head\n$seq\n+\n$qual\n";
    }
}

sub get_headers {
    my ($fqfile) = @_;
    my $headers = {};

    my $next_seq = FastqReader($fqfile);
    my $id = '';
    while ( my $fq = &$next_seq() ) {
        my ( $head, $seq, $qual ) = @$fq;
        $id = (split / /, $head )[0];
        if ($id =~ /(.+)\/\d$/){
            $id = $1;
        }
        $$headers{$id} = '1';
    }

    return $headers;
}

sub FastqReader {
    my ($file) = @_;

    my ( $open_flg, $finished ) = ( 0, 0 );
    my ( $fh, $head, $seq, $qual ) = (undef) x 4;

    if ( $file =~ /^STDIN$/i ) {    # from stdin
        $fh = *STDIN;
    }
    elsif ( ref $file eq '' or ref $file eq 'SCALAR' ) {    # from file
        open $fh, '<', $file or die "fail to open file: $file!\n";
        $open_flg = 1;
    }
    else {    # glob, i.e. given file handler
        $fh = $file;
    }

    return sub {
        return if $finished;

        while (<$fh>) {
            if ( substr( $_, 0, 1 ) ne '@' ) {
                die "bad fq file\n";
            }
            
            $head = $_;
            $head =~ s/\r?\n$//;
            substr( $head, 0, 1, '' );

            $seq = <$fh>;
            $seq =~ s/\r?\n$//;

            <$fh>;
            
            $qual = <$fh;
            $qual =~ s/\r?\n$//;

            return [ $head, $seq, $qual ];
        }

        close $fh if $open_flg;
        $finished = 1;
        return;
    };
}
