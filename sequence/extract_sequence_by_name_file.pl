#!/usr/bin/perl

use strict;

# try to use BioUtil::Seq
if ( eval { require BioUtil::Seq; require BioUtil::Util; 1; } ne 1 ) {
    die "\nPlease install BioUtil::Seq by CPAN:\n"
        . "  cpan install BioUtil\n\n";
}
else {
    BioUtil::Seq->import();
    BioUtil::Util->import();
}

usage() unless @ARGV == 2;
my $name_file = shift;
my $seq_file   = shift;

my $names = get_list_from_file($name_file);
my %data = map {$_ => 1} @$names;

my $next_seq = FastaReader($seq_file, 1);

while ( my $fa = &$next_seq() ) {
    my ( $head, $seq ) = @$fa;
    if ( defined $data{$head} ) {
        print ">$head\n$seq";
    }

}

sub usage {
    die qq(
Usage: $0 <name_file> <sequence_file> 
    
);
}
