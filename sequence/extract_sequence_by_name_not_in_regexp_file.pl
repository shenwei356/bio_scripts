#!/usr/bin/env perl
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
my $regex_file = shift;
my $seq_file   = shift;

my $names = get_list_from_file($regex_file);

my $next_seq = FastaReader($seq_file, 1);
while ( my $fa = &$next_seq() ) {
    my ( $head, $seq ) = @$fa;
    my @in = grep { $_ if $head =~ /$_/i } @$names;
    # print "@in\n" if @in > 0;
    if ( @in == 0){
        print ">$head\n$seq";
    }
}

sub usage{
    die qq(
Usage: $0 <regular expression file> <sequence_file> 
    
);
}
