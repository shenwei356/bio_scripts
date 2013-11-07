#!/usr/bin/env perl

# Function: Delete illegal charactors of head line in fasta file before blast.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    :

die "Usage: $0 fasta_file\n" unless @ARGV > 0;

while (@ARGV) {
    $file = shift @ARGV;
    $file =~ /(.+)\..+?$/;
    $out  = $1;
    open IN, $file or die "in_file $file failed to open\n";
    open OUT, ">$out.rename.fa" or die "out_file $out failed to open\n";
    while (<IN>) {
        if (/^>/) {
            s/[^\>\r\n\w\.\ \-]+/_/g; 
        }
        $buffer1 = $_;
        
        $buffer .= $buffer1;
        $buffer_size += length $buffer1;
        if ( $buffer_size > 30_000_000 ) {      #buffer
            print OUT $buffer;
            $buffer = '';
            $buffer_size = 0;
        }
    }
    close IN;
    print OUT $buffer;
    close OUT;
}