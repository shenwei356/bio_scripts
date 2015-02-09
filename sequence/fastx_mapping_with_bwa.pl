#!/usr/bin/env perl

use strict;
use File::Basename;
use BioUtil::Util;

$0 = basename ($0);
die qq(
    usage: $0 <threads> <refseq> <outprefix> <fasta/q> [<fasta/q>]
    when two fastx file given, they are treated as paired end reads

)
   unless @ARGV == 5 or @ARGV == 4;

my $threads = shift @ARGV;
my $refseq  = shift @ARGV;
my $prefix  = shift @ARGV;
my $read    = shift @ARGV;
my $read2   = shift @ARGV;

check_positive_integer($threads);

# build index
my @suffix      = qw/.amb .ann .bwt .pac .sa/;
my $index_built = 1;
for (@suffix) {
    $index_built = 0 unless -e "$refseq$_";
}
run("bwa index $refseq") unless $index_built;
run("samtools faidx $refseq") unless -e "$refseq.fai";

# =================[ mapping ]===================

print "mapping\n";
if ($read2){
    run("bwa mem -t $threads -M -a $refseq $read $read2 > $prefix.sam");
}else{
run("bwa mem -t $threads -M -a $refseq $read > $prefix.sam");
    }

# =================[ mapping ]===================

print "sam -> bam\n";
run("samtools view  -bS $prefix.sam > $prefix.bam");

print "sort bam\n";
run("samtools sort $prefix.bam $prefix.sorted");

print "index bam\n";
run("samtools index $prefix.sorted.bam");

print "flagstat\n";
run("samtools flagstat $prefix.sorted.bam > $prefix.sorted.bam.flagstat");

run("rm $prefix.bam $prefix.sam");