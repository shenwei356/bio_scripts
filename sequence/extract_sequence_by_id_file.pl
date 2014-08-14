#!/usr/bin/env perl

# Function: Given a id file, extracting records in another file.
#           It works well for super big file.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    : 2013-08-01
# Update  : 2014-08-06
# Docment : http://blog.shenwei.me/extract_records_by_id_file/

use strict;

die "Usage: $0 <id_file> <seq_file> <out_file>\n"
    unless @ARGV == 3;

my $id_file  = shift;
my $seq_file = shift;
my $out_file = shift;

#-------------[ read ids ]-------------

my %ids_hash
    ;    # 用字典（查询效率更高）来存储id及每个id的命中数

open ID, "<", $id_file
    or die "Failed to open file $id_file.\n";
while (<ID>) {
    s/\r?\n//;    # 记得把回车\r和换行符\n删掉
    next if /^\s*$/;
    s/^\s+|\s+$//;

    # 根据具体情况提取id !!!!!!
    # next unless /gi\|(\d+)/; # gi|12313|的情况
    next unless /(.+)/;    # 整个一行作为id的情况

    $ids_hash{$1} = 0;     # 加入字典
}
close ID;

# show number of ids
my @ids = keys %ids_hash;
my $n   = @ids;
print "\nRead $n ids.\n\n";

#-------------[ searching ]-------------

# 显示搜索进度的变量，当目标文件非常大的时候很有用
my $count = 0;    # 当前处理的序列数
my $hits  = 0;    # 匹配到的序列数
local $| = 1
    ; # 输出通道在每次打印或写之后都强制刷新，提高显示进度速度

open OUT, ">", $out_file
    or die "Failed to open file $out_file.\n";

my $next_seq = FastaReader($seq_file);
while ( my $fa = &$next_seq() ) {
    my ( $head, $seq ) = @$fa;

    $count++;

    $seq =~ s/\s+//g;

    # 根据具体情况提取id !!!!!!!!!!!!!!!!!!!!!
    # 取出记录中的id
    # next unless $head =~ /gi\|(\d+)\|/;  # gi|12313|的情况
    # next unless $head =~ /(.+?)_/;       # 我测试的例子，勿套用
    next unless $head =~ /(.+)/;    # 整个一行作为id的情况

    # 在%ids_hash中查询记录
    if ( exists $ids_hash{$1} ) {
        print OUT ">$head\n$seq\n";

# 如果确信目标文件中只有唯一与ID匹配的记录，则从字典中删除，提高查询速度
# delete $ids_hash{$1};

        # record hit number of a id
        $ids_hash{$1}++;
        $hits++;
    }
    print "\rProcessing ${count} th record. hits: $hits";
}
close OUT;

# 显示没有匹配到任何记录的id
my @ids = grep { $ids_hash{$_} == 0 } keys %ids_hash;
my $n = @ids;
print "\n\n$n ids did not match any record in $seq_file:\n";
print "@ids\n";

# FastaReader is a fasta file parser using closure.
# FastaReader returns an anonymous subroutine, when called, it
# will return a fasta record which is reference of an array
# containing fasta header and sequence.
#
# A boolean argument is optional. If set as "true", "return" ("\r") and
# "new line" ("\n") symbols in sequence will not be trimed.
#
# Example:
#
#    # my $next_seq = FastaReader("test.fa", 1);
#    my $next_seq = FastaReader("test.fa");
#
#    while ( my $fa = &$next_seq() ) {
#        my ( $header, $seq ) = @$fa;
#
#        print ">$header\n$seq\n";
#    }
#
sub FastaReader {
    my ( $file, $not_trim ) = @_;

    my ( $last_header, $seq_buffer ) = ( '', '' ); # buffer for header and seq
    my ( $header,      $seq )        = ( '', '' ); # current header and seq
    my $finished = 0;

    open FH, "<", $file
        or die "fail to open file: $file!\n";

    return sub {

        if ($finished) {                           # end of file
            return undef;
        }

        while (<FH>) {
            s/^\s+//;    # remove the space at the front of line

            if (/^>(.*)/) {    # header line
                ( $header, $last_header ) = ( $last_header, $1 );
                ( $seq,    $seq_buffer )  = ( $seq_buffer,  '' );

                # only output fasta records with non-blank header
                if ( $header ne '' ) {
                    $seq =~ s/\s+//g unless $not_trim;
                    return [ $header, $seq ];
                }
            }
            else {
                $seq_buffer .= $_;    # append seq
            }
        }
        close FH;
        $finished = 1;

        # last record
        # only output fasta records with non-blank header
        if ( $last_header ne '' ) {
            $seq_buffer =~ s/\s+//g unless $not_trim;
            return [ $last_header, $seq_buffer ];
        }
    };
}
