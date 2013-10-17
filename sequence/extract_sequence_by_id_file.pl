#!/usr/bin/env perl

# Function: Given a id file, extracting records in another file.
#           It works well for super big file.
# Author  : Wei Shen <shenwei356#gmail.com> http://shenwei.me
# Date    : 2013-08-01
# Docment : http://blog.shenwei.me/extract_records_by_id_file/

use strict;

die "Usage: $0 <id_file> <seq_file> <out_file>\n"
  unless @ARGV == 3;

my $id_file  = shift;
my $seq_file = shift;
my $out_file = shift;

#-------------[ read ids ]-------------

my %ids_hash;                  # 用字典来存储id，后面查询的时候效率更高

open ID, "<", $id_file
  or die "Failed to open file $id_file.\n";
while (<ID>) {
    s/\r?\n//;                 # 记得把回车\r和换行符\n删掉
    
    # 根据具体情况提取id !!!!!!
    # next unless /gi\|(\d+)/; # gi|12313|的情况
    next unless /(.+)/;        # 整个一行作为id的情况
    
    $ids_hash{$1} = 1;         # 加入字典
}
close ID;

# show number of ids
my @ids = keys %ids_hash;
my $n   = @ids;
print "\nRead $n ids.\n\n";

#-------------[ searching ]-------------

# 显示搜索进度的变量，当目标文件非常大的时候很有用
my $count = 0; # 当前处理的序列数
my $hits  = 0; # 匹配到的序列数
local $|  = 1; # 输出通道在每次打印或写之后都强制刷新，提高显示进度速度

open NT, "<", $seq_file
  or die "Failed to open file $seq_file.\n";
open OUT, ">", $out_file
  or die "Failed to open file $out_file.\n";

local $/ = ">";       # 用local是处于安全考虑，怕影响到后续读文件。
                      # $/ 为输入记录分隔符，默认为换行符
<NT>;                 # 略过第一个“>”
my ( $head, $seq );
while (<NT>) {
    $count++;
    s/\r?\n>//;
    ( $head, $seq ) = split "\n", $_, 2;
    $seq =~ s/\s+//g;
    
    # 根据具体情况提取id !!!!!!!!!!!!!!!!!!!!!
    # 取出记录中的id
    # next unless $head =~ /gi\|(\d+)\|/;  # gi|12313|的情况
    # next unless $head =~ /(.+)/;         # 整个一行作为id的情况
    next unless $head =~ /(.+?)_/;         # 我测试的例子，勿套用
    
    # 在%ids_hash中查询记录
    if ( exists $ids_hash{$1} ) {
        print OUT ">$head\n$seq\n";
        
        # 如果确信目标文件中只有唯一与ID匹配的记录，则从字典中删除，提高查询速度
        # delete $ids_hash{$1}; 
        $hits++;
    }
    print "\rProcessing ${count} th record. hits: $hits";
}
$/ = "\n";
close NT;
close OUT;

# 显示没有匹配到任何记录的id
my @ids = keys %ids_hash;
my $n = @ids;
print "\n\n$n ids did not match any record in $seq_file:\n";
print "@ids\n";