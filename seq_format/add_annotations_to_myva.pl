#!/usr/bin/perl
# Fuction     : To add annotations to COG aa file
# Author      : Wei Shen
# Email       : shenwei356@gmail.com
# Date        : 2011-04-08, cost 2 hour.
# Last Update : 2011-04-08

# Annotations are in following files downloaded from
# ftp://ftp.ncbi.nlm.nih.gov/pub/COG/COG
# 
# FILE       DATA(* means important)
# *fun.txt    function_id(one letter) -> function
# *myva      protein_id*, aa sequence*
# *myva=gb   protein_id  -> GI numbers*
# *whog      Organism(three letters)*, protein_id, 
#            function_id, cog_id, protein*
# *org.txt   Organism(three letters) -> detail
#
# Output format
# > <protein_id>_<protein>_<function>_<Organism>_<GI numbers>
# aa sequence
#
# Therefore, five files with sign * will be used
#
# Attention:
# 1. NOT all protein_ids from file myva=gb could be found in file whog
# 2. NOT all protein_ids from file myva    could be found in file whog

use strict;

# parse file fun
my $fun = &parse_file_fun('fun.txt');


# parse file myva=gb
my $pro_gi = &parse_file_myva_gb('myva=gb');


# parse file org
my $org = &parse_file_org('org.txt');


# parse file whog
my $whog = &parse_file_whog('whog');

# my @keys = keys %$whog;
# print "$_\n" unless $_ ~~ @keys for keys %$pro_gi;
# the result showed that not all protein_ids from file myva=gb
# are in file whog

# parse file myva and add annotation
my $file = 'myva';
my $out_file = "$file.full_annotation.txt";
my ($head, $seq, $pro_id_trim);

open IN, $file or die "File $file failed to open sequence.\n";
$/ = '>';<IN>;
open OUT, ">", $out_file or die "File $out_file failed to open sequence.\n";

while ( <IN> ) {
    s/\r?\n>//;
    ( $head, $seq ) = split /\r?\n/, $_, 2;
    ## > <protein_id>_<protein>_<function>_<Organism>_<GI numbers>
    $pro_id_trim = $head;
    $pro_id_trim = $1 if $head =~ /(.+)\_\d+/;    # for gi
    $head = $head
            . " __pro__". $$whog{$head}{protein}
            . "__fun_". $$whog{$head}{fun_id}. "_". $$fun{$$whog{$head}{fun_id}}
            . "__org__". $$org{$$whog{$head}{org_id}}{organism}
            . "__gi__".  $$pro_gi{$pro_id_trim};
    print OUT ">$head\n$seq\n";
}
close IN;
close OUT;

#====================================================================
# out put data structure:
# $hash_ref = {fun_id => function}
sub parse_file_fun($){
    my ($file) = @_;
    my $fun = {};
    
    open IN, $file or die "File $file failed to open\n";
    while (<IN>) {
        next unless /\[(\w)\] (.+) $/;
        $$fun{$1} = $2;
    }
    close IN;
    # print scalar keys %$fun;
    return $fun;
}

# out put data structure:
# $hash_ref = {protein_id => gi}
sub parse_file_myva_gb($){
    my ($file) = @_;
    my $pro_gi = {};
    
    open IN, $file or die "File $file failed to open\n";
    while (<IN>) {
        next unless /^(.+)\s+(.+)$/;
        $$pro_gi{$1} = $2;
    }
    close IN;
    # print scalar keys %$pro_gi;
    return $pro_gi;
}

# out put data structure:
# $hash_ref = {org_id => {kindom => kindom, organism => organism} }
sub parse_file_org($){
    my ($file) = @_;
    my $org = {};
    
    open IN, $file or die "File $file failed to open\n";
    while (<IN>) {
        next unless /^(\w{3})\s+\d+\s+(.+?)\s+(.+)$/;
        $$org{$1}{kindom}   = $2;
        $$org{$1}{organism} = $3;
    }
    close IN;
    # print scalar keys %$org;
    return $org;
}

# out put data structure:
# $hash_ref = {protein_id => {org_id => org_id, cog_id => cog_id,
#                             fun_id => fun_id, protein => protein} }
sub parse_file_whog($){
    my ($file) = @_;
    my $whog = {};
    
    my ($fun_id, $cog_id, $protein, $org_id, $protein_id, @protein_ids);
    open IN, $file or die "File $file failed to open\n";
    while (<IN>) {
        if (/^\[(\w)\] (\w+) (.+)$/) {
            $fun_id  = $1;
            $cog_id  = $2;
            $protein = $3;
        }
        elsif (/^\s+(\w{3})\:\s+(.+)$/) {
            $org_id     = $1;
            @protein_ids = split /\s+/, $2;
            for $protein_id (@protein_ids) {
                $$whog{$protein_id}{fun_id}  = $fun_id;
                $$whog{$protein_id}{cog_id}  = $cog_id;
                $$whog{$protein_id}{protein} = $protein;
                $$whog{$protein_id}{org_id}  = $org_id;
            }
        }
        elsif (/        (.+)/) {
            @protein_ids = split /\s+/, $1;
            for $protein_id (@protein_ids) {
                $$whog{$protein_id}{fun_id}  = $fun_id;
                $$whog{$protein_id}{cog_id}  = $cog_id;
                $$whog{$protein_id}{protein} = $protein;
                $$whog{$protein_id}{org_id}  = $org_id;
            }
        }
        elsif (/_______/) {
        }
        else {
        }
    }
    close IN;
    # print scalar keys %$whog;
    # print "$_\t". ($$whog{'PH0109_1'}{$_}). "\n" for keys %{$$whog{'PH0109_1'}};
    return $whog;
}
