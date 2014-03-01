# https://github.com/shenwei356
# 
# Ths script illustrates how to parse grouped data in multi-line, as below. 
# String of first column is the group ID, and a group may have
# more than one records in multi-line.
# 
#     g1 2 3
#     g1 2 5
#     g2 2 3
#     g2 2 5
#     g3 2 3
#     g3 2 5
#
# Outline
# 
# A flag “last_id” is used to judge first / same / new group (See code bellow). 
# 
# For different situation,
# 
#     1. First record. Initializing container for current group ( id), 
#        and add in this record.   last_id = id
#     2. Same group. Add this record into the container for current group ( id ).
#     2. New group. Do something with previous group ( last_id ). Initializing 
#        container for current group ( id ), and add in this record. last_id = id .
#     2. Last group. Adding last group ( last_id) at the end of file (EOF).
#
# Extension
#
# In previous case, the marker for a new record is a different id. In other cases,
# parsing fasta file for example, the marker is the character “>”.
#
use strict;

my $data = {};  # container for all data
my ( $id, $last_id ) = ( "", "" );
my $record = "";

while (<DATA>) {

    # parse id
    next unless /^(.+?)\s+/;
    $id = $1;

    # parse record. Here is the whole line
    $record = $_;

    if ( $last_id eq "" ) {    # first record
        $$data{$id} = [];      # initialize container for this group
        push @{ $$data{$id} }, $record;     # add this record
        $last_id = $id;                     # restore this id for further use
    }
    else {
        if ( $id eq $last_id ) {            # same group
            push @{ $$data{$id} }, $record; # add this record
        }
        else {                              # new group
            # do something with previous group
            &dosomthing( $$data{$last_id} );

            $$data{$id} = [];
            push @{ $$data{$id} }, $record;
            $last_id = $id;
        }
    }
}

# do something with the last group
&dosomthing( $$data{$id} );

sub dosomthing {
    my ($records) = @_;
    for (@$records) {
        print " $_";
    }
    print "\n";
}

# example data
__DATA__
g1 2 3
g1 2 5
g2 2 3
g2 2 5
g3 2 3
g3 2 5