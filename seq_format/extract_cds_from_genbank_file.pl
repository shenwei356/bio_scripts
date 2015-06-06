#!/usr/bin/env perl
# by Wei Shen <shenwei356@gmail.com>
#
use strict;
use Bio::SeqIO;
use File::Basename;
# use Data::Dumper;

$0 = basename($0);
die "\nUsage: $0 <genbank_file> <seqtype aa/nt>\n\n"
    unless @ARGV == 2;
my $file    = shift @ARGV;
my $seqtype = shift @ARGV;

die "seqtype should be \"aa\" or \"nt\".\n"
    unless $seqtype eq "aa" or $seqtype eq "nt";

my $cdses = get_cds_from_genbank_file($file);
for my $cds (@$cdses) {
    print ">$$cds{seqid}_$$cds{start}..$$cds{end}..$$cds{strand} $$cds{product}\n"
        . $$cds{$seqtype} . "\n";
}

sub get_cds_from_genbank_file {
    my ( $file, $codontable_id ) = @_;
    $codontable_id = 4 unless defined $codontable_id;

    my $gb = Bio::SeqIO->new( -file => $file, -format => 'GenBank' );
    my $cdses = [];
    my ( $sequence, $features, $feat, $location );
    my ( $seqid, $strand, $start, $end, $cds_seq, $seqid );
    my ( $product, $aa );
    while ( my $seq = $gb->next_seq() ) {    # Bio::Seq
        $sequence = $seq->seq();

        # my @tmp = keys %$sequence;
        # print "@tmp\n";

        $features = $$seq{_as_feat};
        for (@$features) {
            $feat = $_;    # Bio::SeqFeature::Generic
            next unless $feat->primary_tag() eq 'CDS';
            $location = $feat->location();     # Bio::Location::Simple
            $strand   = $location->strand();
            $start    = $location->start();
            $end      = $location->end();
            $seqid    = $location->{_seqid};
            $cds_seq = substr( $sequence, $start - 1, $end - $start + 1 );

            $product = @{ $feat->{'_gsf_tag_hash'}->{product} }[0];

            if ( $strand ne '1' ) {    # reversed complement of the sequence
                my $tmp_seq = Bio::Seq->new( -seq => $cds_seq );
                $cds_seq = $tmp_seq->revcom()->seq();
            }

            $aa = Bio::Seq->new( -seq => $cds_seq )
                ->translate( -codontable_id => $codontable_id )->seq();
            $aa =~ s/\*$//g;
            push @$cdses,
                {
                strand  => $strand,
                start   => $start,
                end     => $end,
                nt      => $cds_seq,
                seqid   => $seqid,
                aa      => $aa,
                product => $product
                };
        }
    }
    return $cdses;
}
