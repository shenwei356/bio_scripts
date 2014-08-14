

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

# Function: reading sequence from fasta file
# Input   : sequence file
# Output  : A hash reference
sub read_sequence_from_fasta_file {
    my ( $file, $not_trim ) = @_;
    my $seqs = {};

    my $next_seq = FastaReader( $file, $not_trim );
    while ( my $fa = &$next_seq() ) {
        my ( $header, $seq ) = @$fa;

        $$seqs{$header} = $seq;
    }

    return $seqs;
}

sub shuffle_sequences {
    use List::Util qw(shuffle);
    my ( $file, $not_trim ) = @_;
    my $seqs = read_sequence_from_fasta_file( $file, $not_trim );
    my @keys = shuffle( keys %$seqs );

    my $file_out = "$file.shuffled";
    open OUT, ">$file_out" or die "fail to write file $file_out\n";
    print OUT ">$_\n$$seqs{$_}\n" for @keys;
    close OUT;

    return $file_out;
}

sub write_sequence_to_fasta_file($$) {
    my ( $seqs, $file ) = @_;
    open OUT, ">$file" or die "failed to write to $file\n";
    for ( keys %$seqs ) {
        print OUT ">$_\n", format_seq( $$seqs{$_}, 60 ), "\n";
    }
    close OUT;
}

#here, the seq1 and seq2 have same length
sub is_sequence_similar($$$) {
    my ( $seq1, $seq2, $threshold ) = @_;
    $seq1 = uc $seq1;
    $seq2 = uc $seq2;
    my $n     = length($seq1);
    my $score = 0;
    for ( 0 .. ( $n - 1 ) ) {
        $score++ if substr( $seq1, $_, 1 ) eq substr( $seq2, $_, 1 );
    }

    if ( $score / $n > $threshold ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub revcom {
    my ($s) = @_;
    $s =~ tr/ACGTRYMKSWBDHVNacgtrymkswbdhvn/TGCAYRKMWSVHDBNtgcayrkmwsvhdbn/;
    return reverse $s;
}

sub generate_random_seqence($$) {
    my ( $alphabet, $length ) = @_;
    my $n = @$alphabet;
    my $seq;
    $seq .= $$alphabet[ int rand($n) ] for ( 1 .. $length );
    return $seq;
}

sub format_seq($$) {
    my ( $s, $n ) = @_;
    $n = 60 unless defined $n;
    my $s2 = '';
    my ( $j, $int );
    $int = int( ( length $s ) / $n );
    for ( $j = 0; $j <= $int; $j++ ) {
        $s2 .= substr( $s, $j * $n, $n ) . "\n";
    }
    return $s2;
}

sub edit_str_for_file_name($) {
    my $s = $_[0];
    $s =~ s/[\\\/\:\*\?\"\<\>\|]/\_/g;
    return $s;
}

sub base_content($$) {
    my ( $seq, $bases ) = @_;
    if ( $seq eq '' ) {
        return;
    }
    if ( $bases =~ /[^acgturymkswbdhvn]/i ) {
        warn "wrong bases given!\n";
        return;
    }
    my $sum = 0;
    $sum += $seq =~ s/$_/$_/ig for split "", $bases;
    return sprintf "%.4f", $sum / length $seq;
}

sub validate_sequence {
    my ($seq) = @_;
    return 0 if $seq =~ /[^ATCGRYMKSWHBVDN]/i;
    return 1;
}

# A subroutine to translate DNA sequence into a peptide
sub dna2peptide {
    my($dna) = @_;
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i = 0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }
    return $protein;
}

# A subroutine to translate a DNA 3-character codon to an amino acid
sub codon2aa {
    my($codon) = @_;
    $codon = uc $codon;
    my %genetic_code = (
        'TCA' => 'S',    # Serine
        'TCC' => 'S',    # Serine
        'TCG' => 'S',    # Serine
        'TCT' => 'S',    # Serine
        'TTC' => 'F',    # Phenylalanine
        'TTT' => 'F',    # Phenylalanine
        'TTA' => 'L',    # Leucine
        'TTG' => 'L',    # Leucine
        'TAC' => 'Y',    # Tyrosine
        'TAT' => 'Y',    # Tyrosine
        'TAA' => '_',    # Stop
        'TAG' => '_',    # Stop
        'TGC' => 'C',    # Cysteine
        'TGT' => 'C',    # Cysteine
        'TGA' => '_',    # Stop
        'TGG' => 'W',    # Tryptophan
        'CTA' => 'L',    # Leucine
        'CTC' => 'L',    # Leucine
        'CTG' => 'L',    # Leucine
        'CTT' => 'L',    # Leucine
        'CCA' => 'P',    # Proline
        'CCC' => 'P',    # Proline
        'CCG' => 'P',    # Proline
        'CCT' => 'P',    # Proline
        'CAC' => 'H',    # Histidine
        'CAT' => 'H',    # Histidine
        'CAA' => 'Q',    # Glutamine
        'CAG' => 'Q',    # Glutamine
        'CGA' => 'R',    # Arginine
        'CGC' => 'R',    # Arginine
        'CGG' => 'R',    # Arginine
        'CGT' => 'R',    # Arginine
        'ATA' => 'I',    # Isoleucine
        'ATC' => 'I',    # Isoleucine
        'ATT' => 'I',    # Isoleucine
        'ATG' => 'M',    # Methionine
        'ACA' => 'T',    # Threonine
        'ACC' => 'T',    # Threonine
        'ACG' => 'T',    # Threonine
        'ACT' => 'T',    # Threonine
        'AAC' => 'N',    # Asparagine
        'AAT' => 'N',    # Asparagine
        'AAA' => 'K',    # Lysine
        'AAG' => 'K',    # Lysine
        'AGC' => 'S',    # Serine
        'AGT' => 'S',    # Serine
        'AGA' => 'R',    # Arginine
        'AGG' => 'R',    # Arginine
        'GTA' => 'V',    # Valine
        'GTC' => 'V',    # Valine
        'GTG' => 'V',    # Valine
        'GTT' => 'V',    # Valine
        'GCA' => 'A',    # Alanine
        'GCC' => 'A',    # Alanine
        'GCG' => 'A',    # Alanine
        'GCT' => 'A',    # Alanine
        'GAC' => 'D',    # Aspartic Acid
        'GAT' => 'D',    # Aspartic Acid
        'GAA' => 'E',    # Glutamic Acid
        'GAG' => 'E',    # Glutamic Acid
        'GGA' => 'G',    # Glycine
        'GGC' => 'G',    # Glycine
        'GGG' => 'G',    # Glycine
        'GGT' => 'G',    # Glycine
        );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }
    else{
        print STDERR "Bad codon \"$codon\"!!\n";
        exit;
    }
}

1;