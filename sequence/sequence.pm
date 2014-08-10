# FastaReader is a fasta file parser, which returns a function that
# returns a pair of head and sequence when it was called
#
# Example:
#
#    my $next_seq = FastaReader("test.fa");
#
#    my ( $head, $seq );
#    while (1) {
#        ( $head, $seq ) = &$next_seq();
#        last
#          if $head eq "" and $seq eq "";
#        print ">$head\n$seq\n";
#    }
sub FastaReader {
    my ($file) = @_;
    open IN, "<", $file
        or die "fail to open file: $file!\n";
    local $/ = '>';
    <IN>;
    $/ = '\n';

    my ( $line, $head, $seq );
    return sub {
        local $/ = '>';
        while ( $line = <IN> ) {
            $line =~ s/\r?\n>?$//;
            ( $head, $seq ) = split /\r?\n/, $line, 2;
            $seq =~ s/\s+//g;
            return ( $head, $seq );
        }
        close IN;
        $/ = "\n";
        return ( "", "" );
    };
}

# Function: reading sequence from fasta file
# Input   : sequence file
# Output  : A hash reference
sub read_sequence_from_fasta_file {
    my ( $seq_file, $not_trim ) = @_;
    my $sequences = {};
    my ( $head, $seq );
    open IN, $seq_file
        or die "fail to open sequence file $seq_file!\n";
    my $/ = '>';
    <IN>;
    while (<IN>) {
        s/\r?\n>?$//;
        ( $head, $seq ) = split "\r?\n", $_, 2;
        $seq =~ s/\s+//g unless $not_trim;
        warn "reduplicate sequence: $head\n" if exists $$sequences{$head};
        $$sequences{$head} = $seq;
    }
    close IN;
    my $/ = "\n";
    die "no sequences in $seq_file\n" unless keys %$sequences > 0;
    return $sequences;
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
