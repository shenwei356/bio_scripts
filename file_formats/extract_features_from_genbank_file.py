#!/usr/bin/env python
# https://github.com/shenwei356/bio_scripts

import sys
import argparse
import gzip

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract features from Genbank file",
        epilog="https://github.com/shenwei356/bio_scripts")

    parser.add_argument('gbkfile', type=str, help='Genbank file')
    parser.add_argument(
        '-t',
        '--type',
        type=str,
        default='CDS',
        help='Feature type (CDS tRNA). Multiple values should be separated by comma. "." for any types.')
    outfmt_choices = ['fasta', 'gtf', 'gff']
    parser.add_argument('-f',
                        '--outfmt',
                        type=str,
                        default='fasta',
                        help='Out format, fasta or gtf')

    parser.add_argument('-p',
                        '--peptide',
                        action="store_true",
                        help='Translate the nucleotides to peptides')
    parser.add_argument(
        '--table',
        type=int,
        help='Genetic code table (detail: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi ) [1]')

    args = parser.parse_args()

    if args.outfmt not in outfmt_choices:
        sys.stderr.write('[ERROR] -f | --outfmt should be in {}\n'.format(outfmt_choices))
        sys.exit(1)

    if args.table:
        args.peptide = True

    return args


if __name__ == '__main__':
    args = parse_args()

    types = set(args.type.lower().split(','))
    with gzip.open(args.gbkfile) if args.gbkfile.endswith('.gz') else open(args.gbkfile) as fh:
        records = SeqIO.parse(fh, "genbank")
        for record in records:
            for f in record.features:
                if '.' not in types and f.type.lower() not in types:
                    continue

                start, end = f.location.start, f.location.end
                strand = '+' if f.strand > 0 else '-'

                qualifiers = f.qualifiers
                if 'product' in qualifiers:
                    product = qualifiers['product'][0]
                else:
                    product = ''

                if 'note' in qualifiers:
                    note = qualifiers['note'][0]
                else:
                    note = ''

                if 'gene' in qualifiers:
                    gene_id = qualifiers['gene'][0]
                elif 'locus_tag' in qualifiers:
                    gene_id = qualifiers['locus_tag'][0]
                else:
                    gene_id = ''

                if args.outfmt == 'fasta':
                    seq = None
                    if args.peptide:
                        if args.table:
                            transl_table = args.table
                        elif 'transl_table' in qualifiers:
                            transl_table = qualifiers['transl_table']
                        else:
                            sys.stderr.write('[WARNING] neither translate table given or found in features. set 1\n')
                            transl_table = 1

                        if 'translation' in qualifiers:
                            seq = Seq(qualifiers['translation'][0])
                        else:
                            seq = record.seq[start:end].translate(table=transl_table)
                    else:
                        seq = record.seq[start:end]

                    SeqIO.write(
                        [SeqRecord(seq,
                                   id='{}_{}..{}..{}'.format(record.id, start + 1, end, strand), description=product)],
                        sys.stdout,
                        "fasta")

                elif args.outfmt == 'gtf':
                    frame = int(qualifiers['codon_start'][0]) - 1 if 'codon_start' in qualifiers else 0

                    transcript_id = gene_id

                    attribute = 'gene_id "{}"; transcript_id "{}"'.format(gene_id, transcript_id)

                    if 'protein_id' in f.qualifiers:
                        attribute += '; protein_id "{}"'.format(qualifiers['protein_id'][0])

                    if 'db_xref' in qualifiers:
                        for ext in qualifiers['db_xref']:
                            attribute += '; db_xref "{}"'.format(ext)

                    if 'note' in f.qualifiers:
                        attribute += '; note "{}"'.format(qualifiers['note'][0])

                    attribute += '; product "{}"; '.format(product)

                    sys.stdout.write('\t'.join(
                        [record.id, 'genbank', f.type, str(start + 1), str(end), '.', strand, str(frame), attribute]) + "\n")

                elif args.outfmt == 'gff':
                    if 'codon_start' in qualifiers:
                        frame = int(qualifiers['codon_start'][0]) - 1
                    else:
                        frame = 0
                    sys.stdout.write('\t'.join(
                        [record.id, 'genbank', f.type, str(start + 1),
                         str(end), '.', strand, str(frame),
                         "{},{}".format(gene_id, product)]) + "\n")
