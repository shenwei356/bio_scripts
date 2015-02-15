# Utilities 

## csv_grep

grep CSV file, including tab-delimited file, by excat matching or regluar
expression. The query pattern could be given from command line or file.
Column number of key in target file and pattern file is settable.

### Usage:

	usage: csv_grep [-h] [-v] [-o [OUTFILE]] [-k KEY] [-H] [-F FS] [-Q QC]
					[-p [PATTERN]] [-pf [PATTERNFILE]] [-pk [PK]] [-r] [-s] [-n]
					[csvfile [csvfile ...]]

	Grep CSV file

	positional arguments:
	csvfile               Input file(s)

	optional arguments:
	-h, --help            show this help message and exit
	-v, --verbose         Verbosely print information
	-o [OUTFILE], --outfile [OUTFILE]
							Output file [STDOUT]
	-k KEY, --key KEY     Column number of key in csvfile
	-H, --ignoretitle     Ignore title
	-F FS, --fs FS        Field separator [\t]
	-Q QC, --qc QC        Quote char["]
	-p [PATTERN], --pattern [PATTERN]
							Query pattern
	-pf [PATTERNFILE], --patternfile [PATTERNFILE]
							Pattern file
	-pk [PK]              Column number of key in pattern file
	-r, --regexp          Pattern is regular expression
	-s, --speedup         Delete matched pattern, if you know what it means
	-n, --invert          Invert match (do not match)

### Example

1) Grep table by a reglur expression, the first column is the key.

	csv_grep -p '^\d+$' -r file.txt

2) Find common record of two fasta, with same headers. 
[*fasta2tab*](https://github.com/shenwei356/bio_scripts/blob/master/sequence/fasta2tab) 
 transforms the FASTA fromat to two-column table, fist column is the header and the second is sequence. 
[*tab2fasta*](https://github.com/shenwei356/bio_scripts/blob/master/sequence/tab2fasta) just tranform the
table back to FASTA format.
	
	fasta2tab seq1.fa | csv_grep -pf <(fasta2tab seq.fa) | tab2fasta
	
With same sequence (Second column)

	fasta2tab seq1.fa | csv_grep -pf <(fasta2tab seq.fa) -pk 2  -k 2  | tab2fasta

3) Find common records of two gtf file.
The columns 1,4,5,7 together make up the key of a record,
so before grep, we add a new column as the key.

	awk -F"\t" '{print $0"\t"$1""$4""$5""$7}' c.gtf > c1.gtf
	awk -F"\t" '{print $0"\t"$1""$4""$5""$7}' d.gtf > d1.gtf
	cat c1.gtf | python3 csv_grep -k 10 -pf d1.gtf -pk 10 | awk 'NF-=1' > common.gtf
	
## plot_distribution.py
