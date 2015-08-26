# Utilities 

## intersection

Intersecion of multiple (>=2) files.

## csv_grep

Grepping CSV file, tab-delimited file by default, by exactly matching or query by
regluar expression, multiple keys (indice) supported. The query patterns could be given from command line or file.

### Usage:

    usage: csv_grep [-h] [-v] [-o OUTFILE] [-k KEY] [-H] [-F FS] [-Fo FS_OUT]
                    [-Q QC] [-t] [-p [PATTERN]] [-pf [PATTERNFILE]] [-pk [PK]]
                    [-r] [-d] [-i]
                    [csvfile [csvfile ...]]
    
    Grep CSV file. Multiple keys supported.
    
    positional arguments:
      csvfile               Input file(s)
    
    optional arguments:
      -h, --help            show this help message and exit
      -v, --verbose         Verbosely print information
      -o OUTFILE, --outfile OUTFILE
                            Output file [STDOUT]
      -k KEY, --key KEY     Column number of key in csvfile. Multiple values shoud
                            be separated by comma
      -H, --ignoretitle     Ignore title
      -F FS, --fs FS        Field separator [,]
      -Fo FS_OUT, --fs-out FS_OUT
                            Field separator of ouput [same as --fs]
      -Q QC, --qc QC        Quote char["]
      -t                    Field separator is "\t". Quote char is "\t"
      -p [PATTERN], --pattern [PATTERN]
                            Query pattern
      -pf [PATTERNFILE], --patternfile [PATTERNFILE]
                            Pattern file
      -pk [PK]              Column number of key in pattern file. Multiple values
                            shoud be separated by comma
      -r, --regexp          Pattern is regular expression
      -d, --speedup         Delete matched pattern when matching one record
      -i, --invert          Invert match (do not match)
    
    https://github.com/shenwei356/bio_scripts


### Examples

1) For a table file. Note that the 3rd column of 4th line contains "\t".

	$ cat testdata/data.tab
	column1 column 2        3rd c
	str     123     abde
	123     134     我
	245     135     "string with    tab"

Find lines of which the 2nd column are digitals, ignoring title

	$ cat testdata/data.tab | csv_grep -H  -t -k 2 -r -p '^\d+$'
	str     123     abde
	123     134     我
	245     135     "string with    tab"

	
Find lines that have ID (first column, by default) in (or NOT in) a given ID files.

	$ cat testdata/data.tab | csv_grep -t -pf testdata/data.pattern.tab
    123     134     我
    245     135     "string with    tab"

	
	$ cat testdata/data.tab | csv_grep -H -t -pf testdata/data.pattern.tab -i
    str     123     abde


2) Find common records with same headers in two fasta files. 
[*fasta2tab*](https://github.com/shenwei356/bio_scripts/blob/master/sequence/fasta2tab) 
 transforms the FASTA fromat to two-column table, fist column is the header and the second is sequence. 
[*tab2fasta*](https://github.com/shenwei356/bio_scripts/blob/master/sequence/tab2fasta) just tranform the
table back to FASTA format.
	
	fasta2tab seq1.fa | csv_grep -t -pf <(fasta2tab seq.fa) | tab2fasta
	
Records with same sequence (second column).

	fasta2tab seq1.fa | csv_grep -t -pf <(fasta2tab seq.fa) -pk 2  -k 2  | tab2fasta

3) Find common records of two GTF file.
The columns 1,4,5,7 together make up the key of a record.

    cat a.gff | csv_grep -t -k 1,4,5,7 -pk 1,4,5,7 -pf b.gff > commom.gff
	
## csv_grep2
    
    golang version. Faster with multi-threads.
    
### Usage:
    
    NAME:
       csv_grep - grep for csv format
    
    USAGE:
       csv_grep [global options] command [command options] [arguments...]
       
    VERSION:
       1.0
       
    AUTHOR(S):
       Wei Shen <https://github.com/shenwei356/bio_scripts> 
       
    COMMANDS:
       help, h      Shows a list of commands or help for one command
       
    GLOBAL OPTIONS:
       -k, --key "1"                column number of key in csvfile. Multiple values shoud be separated by comma [1]
       -H, --ignoretitle            ignore title
       --fs ","                     field separator [,]
       --fs-out                     field separator of ouput [same as --fs]
       -t, --tab                    field separator is "\t". Quote char is "\t"
       -p, --pattern                query pattern
       --pf, --patternfile          pattern file
       --pk "1"                     column number of key in pattern file. Multiple values shoud be separated by comma [1]
       --pfs ","                    field separator of pattern file [,]
       -r, --use-regexp             use regular expression
       -d, --speedup                delete matched pattern when matching one record
       -i, --invert                 invert match (do not match)
       -j, --ncpus "4"              CPU number [4]
       -c, --chunksize "1000"       chunk size [1000]
       -o, --outfile                output file [stdout]
       --vv, --verbose              verbosely print information
       --help, -h                   show help
       --version, -v                print the version
    
    
	
## csv_join

### Usage

    usage: csv_join [-h] [-f1 F1] [-q1 Q1] [-f2 F2] [-q2 Q2] [-of OF] [-t1] [-t2]
                    [-to] [-t] [-k]
                    csvfile1 key1 csvfile2 key2
    
    Merge csvfile2 to csvfile1. Multiple keys supported.
    
    positional arguments:
      csvfile1              CSV file 1
      key1                  Column number of key in csvfile1. Multiple values
                            shoud be separated by comma.
      csvfile2              CSV file 2
      key2                  Column number of key in csvfile2. Multiple values
                            shoud be separated by comma.
    
    optional arguments:
      -h, --help            show this help message and exit
      -f1 F1                Field separator in csvfile1 [,]
      -q1 Q1                Quote char in csvfile1 ["]
      -f2 F2                Field separator in csvfile2 [,]
      -q2 Q2                Quote char in csvfile2 ["]
      -of OF                Field separator in csvfile1 [,]
      -t1                   csvfile1 is table file. Quote char is "\t"
      -t2                   csvfile1 is table file. Quote char is "\t"
      -to                   Output quote char is "\t"
      -t                    abbr for "-t1 -t2 -to"
      -k, --keep-unmatched  keep rows in CSV file1 not matching row in file2"
    
    https://github.com/shenwei356/bio_scripts

    