## data

1M reads

    zcat reads_1.fq.gz | fastq2tab | head -n 1000000 | tab2fastq | gzip -c > reads.fq.gz

100K id

    zcat reads.fq.gz | fastq2tab | cut -f 1 > id
    tac id | head -n 100000  > id2

## single pattern

    time zcat reads.fq.gz | fastq2tab | csv_grep.py   -t -p . > /dev/null

    real    0m15.906s                                                                                      
    user    0m21.220s                                                                                      
    sys     0m1.919s  

    time zcat reads.fq.gz | fastq2tab | csv_grep  -t -p .  > /dev/null
    patterns: 1

    real    0m11.473s
    user    0m18.957s
    sys     0m2.068s


## multiple patterns

    time zcat reads.fq.gz | fastq2tab | csv_grep.py -t -pf id2 > /dev/null

    real    0m16.989s
    user    0m22.192s
    sys     0m1.952s

    time zcat reads.fq.gz | fastq2tab | csv_grep  -t  -pf id2 -j 4  > /dev/null
    patterns: 100000

    real    0m12.904s
    user    0m22.811s
    sys     0m2.243s

## multiple patterns, single thread

    time zcat reads.fq.gz | fastq2tab | csv_grep  -t  -pf id2  -j 1 > /dev/null

    real    0m15.823s
    user    0m21.067s
    sys     0m1.911s
