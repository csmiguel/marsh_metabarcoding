#!/bin/bash
# write sequencing depth of samples
for file in data/raw/S*R1.fastq.gz
do seqkit stats $file
done > output/sequencing_report.txt
