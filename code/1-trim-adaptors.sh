#!bin/bash
#trim adaptors
#Cutadapt cutadapt 2.10 with Python 3.6.11 was used to trim F341 (CCTACGGGNGGCWGCAG)
#and R785 (GACTACHVGGGTATCTAATCC) primers (Klindworth et al. 2013) from R1 and R2 reads respectively.
#Run loop to (1) trim forward primer anchored to 5' end in R1 reads and
# reverse primer anchored to 5' end in R2 reads
# write to ouput with only paired sequences for which both primers were found:

find data/raw/*R1.fastq.gz | xargs basename | sed 's/_R1.fastq.gz//' | while read reads
do
  cutadapt -g ^CCTACGGGNGGCWGCAG -G ^GACTACHVGGGTATCTAATCC --trimmed-only \
  -o data/intermediate/"$reads"_R1.cutadapt.fastq.gz \
  -p data/intermediate/"$reads"_R2.cutadapt.fastq.gz \
  data/raw/"$reads"_R1.fastq.gz data/raw/"$reads"_R2.fastq.gz
done > output/cutadatp_filtering.txt
