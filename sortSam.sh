#!/usr/bin/env bash
# sortSam.sh

samtools sort -@ 8 -m 4G SRR6808334.sam -o SRR6808334.sorted.bam \
1>sort.log 2>sort.err

# now that you have a SAM file from bwa, so you you need to sort it, and create
# a binary (sam) = BAM
# -@ 8 number of threads
# -m 4G : Set maximum memory per thread

# If you want these files, do
# symbolic links here

# ln -s /data/METHODS/Spring/Module05/SRR6808334.sam SRR6808334.sam