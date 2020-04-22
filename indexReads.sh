#!/usr/bin/env bash
# index.sh

samtools index -@ 8 -m 4G -b SRR6808334.sorted.bam \
1>index.log 2>index.err

# Do not need to run on Defiance
# now that you have a sorted BAM file, it needs to be indexed
# so other programs have fast access
# -@ 8 number of threads
# -m 4G : Set maximum memory per thread


# If you want these files, do
# symbolic links here

# ln -s /data/METHODS/Spring/Module05/SRR6808334.sorted.bam.bai SRR6808334.sorted.bam.bai