#!/usr/bin/env bash
# getReads.sh

# Retrieve the NGS reads from the NA12878 reference sample

fastq-dump --split-files SRR6808334 1>getReads.log 2>getReads.err

# Do not download on Defiance!!
# If you want these files, do
# symbolic links here

# ln -s /data/METHODS/Spring/Module05/SRR6808334_1.fastq SRR6808334_1.fastq
# ln -s /data/METHODS/Spring/Module05/SRR6808334_2.fastq SRR6808334_2.fastq