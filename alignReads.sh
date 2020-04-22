#!/usr/bin/env bash
# alignReads.sh

bwa mem -M -t 8 -R "@RG\tID:SRR6808334\tSM:bar" \
GRCh38_reference.fa \
SRR6808334_1.paired.fastq SRR6808334_2.paired.fastq \
1>SRR6808334.sam 2>alignReads.err &

# -R just adding read groups (RG) ID = read group identifier (ID)
# read group sample (SM), required if ID is given, note is "bar" here
# -M mark shorter split hits as secondary
# -t number of threads


# Technically, we'd want to do two more alignments too! Using commands abovle

# bwa mem -M -t 8 -R "@RG\tID:SRR6808334\tSM:bar" \
# GRCh38_reference.fa \
# SRR6808334_1.unpaired.fastq \
# 1>SRR6808334_1.unpaired.sam 2>alignReads2.err &

# bwa mem -M -t 8 -R "@RG\tID:SRR6808334\tSM:bar" \
# GRCh38_reference.fa \
# SRR6808334_2.unpaired.fastq \
# 1>SRR6808334_2.unpaired.sam 2>alignReads3.err &

# WHy?  Reads where one mate fails but are otherwise fine are included in a single-end file. 
# It's still good data! Map them the same way you do the PE data, and subsequently 
# combine the three Same files above using Picard's MergeSamFiles or another approach.
# but we'll just do PE reads here