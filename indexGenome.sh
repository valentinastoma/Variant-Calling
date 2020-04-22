#!/usr/bin/env bash
# indexGenome.sh

bwa index -a bwtsw GRCh38_reference.fa \
1>indexGenomeBwa.log 2>indexGenomeBwa.err

# Do not need to run on Defiance
# use algorithm bwt-sw
# symbolic links here

# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa.sa GRCh38_reference.fa.sa
# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa.amb GRCh38_reference.fa.amb
# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa.ann GRCh38_reference.fa.ann
# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa.pac GRCh38_reference.fa.pac
# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa.bwt GRCh38_reference.fa.bwt

samtools faidx GRCh38_reference.fa \
1>indexGenomeSamtools.log 2>indexGenomeSamtools.err

# Many NGS tool require *.fai files, so just create one
# Do not need to run on Defiance
# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa.fai GRCh38_reference.fa.fai