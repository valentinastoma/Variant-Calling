#!/usr/bin/env bash
# getGenome.sh

#!/usr/bin/env bash
# getGenome.sh

# Get the GRCh38 reference genome

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/GRCh38.primary_assembly.genome.fa.gz \
    -O GRCh38_reference.fa.gz \
    1>getGenome.log 2>getGenome.err

gunzip GRCh38_reference.fa.gz

# Do not download on Defiance!!
# you can just soft link this to your working directory on defiance
# symbolic links here

# ln -s /data/METHODS/Spring/Module05/GRCh38_reference.fa GRCh38_reference.fa