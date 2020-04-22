## Introduction

With the rapid development of Whole Genome Sequencing and Whole Exome
sequencing techniques, the bottleneck of research, according to Chen et
al. (2019), has become the ability to accurately call complex mutations.
Utilization of variant calling has allowed us to detect a number of
single-nucleotide polymorphisms, insertions and deletions. Accurate
description of these difference between the genomes allows us to
characterize cells and tissue, depending on whether the samples are
germline (belonging to the same species) or somatic (different tissues
of the same individual). For this exercise, the pipeline repeats steps
taken by the researchers in a comparison description of DeepVariant
accuracy as a variant calling tool (Supernat et al. (2018)).

## Methods

In order to perform variant calling on human genome sequencing, the
human genome sequencing file needs to be downloaded from the internet
which can be obtained with wget from the open internet source
(sequencing database). The reads NGS read which will be assessed in this
pipeline are also obtained from an open source under SRA accession
number SRR6808334. Fastq-dump is used to covert the source files which
are in the .sra format into either fastq or fasta; split-files is used
to dump each read into a separate file and the names of the files will
receive the corresponding suffix.

Once the genome sequencing is obtained, an index can be created using
Burrows – Wheelers alignment tool (H. Li and Durbin (2009)). BWA has
three algorithms for mapping – backtrack, SW, and MEM. Each can be
chosen depending on the length of the reads one is trying to map (
backtrack for shorter sequences up to 100bp, SW and MEM for longer –
70bp- 1Mbp). The genome is indexed using bwtsw algorithm (SW) because of
the length of the mapping done for human genome. The program is
inititated by calling “bwa idex”. The setting -a bwtsw is used to call
the algorythm. Input file needs to be specified, as well as the outputs
for the index log and error.

Reads obtained from SRA need to be quality trimmed. To achieve this, we
use Trimmomatic (Bolger, Lohse, and Usadel (2014)). Using this trimmer
allows to minimize the presence and subsequent bias from low quality
reads, which are a result of normal sequencing procedure (adapter and
PCR fragments). Trimmomatic works well with paired reads because of its
palindrome mode - it is efficient at removing technical sequence and
achieves better quality reads filtering for paired reads without loosing
this relationship between the paired reads. The program is initiated
through calling the application with -n 19 java -jar, specifiying the
input files for reference, and the paired and unpaired reads. Other
specific conditions are set as follows: headcrop is set at 0, leadign
and trailing options are used to remove bases from the beginnign and end
of the sequences if below a certain threshold. This needs to be done
depending on the read length to minimize the low quality reads that are
often found at the end of the sequences. Sliding window size is used to
clip the reads once the quality of the window falls below a certain
threshold - the values used for this option in this case are
SLIDINGWINDOW:4:30. MINLEN is :36 which indicates that the read will be
dropped if it is below this length in bases.

Once the reads are trimmed and the genome is indexed, we can proceed to
align reads using one of the BWA algorithms – MEM (H. Li and Durbin
(2009)). Using maximal exact matches is an efficient and accurate way to
achieve the target alignment. For the number of threads, we indicate 8;
additionally, we utilize -M function which allows to mark short reads as
secondary. We create a read group header line with -R option – this
allows to add an ID and group sample (SM). The read group will be added
to the output alignment .sam files.

Obtained aligned files need to be sorted into binary .sam = .bam format.
That can be achieved using samtool (H. Li et al. (2009)) sort function.
We set the number of sorting and comprehending threads to 8, as well as
out a limit to the maximum possible memory usage for this sorting per
thread – 4G. Next we identify the input .sam file and specify the name
of the output .bam file. Upon creating a binary .sam file, we index the
reads using another function of samtools (H. Li et al. (2009)). Similar
specification of 8 threads and memory limitation to 4GB are set for the
input file of .bam format. Indexing algorithn of samtools if called by
“samtools index” command, threads are indicated through ‘@’ and -bn
option which indicates that we want to create a BAI index as a format.

Lastly, we perform variant calling. In this case, DeepVariant is used
(Supernat et al. (2018)). The program utilizes deep neural network to
call genetic variants. The process can be simplified into three stages –
preprocessing, calling variants, and post-processing. In order to run
the program ,we need to specify the inpute file, directory, reference
file, bam file, number of references to be made in the pre – processing
step, which works through creating “examples” using an internal
TensorFlow format. The number of the examples needs to be specified as
the number of CPU cores one has. Next, specify the output directory,
file name for VCF, and output for the logs. If docker is not found in
PATH, it needs to be downloaded, which is what is done in the example
script. Once the docker is ready to use, DeepVariant is run – input and
ouput directories are previously specified and referenced here,
model\_type needs to be indicated as WGS in this case, reference refers
to the whole genome file used as a reference earlier, reads are
specified as bam file, previously output directory and file name are
references, as well as the number of shards to be used.

Including the code for easier review and reference of specific outputs:

``` bash
#!/usr/bin/env bash
# runDeepVariant.sh

# DeepVariant is an analysis pipeline that uses a deep neural network 
# to call genetic variants from next-generation DNA sequencing data

# if you want to run this on google cloud see the following:
# https://cloud.google.com/life-sciences/docs/tutorials/deepvariant

# code here:
# https://github.com/google/deepvariant


# Running DeepVariant consists of three stages:

# Making examples: DeepVariant pre-processes the input data and saves examples 
# from the data using an internal TensorFlow format. 
# You can run this stage in parallel where all of the input shards are processed independently.

# Calling variants: DeepVariant runs a deep neural network that makes inferences 
# from the examples and saves them into shared files using an internal TensorFlow format.

# Post-processing variants: DeepVariant converts variants from the internal TensorFlow format to VCF 
# or gVCF files. This stage runs on a single thread.

# some slight changes have been made, but pretty much what Google above says to do:

set -euo pipefail

BIN_VERSION="0.9.0"
BASE="/mnt/disks/sdb/binf6309/VariantCalling"
INPUT_DIR="${BASE}/input/data"
REF="GRCh38_reference.fa"
BAM="SRR6808334.sorted.bam"
# How many cores the `make_examples` step uses. Change it to the number of CPU cores you have
N_SHARDS="64"
OUTPUT_DIR="${BASE}/output"
OUTPUT_VCF="SRR6808334.output.vcf.gz"
OUTPUT_GVCF="SRR6808334.output.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

## Create directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"
mkdir -p "${LOG_DIR}"

## Downloads
sudo apt-get -qq -y update

if ! hash docker 2>/dev/null; then
      echo "'docker' was not found in PATH. Installing docker..."
      # Install docker using instructions on:
      # https://docs.docker.com/install/linux/docker-ce/ubuntu/
      sudo apt-get -qq -y install \
          apt-transport-https \
          ca-certificates \
          curl \
          gnupg-agent \
          software-properties-common
      curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
      sudo add-apt-repository \
          "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
          $(lsb_release -cs) \
          stable"
      sudo apt-get -qq -y update
      sudo apt-get -qq -y install docker-ce
fi

# Copy the data
echo "Copying data"
cp SRR6808334.sorted.bam -d "${INPUT_DIR}"
cp SRR6808334.sorted.bam.bai -d "${INPUT_DIR}"
cp GRCh38_reference.fa -d "${INPUT_DIR}"
cp GRCh38_reference.fa.fai -d "${INPUT_DIR}"
#cp GRCh38_reference.fa.gz -d "${INPUT_DIR}"
#cp GRCh38_reference.fa.gz.gzi -d "${INPUT_DIR}"
#cp GRCh38_reference.fa.gz.fai -d "${INPUT_DIR}"

## Pull the docker image.
# https://console.cloud.google.com/gcr/images/deepvariant-docker/GLOBAL/deepvariant?gcrImageListsize=30
sudo docker pull gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"

# --model_type=WGS \ **Replace this string with exactly one of the following [WGS,WES,PACBIO]**
echo "Running DeepVariant..."
sudo docker run \
      -v "${INPUT_DIR}":"/input" \
      -v "${OUTPUT_DIR}:/output" \
      gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --ref="/input/${REF}" \
      --reads="/input/${BAM}" \
      --output_vcf=/output/${OUTPUT_VCF} \
      --output_gvcf=/output/${OUTPUT_GVCF} \
      --num_shards=${N_SHARDS}
echo "Done."
echo
```

## References:

<div id="refs" class="references">

<div id="ref-Bolger">

Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. 2014. “Trimmomatic: A
Flexible Trimmer for Illumina Sequence Data.” *Bioinformatics (Oxford,
England)* 30 (15): 2114–20.
<https://doi.org/10.1093/bioinformatics/btu170>.

</div>

<div id="ref-Chen">

Chen, Jiayun, Xingsong Li, Hongbin Zhong, Yuhuan Meng, and Hongli Du.
2019. “Systematic Comparison of Germline Variant Calling Pipelines Cross
Multiple Next-Generation Sequencers.” *Scientific Reports* 9 (1): 9345.
<https://doi.org/10.1038/s41598-019-45835-3>.

</div>

<div id="ref-Li">

Li, Heng, and Richard Durbin. 2009. “Fast and Accurate Short Read
Alignment with Burrows-Wheeler Transform.” *Bioinformatics (Oxford,
England)* 25 (14): 1754–60.
<https://doi.org/10.1093/bioinformatics/btp324>.

</div>

<div id="ref-HengLi">

Li, Heng, Bob Handsaker, Alec Wysoker, Tim Fennell, Jue Ruan, Nils
Homer, Gabor Marth, Goncalo Abecasis, Richard Durbin, and 1000 Genome
Project Data Processing Subgroup. 2009. “The Sequence Alignment/Map
Format and SAMtools.” *Bioinformatics (Oxford, England)* 25 (16):
2078–9. <https://doi.org/10.1093/bioinformatics/btp352>.

</div>

<div id="ref-Supernat">

Supernat, Anna, Oskar Valdimar Vidarsson, Vidar M. Steen, and Tomasz
Stokowy. 2018. “Comparison of Three Variant Callers for Human Whole
Genome Sequencing.” *Scientific Reports* 8 (1): 17851.
<https://doi.org/10.1038/s41598-018-36177-7>.

</div>

</div>
