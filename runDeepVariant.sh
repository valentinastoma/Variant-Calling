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