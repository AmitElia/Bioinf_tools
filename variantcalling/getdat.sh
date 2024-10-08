#General script to retrieve data for variant calling

#Usage: bash getdat.sh [WD] [ACC] [SRR]

set -uex

# get working dir
WD=$1

#set data dir
DIR=${WD}/data/

#make directories for dataf
mkdir -p ${DIR}
mkdir -p ${DIR}refs/
mkdir -p ${DIR}alignments/
mkdir -p ${DIR}sequences/

# Ref genome accession number for Ebola virus - Mayinga, Zaire, 1976, complete genome.
ACC=$2

#Sequencing reads number
SRR=$3

#file path for fasta genome
REF=${DIR}refs/${ACC}.fa

#file path for BAM file
BAM=${DIR}alignments/${SRR}.bam

#file path for fastq files
SEQ=${DIR}sequences/

# file path for gff files
GFF=${DIR}refs/${ACC}.gffc

# Obtain sequences
efetch -db nuccore -format fasta -id ${ACC} | seqret -filter -sid ${ACC} > ${REF}
# Obtain sequences
efetch -db nuccore -format gff -id ${ACC} > ${GFF}

# Index reference for the aligner.
bwa index ${REF}

# Index the reference genome for IGV
samtools faidx ${REF}

# Get data from an Ebola sequencing run.
fastq-dump -X 100000 --split-files ${SRR} -O $SEQ