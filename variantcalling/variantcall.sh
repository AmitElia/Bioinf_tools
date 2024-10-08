# Variant calling script using bwa aligner and generates genotypes using bcftools, freebayes, and GATK. Visualizing using IGV.

# Usage: bash variantcall.sh [WD] [ACC] [SRR]

set -uex

# get working dir
WD=$1

#set data dir
DIR=${WD}/data/

# Ref genome accession number for Ebola virus - Mayinga, Zaire, 1976, complete genome.
ACC=$2

# Sequencing reads number
SRR=$3

# file path for fasta genome
REF=${DIR}refs/${ACC}.fa

# file path for BAM file
BAM=${DIR}alignments/${SRR}.bam

# file path for fastq files
SEQ=${DIR}sequences/

# shortcut for sequence names
R1=${SEQ}${SRR}_1.fastq
R2=${SEQ}${SRR}_2.fastq

# Tag the alignments. GATK will only work when reads are tagged.
TAG="@RG\tID:$SRR\tSM:$SRR\tLB:$SRR"

# Align and generate a BAM file. example using BWA
bwa mem -R $TAG $REF $R1 $R2 | samtools sort > $BAM

# Index the BAM file
samtools index $BAM

# Compute the genotypes from the alignment file | Call the variants with bcftools.
bcftools mpileup -O v -f $REF $BAM | bcftools call --ploidy 1 -vm -O v > ${WD}/variants_bcftools.vcf

# Compute the genotypes from the alignment file + Call the variants with freebayes.
freebayes -f $REF $BAM > ${WD}/variants_fb.vcf

# Compute the genotypes from the alignment file + Call the variants with GATK.

# Generate Dictionary index.
picard CreateSequenceDictionary REFERENCE=$REF  OUTPUT=${DIR}refs/$ACC.dict

# Generate the variants.
gatk HaplotypeCaller -R $REF -I $BAM -O ${WD}/variants_gatk.vcf