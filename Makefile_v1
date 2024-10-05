#
# A variant calling Makefile made by @AmitElia. inspired by biostar handbook.
#

# Set SRA number
SRR ?= SRR1972739

# The accession number.
ACC ?= AF086833

# Set the project number.
PRJN ?= PRJNA257197

# Number of reads to get
N ?= 10000
NR ?= 6

#runinfo file.
RINF = dat/runinfo/${PRJN}_runinfo.csv
RIDS = dat/runinfo/${PRJN}_ids.txt
LINES := $(shell cat ${RIDS})
# Reference genome.
REF = dat/refs/${ACC}.fa

# The index file for the reference
IDX = ${REF}.amb

# The GFF file
GFF = dat/refs/${ACC}.gff

# Downloaded reads.
R1 = dat/reads/${SRR}_1.fastq
R2 = dat/reads/${SRR}_2.fastq

RT1 = dat/qc/${SRR}_1.trim.fq
RT2 = dat/qc/${SRR}_2.trim.fq

# The resulting alignment file.
BAM = dat/bam/${SRR}_${ACC}.bam

# The resulting variant file.
VCF = out/vcf/${SRR}_${ACC}.vcf.gz

# Merged VCF file
ALL_VCF = out/merged/all_variants_${ACC}_${PRJN}.vcf.gz

# -- Nothing needs to change below --
SHELL := bash
.SHELLFLAGS := -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

help:
    
	@echo "#"
	@echo "# Usage:"
	@echo "#"
	@echo "# ACC=${ACC} SRR=${SRR} N=${N} PRJN=${PRJN} NR=${NR}"
	@echo "#"
	@echo "# make init - initialize the genome"
	@echo "#"
	@echo "# make run - get fastq, trim, align, and call variants"
	@echo "#"
	@echo "# make runbatch - fastq, trim, align, and call variants based on PRJN bioproject number"
	@echo "#"
	@echo "# make runlist - fastq, trim, align, and call variants based on a runlist of SRRs"
	@echo "#"
	@echo "# make clean - remove all files created"
	@echo "#"


# Downloads the reference genome.
${REF}:
	@echo "# Fetching FASTA file for ACC=${ACC}"
	mkdir -p $(dir ${REF})
	bio fetch ${ACC} -format fasta > ${REF}
${GFF}:
	@echo "# Fetching GFF file for ACC=${ACC}"
	mkdir -p $(dir ${REF})
	bio fetch ${ACC} -format gff > ${GFF}

# Triggers the genome download.
genome: ${REF} ${GFF}
	@ls -lh ${REF}
	@ls -lh ${GFF}
	
# Index the genome with bwa and samtools.
${REF}.bwt: ${REF}
	@echo "# Build bwa index of the reference genome ${REF} for the aligner."
	bwa index ${REF}
	@echo "# Build index of the reference genome ${REF} for IGV."
	samtools faidx ${REF}

# Trigger the indexing.
index : ${REF}.bwt
	@ls -lh ${REF}.bwt

# Download the FASTQ files.
${R1}:
	mkdir -p $(dir ${R1})

	@echo "# Obtain the FASTQ sequences for the SRR=${SRR}"
	fastq-dump -X ${N} --outdir $(dir ${R1}) --split-files ${SRR} 

# Trigger the fastq download.
fastq: ${R1}
	@ls -lh ${R1}
	fastqc ${R1} ${R2}

#QC and trimming reads
trim:
	mkdir -p $(dir ${RT1})
	fastp --cut_tail -i ${R1} -I ${R2} -o ${RT1} -O ${RT2}
	fastqc ${RT1} ${RT2}

# Align the reads.
${BAM}:
	mkdir -p $(dir ${BAM})
	bwa mem ${REF} ${RT1} ${RT2} | samtools sort > ${BAM}
	samtools index ${BAM}
	samtools flagstat ${BAM}

# Trigger the alignment.
align: ${BAM} 
	@ls -lh ${BAM}

align!: ${BAM}
	rm -f ${BAM}

# Compute the VCF file.
${VCF}: ${BAM}
	# Create the directory for the VCF file.
	mkdir -p $(dir ${VCF})

	# Add more annotations than the default. Normalize the variants.
	bcftools mpileup \
		--annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' \
		-O v -f ${REF} ${BAM} | \
		bcftools call --ploidy 2 --annotate 'FORMAT/GQ' -mv -O v | \
		bcftools norm -f ${REF} -d all -O z -o ${VCF} 
	
	# Index the VCF file.
	bcftools index -t ${VCF}

# Trigger the VCF file generation.
vcf: ${VCF}
	@ls -lh ${VCF}


# Initialize the genome and index.
init: genome index

# Runs the data related pipeline.
run: fastq trim align vcf

# Undo the run results.
run!:
	rm -rf ${BAM} ${VCF} ${R1} ${R2} ${RT1} ${RT2}

#Download runinfo.csv based on PRJN number
${RINF}:

	mkdir -p $(dir ${RINF})

	@echo "# Fetching runinfo.csv for PRJN=${PRJN}"
	esearch -db sra -query ${PRJN} | efetch -format runinfo > ${RINF}

#extract SRR list from runinfo.csv
${RIDS}: 
	cat ${RINF} | csvcut -c Run | grep -v Run | tail -${NR} > ${RIDS}

#trigger
runinfo: ${RINF} ${RIDS}
	@ls -lh ${RINF}
	@ls -lh ${RIDS}

fastqbatch:

	mkdir -p $(dir ${R1})

	#download data from SRA
	cat ${RIDS} | parallel --progress -j 3 fastq-dump -O $(dir ${R1}) -X ${N} --split-files {}

	#run fastqc
	cat ${RIDS} | parallel --progress fastqc $(dir ${R1})/{}_1.fastq $(dir ${R1})/{}_2.fastq

trimbatch:

	mkdir -p $(dir ${RT1})

	#apply QC and trim
	cat ${RIDS} | parallel --progress fastp --cut_tail -i $(dir ${R1}){}_1.fastq \
	-I $(dir ${R1}){}_2.fastq -o $(dir ${RT1}){}_1.trim.fq -O $(dir ${RT1}){}_2.trim.fq

	#run fastqc on trimmed reads
	cat ${RIDS} | parallel --progress fastqc $(dir ${RT1}){}_1.trim.fq $(dir ${RT1}){}_2.trim.fq

alignbatch:

	mkdir -p $(dir ${BAM})

	#align reads to reference genome
	
	for line in $$(cat ${RIDS}); \
	do \
		echo "Processing $$line"; \
		bwa mem ${REF} $(dir ${RT1})$${line}_1.trim.fq $(dir ${RT1})$${line}_2.trim.fq | \
		samtools sort > $(dir ${BAM})$${line}_${ACC}_${PRJN}.bam; \
		samtools index $(dir ${BAM})$${line}_${ACC}_${PRJN}.bam; \
		samtools flagstat $(dir ${BAM})$${line}_${ACC}_${PRJN}.bam; \
	done

vcfbatch:
	mkdir -p $(dir ${ALL_VCF})


	bcftools mpileup \
		--annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' \
		-O v -f ${REF} $(dir ${BAM})*_${ACC}_${PRJN}.bam | \
		bcftools call --ploidy 2 --annotate 'FORMAT/GQ' -mv -O v | \
		bcftools norm -f ${REF} -d all -O z -o ${ALL_VCF} 

runbatch: runinfo fastqbatch trimbatch alignbatch vcfbatch

# Clean up all results.
clean:
	rm -rf dat out


# Inform make that these targets are not files.
.PHONY: genome index fastq trim align vcf run runinfo fastqbatch trimbatch vcfbatch merge clean