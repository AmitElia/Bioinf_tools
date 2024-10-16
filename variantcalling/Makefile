#
# A variant calling Makefile made by @AmitElia. inspired by biostar handbook.
#

# Set SRA number
SRR ?= SRR1553444

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

SNPEFF = ~/micromamba/envs/bioinfo/share/snpeff-5.0-1/snpEff.config
EXISTS = false

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
BED = dat/bam/${SRR}_${ACC}.bedgraph
BW = dat/bam/${SRR}_${ACC}.bw

# The resulting variant file.
VCF = out/vcf/${SRR}_${ACC}.vcf.gz

# Merged VCF file
ALL_VCF = out/merged/all_variants_${ACC}_${PRJN}.vcf.gz
ALL_VCF_ANN = out/merged/all_variants_${ACC}_${PRJN}_annotated.vcf
ANN_STATS = out/merged/snpEff_summary_${ACC}_${PRJN}.html

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
	@echo "# make init run ACC=${ACC} SRR=${SRR} N=${N}"
	@echo "#"
	@echo "# make init runbatch ACC=${ACC} N=${N} PRJN=${PRJN} NR=${NR}"
	@echo "#"
	@echo "# make annotate ACC=${ACC} PRJN=${PRJN}"
	@echo "#"
	@echo "# make visualize ACC=${ACC} SRR=${SRR}"
	@echo "#"
	@echo "# init - initialize the genome."
	@echo "#"
	@echo "# run - get fastq, trim, align, and call variants."
	@echo "#"
	@echo "# runbatch - fastq, trim, align, and call variants based on PRJN bioproject number."
	@echo "#"
	@echo "# annotate - creates custom annotation based on ACC number."
	@echo "#"
	@echo "# visualize - visualize alignment using IGV"
	@echo "#"
	@echo "# clean - remove all files created."
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

# Merge all VCF files
merge:
	@echo "# Merge all VCF files into ${ALL_VCF}"
	mkdir -p $(dir ${ALL_VCF})
	bcftools merge -o ${ALL_VCF} vcf/*.vcf.gz
	bcftools index -t ${ALL_VCF}


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

runbatch: runinfo
	cat ${RIDS} | parallel --lb make run ACC=${ACC}  SRR={}
	mkdir -p $(dir ${ALL_VCF})
	bcftools merge -o ${ALL_VCF} $(dir ${VCF})*.vcf.gz
	bcftools index -t ${ALL_VCF}

annotate:
	mkdir -p $(dir ${SNPEFF})data/${ACC}
	bio fetch ${ACC} > $(dir ${SNPEFF})data/${ACC}/genes.gbk
	if grep -q ${ACC} ${SNPEFF}; \
	then \
		echo "exists in database"; \
	else \
		echo "${ACC}.genome : ${ACC}" >> ${SNPEFF}; \
		echo "${ACC}.reference : ${ACC}" >> ${SNPEFF}; \
		echo "${ACC}.chromosomes : ${ACC}" >> ${SNPEFF}; \
	fi
	snpEff build -genbank -v ${ACC}
	snpEff dump ${ACC}
	snpEff ${ACC} ${ALL_VCF} -stats ${ANN_STATS} > ${ALL_VCF_ANN}


visualize:
#Visualize
	bedtools genomecov -ibam ${BAM} -g ${REF}.fai -bg > ${BED}
	bedGraphToBigWig ${BED} ${REF}.fai ${BW}
	echo -e "new\ngenome ${REF}\nload ${BAM}\nload ${BED}\nload ${BW}\nload ${GFF}\nexpand" >> igv_commands.txt

	bash ${HOME}/IGV/igv.sh -b igv_commands.txt -p 60151

# Clean up all results.
clean:
	rm -rf dat out
	
	
# Inform make that these targets are not files.
.PHONY: genome index fastq trim align vcf run runinfo runbatch annotate visualize merge clean