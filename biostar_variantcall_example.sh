# Example of Variant calling from the biostar handbook.
# Link to reproduced paper: https://www.science.org/doi/full/10.1126/science.1259657
# SRA BioProject: PRJNA257197

set -uex

# Get current time
#currTime=$(date +"%Y_%m_%d.%H_%M_%S")
currTime=$1

# Set working directory
WD=variantcall/${currTime}

# Ref genome accession number for Ebola virus - Mayinga, Zaire, 1976, complete genome.
ACC=AF086833
#Sequencing reads number
SRR=SRR1553500

#call getdat.sh to get data
#Usage: bash getdat.sh [WD] [ACC] [SRR]
bash getdat.sh ${WD} ${ACC} ${SRR}

#call variantcall.sh to run variantcall
#Usage: bash variantcall.sh [WD] [ACC] [SRR]
bash variantcall.sh ${WD} ${ACC} ${SRR}

DIR=${WD}/data/

# file path for fasta genome
REF=${HOME}/scripts/${DIR}refs/${ACC}.fa


# file path for BAM file
BAM=${HOME}/scripts/${DIR}alignments/${SRR}.bam

# file path for gff files
GFF=${HOME}/scripts/${DIR}refs/${ACC}.gff


#Visualize

echo new >> igv_commands.txt
echo genome ${REF} >> igv_commands.txt
echo load ${BAM} >> igv_commands.txt
echo load ${GFF} >> igv_commands.txt
echo expand >> igv_commands.txt

#cat igv_commands.txt | nc localhost 60151

bash ${HOME}/IGV/igv.sh -b igv_commands.txt -p 60151