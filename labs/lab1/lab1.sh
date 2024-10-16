set -uex

REF=$(refgenie seek hg19/bwa_index)
REF="/home/aelia/data/genome/hg19.fa"

C1=~/labs-cse185/lab1/NA12878_child_1.fq
C2=~/labs-cse185/lab1/NA12878_child_2.fq
F1=~/labs-cse185/lab1/NA12891_father_1.fq
F2=~/labs-cse185/lab1/NA12891_father_2.fq
M1=~/labs-cse185/lab1/NA12892_mother_1.fq
M2=~/labs-cse185/lab1/NA12892_mother_2.fq

ref_path=$(echo ${REF} | sed 's|^/[^/]*||' | sed 's|^/[^/]*||')

mkdir -p alignment/
mkdir -p out/

#bwa index ${REF}
#samtools faidx ${REF}

bwa mem ${REF} ${C1} ${C2} > alignment/NA12878_child.sam
samtools view -S -b alignment/NA12878_child.sam > alignment/NA12878_child.bam
samtools sort alignment/NA12878_child.bam > alignment/NA12878_child.sorted.bam
samtools index alignment/NA12878_child.sorted.bam
samtools flagstat alignment/NA12878_child.sorted.bam

bwa mem ${REF} ${F1} ${F2} > alignment/NA12878_father.sam
samtools view -S -b alignment/NA12878_father.sam > alignment/NA12878_father.bam
samtools sort alignment/NA12878_father.bam > alignment/NA12878_father.sorted.bam
samtools index alignment/NA12878_father.sorted.bam
samtools flagstat alignment/NA12878_father.sorted.bam

bwa mem ${REF} ${M1} ${M2} > alignment/NA12878_mother.sam
samtools view -S -b alignment/NA12878_mother.sam > alignment/NA12878_mother.bam
samtools sort alignment/NA12878_mother.bam > alignment/NA12878_mother.sorted.bam
samtools index alignment/NA12878_mother.sorted.bam
samtools flagstat alignment/NA12878_mother.sorted.bam

samtools view alignment/NA12878_child.sorted.bam chr6:128405804-128605805 | cut -f 9 > alignment/child_template_lengths.txt



#samtools mpileup -r chr6:128405804-128605805 -f ${REF} alignment/*.sorted.bam > alignment/trio.mpileup
ls alignment/*.sorted.bam > alignment/bamfiles.txt
bcftools mpileup --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP' -O v -r chr6:128405804-128605805 -f ${REF} -b alignment/bamfiles.txt > alignment/trio.mpileup
java -jar VarScan.jar mpileup2snp alignment/trio.mpileup --variants --min-var-frequency 0.2 --min-freq-for-hom 0.8 --p-value 0.01 --output-vcf > out/trio2.vcf


#| bcftools call --ploidy 2 --annotate 'FORMAT/GQ' -mv -O v | bcftools norm -f ${REF} -d all -O z -o trio.vcf 

python plottemplatedist.py 
python plotcoveragedist.py