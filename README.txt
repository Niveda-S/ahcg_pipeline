
Download Reference Files
	wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz
This contains hg19 and dbsnp vcf for hg19

Bowtie indexes
	bowtie2-build -f hg19.fa
Assigns reference as hg19.fa 

Install samtools
	sudo apt-get install samtools

Install Java
	sudo apt-get install openjdk-8-jre

Create index file using samtools
###The index file must be in the same location as hg19.fa###
samtools faidx resources/genome/hg19.fa

Create dict file using Picard
###The dict file must be in the same location as hg19.fa###
jre1.8.0_101/bin/java -jar lib/picard.jar CreateSequenceDictionary R=resources/genome/hg19.fa O=hg19.dict


Test Files

wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq

Run the script
##Gunzip dbsnp file##
##In script, change the path to Java##

python ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i Test/test_r1.fastq 
Test/test_r2.fastq -w Bowtie_index/hg19 -d resources/dbsnp/dbsnp_138.hg19.vcf -r resources/genome/hg19.fa -a lib/Trimmomatic-0.36/adapters/TruSeq3-SE.fa -o ./

gitignore



