# ahcg_pipeline
Variant calling pipeline for genomic data analysis

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```

##Workflow##
###Download Reference Files###
wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz
This contains hg19 and dbsnp vcf for hg19

###Bowtie indexes###
bowtie2-build -f hg19.fa
Assigns reference as hg19.fa 

###Install samtools###
sudo apt-get install samtools

###Install Java###
sudo apt-get install openjdk-8-jre

###Create index file using samtools###
samtools faidx resources/genome/hg19.fa

###Create dict file using Picard###
jre1.8.0_101/bin/java -jar lib/picard.jar CreateSequenceDictionary R=resources/genome/hg19.fa O=hg19.dict

###Test Files###
```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
head -100000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq
head -100000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq
```

###Run the script###
```{sh}
python ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i Test/test_r1.fastq Test/test_r2.fastq -w Bowtie_index/hg19 -d resources/dbsnp/dbsnp_138.hg19.vcf -r resources/genome/hg19.fa -a lib/Trimmomatic-0.36/adapters/TruSeq3-SE.fa -o ./
```

##Mapping regions of interest for BRCA1##
Steps for extracting reads mapping to BRCA1 from NA12878 HiSeq Exome dataset:

1. Downloading the NA12878 HiSeq Exome dataset:
     The bam files for the sample can be downloaded from the ftp link mentioned on GIAB GitHub page.
     There are four runs for this sample, we can start by downloading one of them and extracting the reads.

2. Using samtools to subset the bam file to regions corresponding to BRCA1:
     Using the bed file containing the BRCA1 exonic coordinates we can subset the NA12878 sample using samtools

     samtools view -L <bed file> -b -o < outout bam file > < input bam file >

     Note: -b just specifies that the output needs to be a bam file.

3. Using bedtools to convert the bam file to a fastq file:
     From the brca1 bam file we now extract the reads aligning to the region using bedtools

     bedtools bamtofastq -i <bam file> -fq < fastq r1> -fq2 < fastq r2>

## Variant Quality Score Recalibration (VQSR)##
```{sh}
jre1.8.0_101/bin/java -Xmx4g -jar lib/GenomeAnalysisTK.jar 
-T VariantRecalibrator 
-R resources/genome/hg19.fa 
-input NA12878_variants.vcf 
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf 
-resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf 
-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf 
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 resources/dbsnp/dbsnp_138.hg19.vcf 
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP 
-recalFile output.recal -tranchesFile output.tranches -rscriptFile output.plots.R
```

## DCM ##
```{sh}
#shrink clinvar to just DCM genes
bedtools intersect -a clinvar.vcf.gz -b dcm_gene_list.bed -header > clinvar_allfrombed.vcf

#shrink variants to just DCM genes
bedtools intersect -a patient2_variants_recal.vcf -b dcm_gene_list.bed -header > patient2_dcm_final.vcf

#match variants to clinvar
bedtools intersect -b patient2_dcm_final.vcf -a clinvar_allfrombed.vcf -header > patient2_intersect_clinvar.vcf

#generate simple report on findings
python3 parse_clnsig.py -i patient2_intersect_clinvar.vcf.gz 2>&1 | tee patient2_simple_report.txt
cut -c 24- patient2_simple_report.txt
```

