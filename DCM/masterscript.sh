#!/bin/bash

BAM_PATH=
PIPELINE_PATH=
TRIM_PATH=
BOWTIE_PATH=
PICARD_PATH=
GATK_PATH=
BOWTIE_INDEX=
DBSNP_PATH=
REF_PATH=
ADAPTER_PATH=
RECALIBRATE_PATH=
GENE_LIST=

while getopts "b:p:t:o:c:g:w:d:r:a:e:n:" opt; do
	case $opt in
		b) BAM_PATH="$OPTARG";;
		p) PIPELINE_PATH="$OPTARG";;
		t) TRIM_PATH="$OPTARG";;
		o) BOWTIE_PATH="$OPTARG";;
		c) PICARD_PATH="$OPTARG";;
		g) GATK_PATH="$OPTARG";;
		w) BOWTIE_INDEX="$OPTARG";;
		d) DBSNP_PATH="$OPTARG";;
		r) REF_PATH="$OPTARG";;
		a) ADAPTER_PATH="$OPTARG";;
		e) RECALIBRATE_PATH="$OPTARG";;
		n) GENE_LIST="$OPTARG";;
	esac
done

#Convert bam to fastq
bedtools bamtofastq -i $BAM_PATH -fq Patient1_1.fq -fq2 Patient1_2.fq

#Run variant calling pipeline
python $PIPELINE_PATH -t $TRIM_PATH -b $BOWTIE_PATH -p $PICARD_PATH -g $GATK_PATH -i Patient1_1.fq Patient1_2.fq -w $BOWTIE_INDEX -d $DBSNP_PATH -r $REF_PATH -a $ADAPTER_PATH -o ./

#Variant Recalibration
../../ahcg_pipeline/jre1.8.0_101/bin/java -Xmx4g -jar $GATK_PATH -T VariantRecalibrator -R $REF_PATH -input variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $RECALIBRATE_PATH/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 $RECALIBRATE_PATH/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 $RECALIBRATE_PATH/1000G_phase1.snps.high_confidence.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP_PATH -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -mode SNP -recalFile recal.vcf -tranchesFile output.tranches -rscriptFile output.plots.R

#Coverage
samtools view -L $GENE_LIST $BAM_PATH -b > new.bam
bedtools genomecov -ibam new.bam -bga > coverage_output.bed
bedtools intersect -loj -split -a $GENE_LIST -b coverage_output.bed > cov.bed
awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n", $1,$2,$3,$8,$4,$6)}' cov.bed > new_cov.bed
python cov.py new_cov.bed new_cov_depth.txt
awk '{print >> $4}' new_cov_depth.txt
genes=( "LMNA" "MYBPC3" "MYH6" "MYH7" "SCNSA" "TNNT2" )
for i in "${genes[@]}"
do
	Rscript draw_depth.R "$i" "$i".png
done

#Intersect
bedtools intersect -a clinvar.vcf -b dcm_gene_list.bed -header > clinvar_allfrombed.vcf
bedtools intersect -a patient1_variants_recal.vcf -b dcm_gene_list.bed -header > patient1_dcm_final.vcf
bedtools intersect -b patient1_dcm_final.vcf -a clinvar_allfrombed.vcf -header > patient1_intersect_clinvar.vcf
python3 parse_clnsig.py -i patient1_intersect_clinvar.vcf 2>&1 | tee patient1_simple_report.txt
cut -c 24- patient1_simple_report.txt

#Create report as PDF
convert patient1_simple_report.txt LMNA.png MYBPC3.png MYH6.png MYH7.png SCNSA.png TNNT2.png report.pdf
