#!/bin/bash
#Pipeline for mapping reads and calling SNPs from RNAseq data

#Index assembly bowtie2
bowtie2-build ahy.fa ahy

#Map reads to assembly
for i in *.fq.gz; do bowtie2 --rg-id $(basename $i .fq.gz) --rg SM:$(basename $i .fq.gz) --very-sensitive -x ahy -U $i > $(basename $i .fq.gz).sam; done

#Convert .sam to .bam
for i in *.sam; do samtools view -bSq 10 $i > $(basename $i .sam)_UNSORTED.bam; done

#Sort and index your alignments
for i in *UNSORTED.bam; do samtools sort $i > $(basename $i _UNSORTED.bam).bam; samtools index $(basename $i _UNSORTED.bam).bam; done

#Remove intermediate files
rm *UNSORTED.bam
rm *.sam

#Index assembly for freebayes
samtools faidx ahy.fa

#Call SNPS
freebayes --genotype-qualities -f ahy.fa *.bam > ahy_unfiltered.vcf

#Filter SNPs
vcffilter -f "TYPE = snp & QUAL > 30 & AF > 0.05 & AF < 0.95" -g "GQ > 20" ahy_unfiltered.vcf \
| vcfallelicprimitives \
| vcfbiallelic \ 
| vcfnulldotslashdot \
| grep -vF './.' | grep -vF '.|.' \
> biallelic_snps_noNA_minmaf05.vcf

#Create 0,1,2 genotype matrix
vcftools --vcf biallelic_snps_noNA_minmaf05.vcf --012 --out ahy_snps



