##UC Genomics Workshop @ Asilomar: How to call SNPs from transcriptomic data
by the Palumbi Lab, September 2016

About the dataset: 
Reads from a subset of contigs from 38 Acropora hyacinthus samples, 100bp single end reads, from Bay and Palumbi (2014) Current Biology

Important: you need Command Line Tools installed to install and run these programs

```
#are command line tools installed? If so, this command will return where it is installed
xcode-select -p
#if command line tools are not installed, install them
xcode-select --install
```

Programs needed: bowtie2, samtools, freebayes, vcflib, VCFtools, R

Easy install to try pipeline on a personal mac:

- install homebrew, a science package manager

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" 
brew tap homebrew/science

# install packages used for this pipeline
brew install bowtie2
brew install samtools
brew install freebayes
brew install vcflib
brew install vcftools

```

- if you get a permissions error, try 
`sudo chown -R "$USER":admin /usr/local"`


##Set up your workspace
- open Terminal
- move to desktop 
`cd ~/Desktop`
- clone this repository
`git clone https://github.com/bethsheets/SNPcalling_tutorial.git`
- move into directory 
`cd SNPcalling_from_RNAseq`


##Step 1: Map reads to assembly with Bowtie2

1a) Make a Bowtie2 index of your assembly

`bowtie2-build ahy.fa ahy`

1b) Check that the index outputs 6 .bt2 files correctly 
ls

1c) Call program:

```
for i in *.fq.gz; do bowtie2 --rg-id $(basename $i .fq.gz) --rg SM:$(basename $i .fq.gz) --very-sensitive -x ahy -U $i > $(basename $i .fq.gz).sam; done
```
- --rg-id & --rg adds sample ids to your alignments, so you can combine them later but still tell which reads go with which samples
- --very-sensitive is running -D 20 -R 3 -N 0 -L 20 -i S,1,0.50.
- -x is your Bowtie2 index that you made previously
- -U is your input file


##Step 2: convert bam to .sam files using Samtools

2a) use ‘view’ to convert SAM<->BAM

`for i in *.sam; do samtools view -bSq 10 $i > $(basename 	$i .sam)_UNSORTED.bam; done`

2b) sort and index your alignments

`for i in *UNSORTED.bam; do samtools sort $i > $(basename $i _UNSORTED.bam).bam; samtools index $(basename $i _UNSORTED.bam).bam; done`

2c) remove intermediate files

```
rm *UNSORTED.bam
rm *.sam
```


##Step 3: Call SNPs with FreeBayes
3a) Index fasta for freebayes

`samtools faidx ahy.fa`

3b) Call program:
`freebayes --genotype-qualities -f ahy.fa *.bam > ahy_unfiltered.vcf`

- --genotype-qualities : Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
- -f : reference assembly

##Step 4: Filter SNPs with VCFlib

4a) We filter for:

- high quality SNPs (99.9% confident of SNP site, 99% confident of individual genotype) with minimum allele frequency of 5%
- strip away complex extra haplotype information present in some snps 
- discard multi-allelic snps 
- mark snps with missing genotypes
- remove snps with missing genotypes
- pass the filtered information into final vcf

```
vcffilter -f "TYPE = snp & QUAL > 30 & AF > 0.05 & AF < 0.95" -g "GQ > 20" ahy_unfiltered.vcf \
| vcfallelicprimitives \
| vcfbiallelic \ 
| vcfnulldotslashdot \
| grep -vF './.' | grep -vF '.|.' \
> biallelic_snps_noNA_minmaf05.vcf
```

##Step 5: Create 0,1,2 genotype SNP matrix with VCFtools

`vcftools --vcf biallelic_snps_noNA_minmaf05.vcf --012 --out ahy_snps`

##Step 6: Format SNP Matrix in R

6a) Open R program

6b) In a new document, paste the following

```
setwd('~/Desktop/ucasilomar')
snps<-read.delim('ahy_snps.012',header=F,na=-1,row.names=1)
pos<-read.delim('ahy_snps.012.pos',header=F)
indv<-read.delim('ahy_snps.012.indv',header=F)

colnames(snps)<-paste(pos[,1],pos[,2],sep=':')
rownames(snps)<-indv[,1]
snps<-as.matrix(snps)

#read in meta data
meta<-read.delim('~/Desktop/ucasilomar/meta.txt')
```


#How to use genotype data:

##PCA
Plot a PCA to visually identify any clusters within your data. Are these clusters associated with your meta data?

In R

```
pc.out<-prcomp(snps)
summary(pc.out)
plot(pc.out$x[,1],pc.out$x[,2],col=meta$Pool)
```

##Fst

In R

```
install.packages('hierfstat')
library(hierfstat)

#prepare 0,1,2 matrix in hierfstat format
#we use our pca to separate samples into 
#clusters to test for genetic differentiation
hf<-snps
hf[hf==0]<-11
hf[hf==1]<-12
hf[hf==2]<-22
pop=as.numeric(pc.out$x[,1]>2)+1
hf<-as.data.frame(cbind(pop,snps))

#calculate Weir-Cockerham Fst
fst.out<-wc(hf)

#global estimate
fst.out$FST

#look at fst distribution across sites
site.fst<-fst.out$per.loc[['FST']]
hist(site.fst)
```

## Other interesting analyses

- structure/admixture: ngsAdmix takes bam files
- Outliers & environmental data: outFLANK, bayenv, etc…
- local & global linkage: vcftools or R linkage package
- somatic mutations: are your samples high depth from the same individual? if so, you could look at this
- dN/dS: orf prediction with biopython scripts or snpEff
- eQTLs: are snps associated with expression?
