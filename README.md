##UC Genomics Workshop @ Asilomar: How to call SNPs from transcriptomic data
by the Palumbi Lab, September 2016

###About the tutorial and dataset: 

- This Github repository contains:
	-  a subset of an Acropora hyacinthus transcriptome assembly (ahy.fa)
	-  a subset of contigs from 38 A. hya samples (100bp single end) in their raw format (rab-field###.fq.gz)
	-  a meta data file with information about collection location (meta.txt)
	- BASH (Map\_and\_callSNPs\_RNAseq.sh) and R (Basic\_SNP\_analyses.R) scripts with the full pipeline we follow below.
	- this data from Bay & Palumbi (2014) Current Biology 
- The pipeline uses the following programs: bowtie2, samtools, freebayes, vcflib, VCFtools, and R
- Because the tutorial files are small, you can easily run this tutorial on a personal computer. If you're interested in using these scripts for your own data, check out the last section of the tutorial for more information. 


###How to install programs on a personal mac to try the tutorial
1) Check to see if Command Line Tools is installed on your computer. 

```
# Open terminal
# are command line tools installed? Type the command below. If it is installed, it should return the path. 
#Type the following into terminal and hit enter:
xcode-select -p


#if command line tools are not installed, install them
xcode-select --install
``` 

2) Install homebrew, a software package manager

```
# make sure the permissions on your computer are correct
sudo chown -R "$USER":admin /usr/local

# install homebrew 
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)" 

# tell homebrew we want to use their science software package manager
brew tap homebrew/science
```

3) Install packages used for this pipeline with homebrew

```
brew install bowtie2
brew install freebayes
brew install samtools 
brew install vcftools

#if you get error "brew link step did not complete successfully" for vcftools or samtools, type:
brew link vcftools
#or
brew link samtools
```

4) Check to see that each program was installed correctly by calling each of them, for example:

```
bowtie2
#the manual should pop up on your screen 
```

5) Install the R program 

- Visit: https://cloud.r-project.org
- this program has its own console outside of the terminal we used above

### Explanation of bash utilities used in our scripts

- wildcards: `*.txt`
	- you can reference a set of files that have parts of their names in common (for instance all files that end in .txt) using the * character, which refers to any number of characters (excluding things like spaces and tabs). For instance \*.txt refers to all text files in the current directory. 
- for loops: `for i in *.txt; do <command> $i; done`
	- This loops through each text file in the current directory in alphabetical order; each text file is given the temporary name "i" inside the loop, which we can use to carry out commands on several files. The syntax "$i" allows us to reference the current file inside the loop. 
- basenames: `$(basename $i .txt)`
	- This extracts the name of a file without the extension. You can append a new filename extension by adding it after the command, like this: `$(basename $i.txt).fa` This is useful for naming the output of commands on files that you are processing in a for loop.

###Set up your workspace
1) Open Terminal

2) move to your Desktop 

`cd ~/Desktop`

3) clone the SNPCalling_tutorial repository from Github

`git clone https://github.com/bethsheets/SNPcalling_tutorial.git`

4) move into the downloaded directory 

`cd SNPcalling_tutorial`


##Step 1: Map reads to assembly with bowtie2

1a) Make a Bowtie2 index of your assembly.

`bowtie2-build ahy.fa ahy`

- Indexing your assembly creates an inventory of all the places in the genome that short sequences might match, so that when we align reads we know where to look.

1b) Check that the index outputs 6 .bt2 files correctly 

`ls`

1c) We currently use the raw (non-trimmed/clipped) input files (.fq.gz) straight from the sequencer to map to your assembly. This is because newer mapping and SNP-calling software is able to model low quality sequence at the end of reads and take into account the uncertainty; modeling is preferable to throwing away information. 

1d) Call program:

```
for i in *.fq.gz; do
bowtie2 --rg-id $(basename $i .fq.gz) \
--rg SM:$(basename $i .fq.gz) \
--very-sensitive -x ahy -U $i \
> $(basename $i .fq.gz).sam
done
```

- bowtie2 ouputs a .sam file (sequence alignment/map format) for each sample that contains information about where reads aligned and how similar they are to the reference sequence, and also which reads did not align to the assembly. This is a relatively simple, very large text file that you can open up and read yourself, for instance using the command `less <input>.sam`

- The innards of the script above are:
	- --rg-id & --rg adds sample ids to your alignments, so you can combine them later but still tell which reads go with which samples
	- --very-sensitive is running -D 20 -R 3 -N 0 -L 20 -i S,1,0.50
		- -D give up extending after <int> failed extends in a row
		- -R for reads w/ repetitive seeds, try <int> sets of seeds
		- -N max # mismatches in seed alignment; can be 0 or 1
		- -L length of seed substrings; must be >3, <32
		- -i interval between seed substrings w/r/t read length 
	- -x is your Bowtie2 index that you made previously
	- -U is your input file


##Step 2: convert .sam to .bam files, sort, & index them using samtools 

2a) use ‘view’ to convert .sam to .bam

`for i in *.sam; do samtools view -bSq 10 $i > $(basename 	$i .sam)_UNSORTED.bam; done`

- Here, we are converting .sam files to a machine-readable binary format (.bam). This saves space on your computer and also allows us to organize the information in a way where we can process it faster.

2b) sort and index your alignments

```
for i in *UNSORTED.bam; do
samtools sort $i > $(basename $i _UNSORTED.bam).bam
samtools index $(basename $i _UNSORTED.bam).bam
done
```
- Because reads off of a sequencer are in no particular order, the sam and bam files that result from alignment are also in no particular order. If we want to look at all the reads that map to a contig, for instance to look for genetic variants or to count expression levels, we first need to sort reads so that all of the reads that map to a given contig are grouped together in the bam file. Then we index the sorted bam file, which just creates a summary of which reads go with which contigs and how many reads are associated with each contig; many programs use this summary to guide sequence processing.

2c) remove intermediate files

```
rm *UNSORTED.bam
rm *.sam
```

- We don't need these intermediate files, since we have all of the mapped reads in our sorted bam files. We will hold on to the original fq.gz files in case we ever need them to carry out a different pipeline.


##Step 3: Call SNPs with freebayes
3a) Index your assembly for freebayes

`samtools faidx ahy.fa`

- We make a list of all contigs and how long they are so freebayes knows where to look for SNPs

3b) Call program:

`freebayes --genotype-qualities -f ahy.fa *.bam > ahy_unfiltered.vcf`

- --genotype-qualities : Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output
- -f : reference assembly

##Step 4: Filter SNPs with VCFlib

4a) We filter for/to:

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

- 012 format is a simple way to represent SNP information that can be used in many different analyses. In this format, we have a row for each sample and a column for each SNP. The integer in each entry represents how many copies of the alternate (i.e. non-reference) allele are present in each sample (0: both copies reference allele, 1: heterozygous, 2: both copies are the alternate allele). 

- The only drawback of this format is that it does not represent the uncertainty we have about genotype. For instance, we might have just a couple of reads mapped and both have the alternate allele at a locus; this could be a homozygote for the alternate allele, or it could be that we just haven't sampled one allele of a heterozygote. For this reason, many newer programs use genotype likelihoods, instead of genotype calls, to carry out population genomic analyses. Good programs that use this information include GPAT++ and Angsd. An alternative to using genotype likelihoods is to have strong filters for genotype likelihood to call a SNP (i.e. only call SNPs that we are very confident about). That's the strategy we are using here.

##Step 6: Format SNP Matrix in R

6a) Open R program

6b) In a new document, paste the following

```
setwd('~/Desktop/SNPcalling_tutorial')
snps<-read.delim('ahy_snps.012',header=F,na=-1,row.names=1)
pos<-read.delim('ahy_snps.012.pos',header=F)
indv<-read.delim('ahy_snps.012.indv',header=F)

colnames(snps)<-paste(pos[,1],pos[,2],sep=':')
rownames(snps)<-indv[,1]
snps<-as.matrix(snps)

#read in meta data
meta<-read.delim('meta.txt')
```


#How to use genotype data:

##PCA
Plot a PCA to visually identify any clusters within your data. Are these clusters associated with your meta data?

In R

```
pc.out<-prcomp(snps)
summary(pc.out)
plot(pc.out$x[,1],pc.out$x[,2],col=meta$Pool,xlab='PC1',ylab='PC2',pch=19)
legend('topright',legend=unique(meta$Pool),fill=c('black','red','green'))
```

##Fst

In R

```
install.packages('hierfstat')
library(hierfstat)

#prepare 0,1,2 matrix in hierfstat format
#we use our pca to separate samples into clusters to test for genetic differentiation
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
hist(site.fst,xlab='Fst',ylab='Counts',main='Distribution of Fst between PC1 clusters',col='grey')
```

## Other interesting analyses

- structure/admixture: ngsAdmix takes bam files
- Outliers & environmental data: outFLANK, bayenv, etc…
- local & global linkage: vcftools or R linkage package
- somatic mutations: are your samples high depth from the same individual? if so, you could look at this
- dN/dS: ORF prediction with biopython scripts or snpEff
- eQTLs: are SNPs associated with expression?

## How to use your own data with these scripts

- This pipeline runs quickly because we subset the data into 10 contigs. With your data, it may be difficult to run through this on a single computer, but it is possible if you have enough memory. We recommend accessing a computer cluster (most universities have a cluster). 
	- There are a few accessible computer clusters online: Amazon cloud, NSF-XSEDE 
- The Palumbi lab's internal use parallelized bash scripts for transcriptome assembly, mapping, gene expression, and SNP calling are available formatted for a SLURM scheduler at https://github.com/bethsheets/Population-Genomics-via-RNAseq, but they may not be applicable to your specific computer system or needs.
	