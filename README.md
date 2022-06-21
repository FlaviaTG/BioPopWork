# BioPopWork
 
 
## Workshop: fast bioinformatics workflow for population genomics
 
This project is the workflow for a workshop in Bioinformatics for population genomics analyses with resequencing data and the reference genome.
The workshop is for practicing how to obtain formats and analyses for fast-and-first exploration of resequencing data of populations. We will reflect on the methods and alternatives used in population genomics analyses during a paper discussion on each workshop session.

## Aims

- Familiarize the student with programming languages used repeatedly in population genomics analyses.
- Provide alternatives to solve recurring formatting problems or develop databases for analyses. 
- Promote creative thinking on how to prepare analyses to answer specific questions.
- Promote conscientious use of bioinformatics resources such as programs and models to maintain critical thinking regarding biological interpretations and recognize the limitations of each tool.

## Programs:

We will create two conda environments and we will also use the module system.

### Ask for an interactive session in TACC
```
idev -p normal -m 180 -A Bioinformatics-Works -t 02:00:00 -N 1 -n 68
```
### Create folders on your $SCRATCH path
```
cd $SCRATCH
mkdir Data
mkdir Alignment
mkdir Ancestry
mkdir Demography
mkdir Genotypes
mkdir ref-genome
mkdir Selection
```
### Create conda env.

Download miniconda from here
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```
Follow the instructions, type yes when needed, and enter
```
bash Miniconda3-latest-Linux-x86_64.sh
```
Logout and log back in using the ssh command.
```
conda config --set auto_activate_base false
```
### Install ANGSD, bcftools, samtools, picard and others
```
conda create --name ANGSD
conda activate ANGSD
conda install -c bioconda angsd
conda install -c bioconda samtools openssl=1.0
conda install -c bioconda/label/cf201901 bcftools 
conda install -c bioconda picard
```
For vcftools we need another environment.
```
conda create --name SUMMARY
conda activate SUMMARY
conda install -c bioconda/label/cf201901 vcftools
conda install -c bioconda pgdspider
conda install -c genomedk psmc
conda install -c bioconda admixture
conda install -c bioconda plink2 
```
### Install BayPass from : `http://www1.montpellier.inra.fr/CBGP/software/baypass/`
```
# Use gcc
module load gcc/9.1.0

# Download and extract to SCRATCH
curl http://www1.montpellier.inra.fr/CBGP/software/baypass/files/baypass_2.3.tar.gz | tar --directory=$SCRATCH -xzf -

# Find the sources
cd $SCRATCH/baypass_2.3/sources

# Build the source
# Note the optimizations used by the developers are ill-advised
# This code has a lot of red flags
make FC=gfortran

# You could try these commands instead
# The use of '-fno-range-check' and other flags are probably hiding bugs
gfortran -O2 -fno-range-check -c mt_kind_defs.F90 -o mt_kind_defs.o
gfortran -O2 -fno-range-check -c mt_stream.F90 -o mt_stream.o
gfortran -O2 -fno-range-check -c gf2xe.F90 -o gf2xe.o
gfortran -O2 -fno-range-check -c f_jump_coeff.F90 -o f_jump_coeff.o
gfortran -O2 -fno-range-check -c f_get_coeff.F90 -o f_get_coeff.o
gfortran -O2 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -fopenmp    baypass.f90 mt_stream.o f_get_coeff.o gf2xe.o -o g_baypass

# Stash the binary in a convenient location
mkdir -p $WORK/apps/baypass_2.3
cp -p g_baypass $WORK/apps/baypass_2.3
cd $WORK/apps/baypass_2.3
./g_baypass

# Clean up (optional)
rm -r $SCRATCH/baypass_2.3
```
This installation will generate an executable script. You will need the path of this executable script to call the program. See below the procedure.
### Install BayeScan from : `http://cmpg.unibe.ch/software/BayeScan/download.html`
```
# Use gcc
module load gcc/9.1.0

# Lets build on scratch
cd $SCRATCH
wget http://cmpg.unibe.ch/software/BayeScan/files/BayeScan2.1.zip
unzip BayeScan2.1.zip
cd BayeScan2.1/source

# Build it
make

# If the above does not work, type
g++ -fopenmp -lpthread -o bayescan_2.1 start.o beta.o dirichlet.o RJupdates.o MHupdates.o likelihood.o read_write.o anyoption.o 

# Install
mkdir -p $WORK/apps/BayeScan2.1
cp -p bayescan_2.1 $WORK/apps/BayeScan2.1

# Clean up (optional)
cd $SCRATCH
rm -r BayeScan2.1.zip BayeScan2.1
```
You will use the binaries files BayeScan2.1_linux64bits in BayeScan2.1/binaries folder
### Prepare your R environment
Load the modules and run R

```
module load intel/18.0.0  impi/18.0.0
module load Rstats/3.5.1
R
```
### Install all R packages needed
```
# Compile a little faster (8 cores)
Sys.setenv(MAKE = 'make -j8')

# Don't ask for repo
options(repo = "https://cloud.r-project.org/")

# To make this permanent use
# cat('options(repo = "https://cloud.r-project.org/")\n', file = "~/.Rprofile", append = TRUE)

install.packages("BiocManager")
BiocManager::install("SNPRelate")
install.packages("tidyverse")
install.packages("ggplot2")
BiocManager::install("gdsfmt")
install.packages("devtools")
install.packages("data.table")
install.packages("qqman")
install.packages("filesstrings")
install.packages("cowplot")
install.packages("GenWin")
install.packages("reshape2")
install.packages("gghighlight")
install.packages("corrplot")
install.packages("ape")
install.packages("geigen")
install.packages("mvtnorm")
install.packages("magrittr")
install.packages("dplyr")
install.packages("boa")
```
## Get the files and storage in the appropriate directory just created
Download reference genome from NCBI by searching for Corvus hawaiiensis genome
```
cd $SCRATCH/ref-genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/740/725/GCF_020740725.1_bCorHaw1.pri.cur/GCF_020740725.1_bCorHaw1.pri.cur_genomic.fna.gz
```
Make softlinks of fastq files and bam files. This way, you copy the file as a link to avoid duplicating the entire file.
```
cd $SCRATCH/Data
ln -s /work2/08209/brian97/shared_reseq/*gz .
cd $SCRATCH/Alignment
ln -s /work2/08752/ftermig/shared_workshop/bamfiles/*bam .
```
## INPUT: fastq 

Pair-end sequencing data, R1.fq and R2.fq files per samples

The first step is to check for read quality using Fred scores by running fastqc on each file. After evaluating if the reads do not drop into the "red zone" in the final plots, you could decide if the reads need to be trimmed.

### Run the modules or environment
```
module load fastqc/0.11.5
fastqc 01.180603.B03.S1.1.fastq.gz -o ./QUALITY-reads
```
In a loop and a slurm job, see bash script for an example `Job-FastqC.sh`. This job could take 2 days
```
for i in *.fastq.gz; do fastqc $i -o ./QUALITY-reads;done
```
## 1) Alignment
Example of BWA alignment on one individual. 
- Sequencing plataform ID: @D00742:CCJ7KANXX
- Individual barcode: ACGGTC. This individual barcode will be the ID of that sample on the final .sam .bam file
Example of a complete fastq header:

```
@A00672:51:HWJMYDSXY:1:1101:1669:1000 2:N:0:NCACGGAC+TGCGAGAC
```

You need to rename files to remove the unnecessary extensions and have the pair extension R1/R2 reachable to make loops or call multiple files simultaneously.
e.

```
rename L001.R1.001.fastq.gz '1.fastq.gz' *.fastq.gz
rename L001.R2.001.fastq.gz '2.fastq.gz' *.fastq.gz
```
### For alignments with BWA. 

This alignment could take 9 hours per sample. You could do many alignments jobs per sample to run simultaneously or a loop to make one by one. See example job `Job-FastqC.sh`

### Explore the header of your fastq file
Look into the platform ID and the individual barcode used for the sequencing.

```
for i in *.fastq.gz; do zcat $i | head -n 4 | grep '@'; done
```
### Make your genome path reachable as a variable

```
ref=/scratch/08752/ftermig/ref-genome
```
### Run the modules or environment

```
module load intel/17.0.4
module load bwa/0.7.16a
```

### Make your reference genome available in path and index with bwa. 
Creating the index will take 4 hours. Run it in a job. See example `job-refgenome-index.sh`

```
GB=/scratch/08752/ftermig/ref-genome/GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna
bwa index -p CorHaw $GB
```
### Now, you can run the alignments
The alignment of DNA resequencing to a reference genome could take 9 hours.

```
bwa mem -t 30 -M -R '@RG\tID:@A00672:51:HWJMYDSXY\tSM:CTAGCGCT' $ref/CorHaw 01.180603.B03.S1.1.fastq.gz 01.180603.B03.S1.2.fastq.gz > sample.CTAGCGCT.aln.sam

```
After alignment, need to sort, index, and compress sam to bam.
### Run the modules or environment
```
conda activate ANGSD
```
### Run Picard to create index, compress, and sort by coordinates

```
picard -Xmx128g -XX:ParallelGCThreads=32 SortSam I=sample.ACAGGCGC.aln.sam O=sample.ACAGGCGC.aln.sam.sort.bam SORT_ORDER=coordinate CREATE_INDEX=true
```
All in a loop and run it in a job. See example `Job-sam-index-sort-bam-1.sh`. This could take 4 hours per sample

```
for sample in *.aln.sam;do picard -Xmx128g -XX:ParallelGCThreads=32 SortSam -I $sample -O $sample.sort.bam SORT_ORDER=coordinate CREATE_INDEX=true;done
```
### Check if your bam file is sorted with samtools
#### Load samtools first

```
module load intel/17.0.4
module load samtools/1.5
samtools view -h sample.TGCGAGAC.aln.sam.sort.bam.dedup.bam | head
```

After sort and index need to mark/remove optical duplicates. This could take 6-9 hours. Run it in a Job. See example `Job-mark duplicates-1.sh`
### Mark/remove optical duplicates
```
picard -Xmx2g -Xms1g -XX:ParallelGCThreads=3 MarkDuplicates TMP_DIR=tmp I=sample.GAACCGCG.aln.sam.sort.bam O=sample.GAACCGCG.aln.sam.sort.bam.dedup.bam METRICS_FILE=sample.GAACCGCG.aln.sam.sort.bam.dedup.bam.metrics.txt MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true TAGGING_POLICY=All
```
### Basic stats
Obtain basic stats of the alignments with samtools, total reads mapped and unmapped, and the number of reads with good quality. The stats could take longer. Run all in a loop and create a new file with all the info. This could take about a day. See example job in: `Job-countReads.sh` and `Job-Depth.sh`.
For plots see `Plot-mappedReads.R`

```
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam
samtools view -c -f 4 sample.TGCGAGAC.aln.sam.sort.bam
samtools stats sample.TGCGAGAC.aln.sam.sort.bam
```
### Explore the alignments
Make index for reference
```
samtools faidx GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna
```
### See the alignments
```
samtools tview sample.TGCGAGAC.aln.sam.sort.bam $GB
```
You can see a specific region, chromosomes CM036346.1. you can look into the index file of the reference to get the name of the chromosomes you would like to look at
```
samtools tview -d T -p CM036346.1:300 sample.TGCGAGAC.aln.sam.sort.bam $GB
```
### Screen results:
- dot: means a base that matched the reference on the forward strand
- comma: means a base that matched the reference on the reverse strand
- asterisk: is a placeholder for a deleted base in a multiple bases
- upper case: denotes a base that did not match the reference on the forward strand
- lower case: denotes a base that did not match the reference on the reverse strand

## 2) Prep-alignments for downstream analyses
#### Make the list of autosomes and sex chromosomes, and unplaced scaffolds.
Here you need to select one or two chromosomes for all the downstream analyses.
```
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep 'JAJGSY' >> unplaced-scaffold.txt
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep 'CM' >> placed-scaffold.txt
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep -e 'CM036387.1' -e 'CM036388.1' >> ZW-chr.txt
cat placed-scaffold2.txt | grep -v -e 'CM036387.1' -e 'CM036388.1' >> AUTOSOMES-chr.txt
```
#### Use the list to split the bam files. 

Try one and then make folders and a loop to split the bam files into different genomic regions or chromosomes. 

```
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep -w -f AUTOSOMES-chr.txt | xargs samtools view -b sample.TGCGAGAC.aln.sam.sort.bam > ./AUTOSOMES-CHR/sample.TGCGAGAC.aln.sam.sort.bam.AUTOSOMES.bam
```
In a loop could take 4 hours.Run it in a job. See example `Job-samtools-split.sh`
```
for file in *.aln.sam.sort.bam ; do samtools idxstats $file | cut -f 1 | grep -w -f AUTOSOMES-chr.txt | xargs samtools view -b $file > ./AUTOSOMES-CHR/$file.AUTOSOMES.bam; done
```
#### Need to sort and index 
Every time you create a new bam, it must be indexed and sorted for downstream analysis. See job example `Job-sort-index.sh`
```
for i in *bam; do samtools sort $i -o $i.sorted.bam;done
for i in *sorted.bam; do samtools index $i;done
```
## 3) Generate the genotypes by likelihood
 
We are going to perform genotype likelihoods in a fast way using ANGSD. This method is the quickest and best approach when samples have variable x coverage. Some information will be negleckt in the final vcf file with this method, such as indels. The mapping quality of reads needs to be filtering during the ANGSD run.

#### Call de environment ANGSD
```
conda activate ANGSD
```
#### Create a list file of your bam 
The list needs to have the files for each sample with the path. See example in `bam-list-unplaced.txt`
This is an old version GATK like genotype likelihoods by using the following flags
- -doGlf 2: binary glf
- -doMajorMinor 1: Infer major and minor from GL
- -SNP_pval 1e-6: If we are interested in looking at allele frequencies only for sites that are variable in our sample.
- -minMaf         0.05        (Remove sites with MAF below) You want to remove major alleles with very low frequency likely to be sequencing errors. You could change this value if you expect to have a variant with shallow frequency.
- -SNP_pval       0.000001       (Remove sites with a pvalue larger) We can consider assigning SNPs sites whose estimated allele frequency is above a certain threshold (e.g., the frequency of a singleton) or whose probability of being variable is above a specified value.
- -doMaf 1: Frequency (fixed major and minor)
- -minInd: you want at least SNP present in 75% of the individuals  75% of 20 individuals is = 15 individuals
- -minMapQ 30: here is the mapping quality filtering. 
- -minQ 20: here is the minimum base quality score filter. 
```
angsd -GL 2 -doBcf 1 -out genolike2-greenjay_ZW -nThreads 68 -doPost 1 -docounts 1 -dogeno 1 -minInd 15 -doMajorMinor 1 -SNP_pval 1e-6 -minMaf 0.05 -doMaf 1 -minMapQ 30 -minQ 20 -bam bam-list-ZW.txt
```
If the above line works well run it in a job. See job example `Job-ANGSD-genotypes.sh`
The output is a bcf format, the compressed version of vcf files.
``
### Explore your vcf file
First convert the compressed form .bcf into vcf decompressed form
```
bcftools convert -O v genolike2-greenjay_ZW.bcf -o genolike2-greenjay_ZW-test.vcf
```
Look how the header looks like in the vcf file
```
bcftools view --header-only genolike1-greenjay-UNPLACED.vcf
```
### Reheader with bcftools
Reheader because ANGSD place individual names into the file format according to the list order in the bam-file-list.txt. It needs to change header name because header in this vcf file is not detected by vcftools. Make a header file `header-order-genotype-names.txt` according to the order in `bam-file-list.txt` and use the name format for the samples you prefer most. But keep it simple without special characters.

```
bcftools reheader -s header-order-genotype-names.txt  genolike1-greenjay_AUTOSOMES.vcf > newH-genolike1-greenjay_AUTOSOMES.vcf
```
### Check if the header added looks good
```
bcftools view --header-only newH-genolike1-greenjay_AUTOSOMES.vcf
```
## 4) Summary statistic calculations
### Calculate population genomics summary statistic genome-wide
Such as Fst, Pi, and relatedness statistics.
### Make variant database without missing data
Select just variants with information present in all individuals. Keep a database 100% complete.
First call the environment SUMMARY.

```
conda activate SUMMARY
vcftools --vcf newH-genolike1-greenjay_AUTOSOMES.vcf --max-missing 1.0 --recode --recode-INFO-all --out newH-genolike1-greenjay_AUTOSOMES.NOmissing
```
Run relatedness on the method of Manichaikul et al., BIOINFORMATICS 2010.
```
vcftools --vcf newH-genolike1-greenjay_AUTOSOMES.vcf --relatedness2
```
### Create the population map for Fst calculation. 
Create a `pop-map-north.txt` and `pop-map-south.txt` file per population listing the sample names. The resulting output file has the suffix ".fst". Use a small window size to calculate Fst. Use the spline window technique in R to visualize the Fst or Pi distribution in a manhattan-like plot. Run the summary statistics Fst calculations. This Fst estimate is from Weir and Cockerham’s 1984 paper. The preferred calculation of Fst. The provided file must contain a list of individuals per line from the VCF file corresponding to one population. The flag "--keep" can be used to provide a list of samples per population. By default, calculations are on a per-site basis. The output file has the suffix ".weir.fst".
```
vcftools --vcf newH-genolike1-greenjay_AUTOSOMES.NOmissing.recode.vcf --weir-fst-pop pop-map-north.txt --weir-fst-pop pop-map-south.txt --fst-window-size 100 --out newH-genolike1-greenjay_AUTOSOMES-NOmissing-100bp
```
With the output you can perform the spline window technique for a Fst genome-wide visualization. See step 6).
### Calculate nucleotide diversity per bin of size 10bp
```
vcftools --vcf newH-genolike1-greenjay_AUTOSOMES.NOmissing.recode.vcf --window-pi-step 10bp --out Pi10bp 
```
### Calculate ts/tv 
With the ts/tv ratio, you can compare substitutions overall between populations, or you could see if the genotype went well if the organism ratio is already known. Calculated with a database 100% complete without missing data
```
vcftools --vcf newH-genolike1-greenjay_AUTOSOMES.NOmissing.recode.vcf --TsTv-summary --out TSTV-AUTOSOMES.NOmissing
```
### Prune SNPs every 10,000kb
For some analyses, you will need to select putative unlinked variants. PCA, scan for selection, and ancestry proportion analyses need a non-linked SNPs database.
```
vcftools --vcf newH-genolike1-greenjay_AUTOSOMES.NOmissing.recode.vcf --thin 10000 --recode --recode-INFO-all --out newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN
```
It is possible to create many formats for downstream analysis from the vcf file and with the appropriate variant filter.
### Create plink format for admixture program and others
```
vcftools --vcf newH-genolike2-greenjay_ZW-test.NOmissing.THIN.recode.vcf --plink 
```
### Creat formats for BayPass
This program needs a specific database format. We are going to use some simple steps to get into that format. First, you need to create the .GESTE file. It is the same as for the BayeScan. It does population counts per allele, and then we will join both alleles in the same file with all populations.
We will use the PGDspider program to format the vcf file to GESTE format. For that, we need to create a .spid file for the vcf format to input into the program. Create the .spid by giving all the arguments without the -spid file. The resulted .spid file template needs to be edited by answering questions related to the formats. See file `template_VCF_GESTE_BAYE_SCAN.spid`.
```
PGDSpider2-cli -inputfile newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN.recode.vcf -inputformat VCF -outputfile newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN.GESTE -outputformat GESTE_BAYE_SCAN
```
Run PGDspider again with the new edited .spid file `template_VCF_GESTE_BAYE_SCAN.spid`
```
PGDSpider2-cli -inputfile newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN.recode.vcf -inputformat VCF -outputfile newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN.GESTE -outputformat GESTE_BAYE_SCAN -spid template_VCF_GESTE_BAYE_SCAN.spid
```
Now we do some bash commands step by step to check one by one we are formating the data correctly
##### Remove the first lines with info not necessary for the BayPass format
```
sed -e '1,4d' newH-genolike2-greenjay_ZW-test.NOmissing.THIN.GESTE > test2.GESTE
```
##### Split files in populations
We have two populations that are separated by new lines.
```
sed '/^$/q' test2.GESTE > test2-pop1.txt
sed '1,/^$/d' test2.GESTE > test2-pop2.txt
```
##### Remove the empty lines at the end of each file
```
sed -i '/^$/d' test2-pop2.txt
sed -i '/^$/d' test2-pop1.txt
```
##### Concatenate by column files of pop1 and pop2 in columns
First, select just the column 4 with two alleles info
```
cut -f4 test2-pop1.txt > format-test2-pop1.txt
cut -f4 test2-pop2.txt > format-test2-pop2.txt
```
##### Paste the two populations files into columns and separate them by space
```
paste -d '' format-test2-pop1.txt format-test2-pop2.txt > BayPass-format-test-ZW-thin.txt
```
The first line with the headers needs to be removed for the final format.
```
sed -i '1,1d' BayPass-format-test-ZW-thin.txt
```
Now you have ready the genotype format per population to be used as input for BayPass. Consult the manual of the programs for more details.

### Index vcf file
For some other programs the vcf file needs to be indexed. First, we need to compress the vcf file and later index with bcftools.
```
bgzip newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN
bcftools index newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN
```
### Count SNPs per chromosome
Visualization of this data can be done with R code example `Plot-SNP-per.CHR.R`
```
zcat newH-genolike1-greenjay_AUTOSOMES.NOmissing.THIN | grep -v "^#"  | cut -f 1 | sort | uniq -c > SNP-per-CHR-Nomissing.txt
```

## 5) Exploring data with PCA analysis. 
We will use the SNPrelate R packages for PCA and other analyses. First, the vcf file needs to be formated to gds in R. See Job `Job-gdsFormat.sh` and R code/script `PCA-SNPrelate-example.R` to make the gds format and run PCA, kinship, and other analyses
## 6) Spline-window technique for genome-wide stats visualization
It is commonly used in population genomics with the window size of 10kb-100kb to visualize the distribution of summary statistics such as Fst. This is helpful to make a manhattan plot and see the genomic regions where some selective pressures could be operating in the populations. Choosing a window size i arbitrary and depends on the density of the data. The spline-window technique uses a statistical value to find boundaries on the data. It gives - depending on the data -a variable window size according to the values of the statistic in use.
#### Calculate window size for each chromosome. 
Run the spline technique per chromosome/region: after calculating fst in vcf tools use your *.weir.fst and split it in chromosomes. You can do that on R or bash. For an R version see first lines on file `GeneWin-example.R`. In bash is also fast, you already have a list file with chromosomes names. Used it similar to as above with the bam file but this time is easier you don't need samtools because the file is just a *txt file. Create a new directory GeneWin to save all data for GeneWin analysis
```
for chr in $(cat unplaced-scaffold.txt); do grep -w $chr 1Kkb.windowed.weir.fst > ./GeneWin/$chr.Kkb.windowed.weir.fst; done
```
#### Plot in a manhattan-like plot Fst genome-wide (all chromosomes). see R code `GeneWin-example.R` and `GeneWin-plots.R`
For visualization of the generated windows size you need to create a file for R that indicates the file name of each chromosome and the chromosome number that corresponds to that file. See file `files_MATCH-CHR_example.txt` for example. In this file, you need to change the first column for the file's complete name for that chromosome. You will need this file to plot within R with the conde in file `GeneWin-plots.R`. It is a straightforward file; create that file in the way you feel more comfortable. It could be excel or in a terminal.

## 7) Coalescent inferences of population demography
PSMC takes the consensus fastq file, and infers the history of population sizes. The first step starts from mapped reads and is to produce a consensus sequence in FASTQ format. We will use the samtools/bcftools, following the methods described in the paper of Palkopoulou et al., 2015, with default parameters for model fitting.

#### Load version samtools 1.5. 
Don't use a higher version of samtools. Load the reference genome.

```
module load intel/17.0.4
module load samtools/1.5
GB=/scratch/08752/ftermig/ref-genome/GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna
```
#### Make your list of chromosomes an array
Check if the array is working by making a small loop into the array
```
mapfile -t CHR < unplaced-scaffold.txt

for str in ${CHR[@]}; do
  echo $str
done
```
#### Produce a consensus sequence per chromosome in one sample
If you want to run all samples, make a loop in a job. See job example `Job-all-mpileup.sh`
```
for str in ${CHR[@]}; do
samtools mpileup -Q 30 -q 20 -u -v \
-f $GB -r $str sample.TGCGAGAC.aln.sam.sort.bam.UNPLACED.bam.sorted.bam |  
bcftools call -c |  
vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > sample.TGCGAGAC.aln.sam.sort.bam.UNPLACED.bam.$str.fq
done
```
#### Concatenate all concensus sequences 
All chromosome/scaffold in one fastq per sample
```
cat sample.TGCGAGAC*JAJGSY0*.fq > sample.TGCGAGAC.aln.sam.sort.bam.UNPLACED.consensus.fq
```
If one sample works well, try to make another loop to create the consensus sequences per sample with all chromosomes.
```
mapfile -t samples < samples-barcode.txt

for str in ${samples[@]}; do
  echo $str
done

for str in ${samples[@]}; do
cat sample.$str.*CM0*.fq > sample.$str.aln.sam.sort.bam.ZW.consensus.fq
done
```
#### Create format for PSMC
If you want to run all samples, make a loop in a job. See job example ``
```
fq2psmcfa sample.TGCGAGAC.aln.sam.sort.bam.UNPLACED.consensus.fq > sample.TGCGAGAC.UNPLACED.consensus.psmcfa
```
How would you do it in a loop for all samples?

#### Run PSMC
Run the models in one sample. If you want to run all samples make a loop in a job. See job example ``
```
psmc -p "4+25*2+4+6" -o sample.TGCGAGAC.UNPLACED.consensus.psmc sample.TGCGAGAC.UNPLACED.consensus.psmcfa
```
for
#### Make the plot
With the generated data, its possible to make your costume plot using -R flag. We will try the plot that comes with PSMC program and generate the data for a custom plot in R using the generated file with extension ".0.txt." See R code/script `PSMC-costum-plot.R` for custom plots of all individuals in one plot.
```
psmc_plot.pl -R -u 0.221e-8 -g 1 Green-jay_TGCGAGAC_UNPLACED_plot sample.TGCGAGAC.UNPLACED.consensus.psmc
```

## 8) Genome-Environmental-Association analysis
Look for regions in the genome that are associated to environmental conditions. Above you have generated the BayPass format file and you also have the covariance environmental matrix file per population ``
### load requirements for ByPass
```
module load gcc/9.1.0
BayPass=/work2/08752/ftermig/stampede2/apps/baypass_2.3/sources
```
First, you need to scale your variables. This could take 40min hour for 2 chromosomes
```
$BayPass/g_baypass -npop 2 -gfile ./BayPass-format-test-ZW-thin.txt -efile ./cov-Bio1-pop.txt -scalecov -outprefix Scale-var-Bio1  
```
Now you need to create an omega file obtained by a first analysis under the core model or the IS covariate mode with your scaled variables
```
$BayPass/g_baypass -npop 2 -gfile ./BayPass-format-test-ZW-thin.txt -efile ./Scale-var-Bio1_covariate.std -omegafile Scale-var-Bio1_mat_omega.out -outprefix anacoreZW
```
Try to run the core model to compare it with the AUX model and you need to use the omega file generated above
```
$BayPass/g_baypass -nthreads 8 -npop 2 -gfile BayPass-format-test-ZW-thin.txt -efile Scale-var-Bio1_covariate.std -covmcmc -omegafile Scale-var-Bio1_mat_omega.out -outprefix auxcoveZW
```
The above commands will generate files for the final GEA analysis under model AUX. These are going to be the input for the last Run.
```
$BayPass/g_baypass -nthreads 8 -npop 2 -gfile BayPass-format-test-ZW-thin.txt -efile Scale-var-Bio1_covariate.std -auxmodel -omegafile Scale-var-Bio1_mat_omega.out -outprefix Aux-var-Bio1
```
#### Load R scripts to evaluate models and plot final results
To assess the models and plot final results, you need to use the source R code `baypass_utils.R`. Examples on how to use it are in the file `BayPass-PLOTS.R`
You can get the variant outliers from the betai.out output by selecting the columns like this
```
less 2anaux1B_summary_betai.out | tr -s '\ '| awk -F ' ' '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}'| sort -m | awk '{ if ($6 > 3) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' > Outliers-Bio1-aux.txt
```
you can be more specific and select outliers with a certain BF value > 10 for STRONG STRENGTH SELECTION
```
less Outliers-Bio1-aux.txt |awk '{ if ($6 > 10) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}}' > Strong-Outliers-Bio1-aux.txt
```
Now you need to track back the exact variant under selection and position by going into the pre-formating input into BayPass.

## 9) Genome scan for selection (multinomial-Dirichlet model)
### run Bayescan
Now that you have generated the input for BayeScan, you can run the program to give it a try.
BayeScan could take a while to run. Run it in a job see example `Job-Bayescan.sh`
This is an old bayesian method based on the multinomial-Dirichlet model. It looks for divergent selection under an island model in which subpopulation allele frequencies (measured by Fst coefficient) are correlated through a common migrant gene pool. This program formulation can consider realistic ecological scenarios where the effective size and the immigration rate may differ among subpopulations. This program can be used as a base line of analyses to explore the data because the format is simple and fast to get. See details in: http://cmpg.unibe.ch/software/BayeScan/index.html.

```
module load gcc/9.1.0
SCAN=/scratch/08752/ftermig/programs/BayeScan2.1/binaries
$SCAN/BayeScan2.1_linux64bits newH-genolike2-greenjay_ZW-test.NOmissing.THIN.GESTE -n 5000 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 -threads 40 -out_freq -od ./BAYESCAN/
```
The fst output can be plotted by using the bayesian R plot function. See examples on how to use it in R code `BayeScan-plot.R`
## 10) Ancestry proportions (population structure)
### run ADMIXTURE
But first you need to creat a .bed file format for that program run with plink2. Plink and admixture are in your conda environment SUMMARY
```
conda activate SUMMARY
plink2 --vcf newH-genolike2-greenjay_ZW-test.NOmissing.THIN.recode.vcf --geno 0.9 --recode --no-fid --no-parents --no-sex --no-pheno --out newH-genolike2-greenjay_ZW-test.NOmissing.THIN --make-bed --allow-extra-chr 0
```
Now you can run admixture test. If it works well run more K
```
admixture newH-genolike2-greenjay_ZW-test.NOmissing.THIN.bed 5
```
Now in a loop to run more K and validations file values. You can run this in a slurm job, see job `Job-ADMIXTURE.sh`
```
for K in 1 2 3 4 5 6 7 8 9 10; do admixture -B2000 -j40 --cv newH-genolike2-greenjay_ZW-test.NOmissing.THIN.bed $K | tee Bootlog${K}.out; done
```
With the validation file and the proportions of ancestry Q files, you can plot the best K value by cross-validation and a barplot, see R code in `ADMIXTURE-plot.R`
You need to format your file to make it easy for R to plot, by putting together all K validation files like this:
```
grep -h CV Bootlog*.out > cross.val.txt
```
Check the `cross.val.txt` for spaces or extra characters that could interfere with R plotting. After that, you can use it to plot in R using code `ADMIXTURE-plot.R`
