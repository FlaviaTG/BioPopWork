# BioPopWork
 
 
## Workshop: fast bioinformatics workflow for population genomics
 
This project is the workflow for a workshop in Bioinformatics for population genomics analyses with resequencing data and the reference genome.
The workshop is for practicing how to obtain formats and analyses for fast-and-first exploration of resequencing data of populations. We will reflect on the methods and alternatives used in population genomics analyses during a paper discussion on each workshop session.

## Aims

- Familiarize the student with programming languages used repeatedly in population genomics analyses.
- Provide alternatives to solve recurring problems in formatting or developing databases for analyses. 
- Promote creative thinking on how to prepare analyses to answer specific questions.
- Promote conscientious use of bioinformatics resources such as programs and models. To maintain critical thinking when it comes to biological interpretations and recognize the limitations of each tool.

## Programs:

We will create an environment and will use the module system.

### Create conda env.

Download miniconda from here
```
https://docs.conda.io/en/latest/miniconda.html#linux-installers
```
Follow the instructions, typing yes when needed, and enter
```
bash Miniconda3-py39_4.9.2-Linux-x86_64.sh
```

Logout and log back in using the ssh command.
```
conda config --set auto_activate_base true
conda create --name ANGSD
```
### Install ANGSD, bcftools, vcftools, samtools, picard and others
```
conda activate ANGSD
conda install -c bioconda angsd
conda install -c bioconda/label/cf201901 bcftools 
conda install -c bioconda samtools openssl=1.0
conda install -c bioconda picard
conda install -c bioconda pgdspider 

```
### install glactools from : `https://github.com/grenaud/glactools` in your $WORK directory

```
module load gcc/9.1.0
git clone --depth 1 https://github.com/grenaud/glactools.git
cd glactools
make
```
### Prepare your R environment
Load the modules and run R

```
module load intel/18.0.0  impi/18.0.0
module load Rstats/3.5.1
R
```
### Install all R packages needed
```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPRelate")
install.packages("tidyverse")
install.packages("ggplot")
```
## INPUT: fastq 

Pair-end sequencing data, R1.fq and R2.fq files per samples

The first step is to check for read quality using Fred scores by running fastqc on each file. After evaluating if the reads do not drop into the "red zone" in the final plots, you could decide if the reads need to be trimmed.

### Run the modules or environment
```
module load fastqc/0.11.5
fastqc 01.180603.B03.S1.1.fastq.gz -o ./QUALITY-reads
```
In a loop and in a slurm job see bash script for an example `Job-FastqC.sh`. This job could take 2 days
```
for i in *.fastq.gz; do fastqc $i -o ./QUALITY-reads;done
```
## Alignment
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

### Explore header of your fastq file
Look into the platform ID and the individuals barcode used for the sequencing.

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
### Now you can run the alignments
The alignment of DNA resequencing to a reference genome could take 9 hours

```
bwa mem -t 30 -M -R '@RG\tID:@A00672:51:HWJMYDSXY\tSM:CTAGCGCT' $ref/CorHaw 01.180603.B03.S1.1.fastq.gz 01.180603.B03.S1.2.fastq.gz > sample.CTAGCGCT.aln.sam

```
After alignment, need to sort, index, and combpress sam to bam.
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
### Check if your bam file is really sorted
```
samtools view -h sample.TGCGAGAC.aln.sam.sort.bam.dedup.bam | head
```

After sort and index need to mark/remove optical duplicates. This could take 6-9 hours. Run it in a Job. See example `Job-mark duplicates-1.sh`
### Mark/remove optical duplicates
```
picard -Xmx2g -Xms1g -XX:ParallelGCThreads=3 MarkDuplicates TMP_DIR=tmp I=sample.GAACCGCG.aln.sam.sort.bam O=sample.GAACCGCG.aln.sam.sort.bam.dedup.bam METRICS_FILE=sample.GAACCGCG.aln.sam.sort.bam.dedup.bam.metrics.txt MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true TAGGING_POLICY=All
```
### Basic stats

Obtain basic stats of the alignments with samtools, total reads mapped and unmapped, and the amount of reads with good quality. The stats could take longer. Run all in a loop and create a new file with all the info. This could take about a day. See example job in: `Job-countReads.sh` and `Job-Depth.sh`

```
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam
samtools view -c -f 4 sample.TGCGAGAC.aln.sam.sort.bam
samtools stats sample.TGCGAGAC.aln.sam.sort.bam > depth-TGCGAGAC.txt
```
### Explore the alignments
Make index for reference
```
samtools faidx GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna
```
### See the alignments in general as a whole
```
samtools tview sample.TGCGAGAC.aln.sam.sort.bam $GB
```
You can see a specific region, chromosomes CM036346.1. you can look into the index file of the reference to get the name of the chromosomes you will like to look at
```
samtools tview -d T -p CM036346.1:300 sample.TGCGAGAC.aln.sam.sort.bam $GB
```
### Screen results:
- dot: means a base that matched the reference on the forward strand
- comma: means a base that matched the reference on the reverse strand
- asterisk: is a placeholder for a deleted base in a multiple bases
- upper case: denotes a base that did not match the reference on the forward strand
- lower case: denotes a base that did not match the reference on the reverse strand

## Prep-alignments for downstream analyses
### Make the list of autosomes and sex chromosomes and unplaced scaffolds.
```
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep 'JAJGSY' >> unplaced-scaffold.txt
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep 'CM' >> placed-scaffold.txt
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep -e 'CM036387.1' -e 'CM036388.1' >> ZW-chr.txt
cat placed-scaffold2.txt | grep -v -e 'CM036387.1' -e 'CM036388.1' >> AUTOSOMES-chr.txt
```
### Use the list to split the bam files. 

Try one and then make folders and a loop to split the bam files into different genomic regions or chromosomes. 

```
samtools idxstats sample.TGCGAGAC.aln.sam.sort.bam | cut -f 1 | grep -w -f AUTOSOMES-chr.txt | xargs samtools view -b sample.TGCGAGAC.aln.sam.sort.bam > ./AUTOSOMES-CHR/sample.TGCGAGAC.aln.sam.sort.bam.AUTOSOMES.bam
```
### In a loop could take 4 hours, run a job and see example `Job-samtools-split.sh`
```
for file in *.aln.sam.sort.bam ; do samtools idxstats $file | cut -f 1 | grep -w -f AUTOSOMES-chr.txt | xargs samtools view -b $file > ./AUTOSOMES-CHR/$file.AUTOSOMES.bam; done
```
## GENERATE THE GENOTYPES BY LIKELIHOODS
 
We are going to perform a genotype likelihoods in a fast way using ANGSD. This method is the fastes and best approach for variable x coverage in samples. There is some information that will not beeing recovered in the final vcf file with this method. Mapping quality of reads needs to be filtering during the ANGSD run.

### Call de environment ANGSD
```
conda activate ANGSD
```
### Creat a list file of your bam files for each sample with the path. See example in `bam-list-unplaced.txt`
This is an old version GATK like genotype likelihoods by using the following flags
- -doGlf 2: binary glf
- -doMajorMinor 1: Infer major and minor from GL
- -SNP_pval 1e-6: If we are interested in looking at allele frequencies only for sites that are actually variable in our sample.
- -minMaf         0.05        (Remove sites with MAF below) You want to remove variants with major alleles very low that could be sequencing errors. If you expect having variant with very low frequency you could change this value.
- -SNP_pval       1.000000        (Remove sites with a pvalue larger) We can consider assigning as SNPs sites whose estimated allele frequency is above a certain threhsold (e.g. the frequency of a singleton) or whose probability of being variable is above a specified value.
- -doMaf 1: Frequency (fixed major and minor)
- -minInd: you want at least SNP present in 75% of the individuals  75% of 20 individuals is = 15 individuals
- -minMapQ 30: here is where the mapping quality filtering is happening. After that SNPs does not need to keep that info in the vcf file
- -minQ 20: here is the filter for the minimum base quality score. 
```
angsd -GL 2 -out genolike1-greenjay_ZW -nThreads 68 -doGlf 2 -minInd 15 -doMajorMinor 1 -SNP_pval 1e-6 -minMaf 0.05 -doMaf 1 -minMapQ 30 -minQ 20 -bam bam-list-ZW.txt
```
If the above line works well run it in a job. See job example `Job-ANGSD-genotypes.sh`
The output is a beagle format, positions and major allele frequencies. The beagle format need to be formated to a simple vcf file with glactools.
### Format beagle genotypes likelihoods 
This vcf format is simplified, and not much SNP filter can be done after. That step was done during the genotype likelihoods calculations.
Run each step one by one. Check that the output file from one step needs to be the input of the next step.
#### load paths to reference and to  glactools
```
gtool=/PATH/to/glactools
REFPATH=GCA_020740725.1_bCorHaw1.pri.cur_genomic.fai
module load gcc/9.1.0
```
##### step 1
Keep the same name for all files
```
make postion file .pos
zcat genolike1-greenjay.mafs.gz | cut -f 1,2 > genolike1-greenjay.mafs.pos
gzip genolike1-greenjay.mafs.pos
```
##### step 2
Creat glf file format
```
$gtool beagle2glf --fai $REFPATH genolike1-greenjay > genolike1-greenjay-UNPLACED.glf.gz
```
##### step 3
Create acf format
```
$gtool glf2acf genolike1-greenjay-UNPLACED.glf.gz > genolike1-greenjay-UNPLACED.acf.gz
```
##### step 4
create simplified vcf format
```
$gtool glac2vcf genolike1-greenjay-UNPLACED.acf.gz > genolike1-greenjay-UNPLACED.vcf
```
### Explore your vcf file
Look how the header looks like in the vcf file
```
bcftools view --header-only genolike1-greenjay-UNPLACED.vcf
```
### Reheader with bcftools
Reheader because ANGSD place individual names into the file format acoording to the list order in the bam-file-list.txt. Needs to change header name because header in this vcf file is not detected by vcftools. Make a header file `header-order-genotype-names.txt` according to the order in `bam-file-list.txt` and use the name format for the samples you prefere most. But keep it simple without special characters.

```
bcftools reheader -s header-order-genotype-names.txt  genolike1-greenjay-UNPLACED.vcf > newH-genolike1-greenjay-UNPLACED.vcf
```
### Check if the header added looks good
```
bcftools view --header-only newH-genolike1-greenjay-UNPLACED.vcf
```
### Calculate population genomics summary statistic genome-wide
Such as Fst, Pi and relatedness statistic. The last one based on the method of Manichaikul et al., BIOINFORMATICS 2010.
```
vcftools --vcf newH-genolike1-greenjay-UNPLACED.vcf --relatedness2
```
### Create the population map for Fst calculation. 
Creat a `pop-map-north.txt` and `pop-map-south.txt` file per population listing the sample names. The resulting output file has the suffix ".fst". A window size could be use but use the smalles one because we are going to use the spline window technique in R to visualize the Fst or Pi distribution in a manhattan-like plot.
```
vcftools --vcf newH-genolike1-greenjay-UNPLACED.vcf --weir-fst-pop pop-map-north.txt --weir-fst-pop pop-map-south.txt --fst-window-size 1000 --out 1Kkb
```
This Fst estimate from Weir and Cockerham’s 1984 paper. This is the preferred calculation of Fst. The provided file must contain a list of individuals - one individual per line - from the VCF file that correspond to one population. This option can be used multiple times to calculate Fst for more than two populations. These files will also be included as "--keep" options. By default, calculations are done on a per-site basis. The output file has the suffix ".weir.fst".

### Calculate nucleotide diversity per bin of size 1kb
```
vcftools --vcf newH-genolike1-greenjay-UNPLACED.vcf --window-pi-step 1000 --out Pi1Kkb 
```

### Prune SNPs every 10,000kb to select unlinked for PCA and other analyses
```
vcftools --vcf newH-genolike1-greenjay-UNPLACED.vcf --thin 10000 --recode --recode-INFO-all --out newH-genolike1-greenjay-UNPLACED-10Kkb.vcf
```
From this vcf file it is possible to create many formats 
### Create plink format for admixture program and others
```
vcftools --vcf newH-genolike1-greenjay-UNPLACED.vcf --plink 
```
