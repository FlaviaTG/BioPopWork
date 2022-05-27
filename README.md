# BioPopWork
 
 
## Workshop : fast bioinformatic workflow for population genomics
 
This project is the workflow for a workshop in Bioinformatics for population genomics analyses with resequencing data and the reference genome.
The aim of the workflow is to show how to obtain formats and analyses for fast-and-first exploration of resequencing data of populations. During a paper
discussion on each session of the workshop, we will reflect on alternative methods to the ones used during the on-hand sessions.

In the future this workflow would have the flow to run in one-click script with many plots for biological interpretation as ouputs.

### INPUT : fastq 

Pair-end sequencing data, R1.fq and R2.fq files per samples

The first step is tio check for read quality with using Fred scores by runing fastqc on each file.

#### run the modules or environment
`module load fastqc/0.11.5`
`fastqc 01.180603.B03.S1.1.fastq.gz -o ./QUALITY-reads`
#### In a loop but in a slurm job see bash script for an example `Job-FastqC.sh`
for i in *.fastq.gz; do fastqc $i -o ./QUALITY-reads;done 
