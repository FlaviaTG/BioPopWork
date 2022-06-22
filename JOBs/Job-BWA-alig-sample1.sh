#!/bin/bash
#SBATCH -J bwa_align-CTAGCGCT           # Job name
#SBATCH -o bwa_align-CTAGCGCT.o%j       # Name of stdout output file
#SBATCH -e bwa_align-CTAGCGCT.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Bioinformatics-Works       # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load bwa/0.7.16a

cd /scratch/08752/ftermig/data

ref=/scratch/08752/ftermig/ref-genome
bwa mem -t 30 -M -R '@RG\tID:@A00672:51:HWJMYDSXY\tSM:CTAGCGCT' $ref/CorHaw 01.180603.B03.S1.1.fastq.gz 01.180603.B03.S1.2.fastq.gz > sample.CTAGCGCT.aln.sam
