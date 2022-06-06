#!/bin/bash
#SBATCH -J mpileup           # Job name
#SBATCH -o mpileup.o%j  # Name of stdout output file
#SBATCH -e mpileup.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load samtools/1.5

GB=/scratch/08752/ftermig/ref-genome/GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna
mapfile -t CHR < AUTOSOMES-chr.txt
#
for str in ${CHR[@]}; do for file in *.sorted.bam; do
samtools mpileup -Q 30 -q 20 -u -v \
-f $GB -r $str $file |  
bcftools call -c |  
vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > /scratch/08752/ftermig/Demography/AUTOSOMES/$file.$str.fq
done;done
