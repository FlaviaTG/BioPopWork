#!/bin/bash
#SBATCH -J bwa_index           # Job name
#SBATCH -o bwa_index.o%j       # Name of stdout output file
#SBATCH -e bwa_index.e%j       # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 01:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Bioinformatics-Works       # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load bwa/0.7.16a

cd /scratch/08752/ftermig/ref-genome/

GB=GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna.gz

bwa index -p CorHaw $GB
