#!/bin/bash
#SBATCH -J FastQC           # Job name
#SBATCH -o FastQC.o%j  # Name of stdout output file
#SBATCH -e FastQC.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Bioinformatics-Works	   # Allocation name (req'd if you have more than 1)


module load fastqc/0.11.5

for i in *.fastq.gz; do fastqc $i -o ./QUALITY-reads;done
