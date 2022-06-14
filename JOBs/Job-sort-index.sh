#!/bin/bash
#SBATCH -J Sortsamtools           # Job name
#SBATCH -o Sortsamtools.o%j  # Name of stdout output file
#SBATCH -e Sortsamtools.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 5:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load samtools/1.5


#for i in *bam; do samtools sort $i -o $i.sorted.bam;done
for i in *sorted.bam; do samtools index $i;done
