#!/bin/bash
#SBATCH -J psmc           # Job name
#SBATCH -o psmc.o%j  # Name of stdout output file
#SBATCH -e psmc.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 48:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Bioinformatics-Works	   # Allocation name (req'd if you have more than 1)

source /home1/08752/ftermig/miniconda3/etc/profile.d/conda.sh
conda activate SUMMARY

for i in *.consensus.fq; do fq2psmcfa $i > $i.psmcfa;done
