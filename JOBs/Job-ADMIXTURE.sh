#!/bin/bash
#SBATCH -J ADMIX           # Job name
#SBATCH -o ADMIX.o%j  # Name of stdout output file
#SBATCH -e ADMIX.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 24:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Bioinformatics-Works	   # Allocation name (req'd if you have more than 1)

source /home1/08752/ftermig/miniconda3/etc/profile.d/conda.sh
conda activate SUMMARY

for K in 1 2 3 4 5 6 7 8 9 10; do admixture -B2000 -j40 --cv newH-genolike2-greenjay_ZW-test.NOmissing.THIN.bed $K | tee Bootlog${K}.out; done

