#!/bin/bash
#SBATCH -J genotypeANGSD           # Job name
#SBATCH -o genotypeANGSD.o%j  # Name of stdout output file
#SBATCH -e genotypeANGSD.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

source /home1/08752/ftermig/miniconda3/etc/profile.d/conda.sh
conda activate ANGSD

angsd -GL 2 -doBcf 1 -out genolike2-greenjay_ZW -nThreads 68 -doPost 1 -docounts 1 -dogeno 1 -minInd 15 -doMajorMinor 1 -SNP_pval 1e-6 -minMaf 0.05 -doMaf 1 -minMapQ 30 -minQ 20 -bam bam-list-ZW.txt
