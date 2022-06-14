#!/bin/bash
#SBATCH -J BayeScan           # Job name
#SBATCH -o BayeScan.o%j  # Name of stdout output file
#SBATCH -e BayeScan.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 6:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load gcc/9.1.0
export OMP_NUM_THREADS=68

SCAN=/scratch/08752/ftermig/programs/BayeScan2.1/binaries
$SCAN/BayeScan2.1_linux64bits newH-genolike2-greenjay_ZW-test.NOmissing.THIN.GESTE -n 5000 -nbp 20 -pilot 5000 -burn 50000 -pr_odds 100 -threads 40 -out_freq -od ./BAYESCAN/
