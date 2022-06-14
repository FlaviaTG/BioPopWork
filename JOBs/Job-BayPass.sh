#!/bin/bash
#SBATCH -J BayPass           # Job name
#SBATCH -o ADMIX.o%j  # Name of stdout output file
#SBATCH -e ADMIX.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 2:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load gcc/9.1.0

BayPass=/scratch/08752/ftermig/programs/baypass_2.3/sources


$BayPass/g_baypass -nthreads 8 -npop 2 -gfile BayPass-format-test-ZW-thin.txt -efile Scale-var-Bio1_covariate.std -auxmodel -omegafile Scale-var-Bio1_mat_omega.out -outprefix Aux-var-Bio1
