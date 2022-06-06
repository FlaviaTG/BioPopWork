#!/bin/bash
#SBATCH -J gdsFormat           # Job name
#SBATCH -o gdsFormat.o%j  # Name of stdout output file
#SBATCH -e gdsFormat.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 1:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)


module load intel/18.0.0  impi/18.0.0
module load Rstats/3.5.1

R CMD BATCH --quiet --no-restore --no-save gds-Format.R outputfile_SNPrelate

