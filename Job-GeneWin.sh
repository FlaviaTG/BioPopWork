#!/bin/bash
#SBATCH -J Genewin           # Job name
#SBATCH -o Genewin.o%j  # Name of stdout output file
#SBATCH -e Genewin.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 02:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load intel/18.0.0  impi/18.0.0
module load Rstats/3.5.1


for i in GeneWin-CHR*.R;do R CMD BATCH --quiet --no-restore --no-save $i outputfile_GeneWin$i;done
