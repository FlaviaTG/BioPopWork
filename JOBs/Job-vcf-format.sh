#!/bin/bash
#SBATCH -J glacFormat           # Job name
#SBATCH -o glacFormat.o%j  # Name of stdout output file
#SBATCH -e glacFormat.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 4:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

gtool=/scratch/08752/ftermig/programs/glactools/glactools
REFPATH=GCA_020740725.1_bCorHaw1.pri.cur_genomic.fna.fai
module load gcc/9.1.0


#$gtool beagle2glf --fai $REFPATH genolike1-greenjay_ZW > genolike1-greenjay_ZW.glf.gz
#$gtool glf2acf genolike1-greenjay_ZW.glf.gz > genolike1-greenjay_ZW.acf.gz
$gtool glac2vcf genolike1-greenjay_ZW.acf.gz > genolike1-greenjay_ZW.vcf
