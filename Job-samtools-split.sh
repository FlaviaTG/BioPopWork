#!/bin/bash
#SBATCH -J SPLITsamtools           # Job name
#SBATCH -o SPLITsamtools.o%j  # Name of stdout output file
#SBATCH -e SPLITsamtools.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 26:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load samtools/1.5

for file in *.aln.sam.sort.bam ; do samtools idxstats $file | cut -f 1 | grep -w -f AUTOSOMES-chr.txt | xargs samtools view -b $file > ./AUTOSOMES-CHR/$file.AUTOSOMES.bam; done
