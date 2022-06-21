#!/bin/bash
#SBATCH -J samtools           # Job name
#SBATCH -o samtools.o%j  # Name of stdout output file
#SBATCH -e samtools.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 15:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Bioinformatics-Works	   # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load samtools/1.5

# file header  
echo "bam_file total mapped unmapped good"

# loop to count reads in all files and add results to a table  
for bam_file in *.bam  
do  
total=$(samtools view -c $bam_file)  
mapped=$(samtools view -c -F 4 $bam_file)  
unmapped=$(samtools view -c -f 4 $bam_file)
good=$(samtools view -c -q 20 $bam_file)
echo "$bam_file $total $mapped $unmapped $good" >> MappedReads.txt 
done  
# ultimately the data is saved in a space-separated file
