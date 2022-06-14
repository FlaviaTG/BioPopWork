#!/bin/bash
#SBATCH -J samtoolsD           # Job name
#SBATCH -o samtoolsD.o%j  # Name of stdout output file
#SBATCH -e samtoolsD.e%j  # Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 35:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

module load intel/17.0.4
module load samtools/1.5

# file header  
echo "bam_file total ref"

# loop to count reads in all files and add results to a table  
for bam_file in *.bam  
do  
total=$(samtools depth $bam_file | awk '{sum+=$3} END { print "Average = ",sum/NR}')  
ref=$(samtools view -H $bam_file | grep -P '^@SQ' | cut -f 3 -d ':' | awk '{sum+=$1} END {print sum}')  
echo "$bam_file $total $ref" >> coverage.txt 
done  
# ultimately the data is saved in a space-separated file
