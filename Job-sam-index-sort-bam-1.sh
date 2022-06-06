#!/bin/bash
#SBATCH -J picard           # Job name
#SBATCH -o picard-CTAGCGCT.o%j	# Name of stdout output file
#SBATCH -e picard-CTAGCGCT.e%j	# Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 06:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

export OMP_NUM_THREADS=32    # 68 total OpenMP threads (1 per KNL core)

module load intel/17.0.4
module load picard/2.11.0

java -Xmx128g -XX:ParallelGCThreads=32 -jar picard/build/libs/picard.jar SortSam I=sample.CTAGCGCT.aln.sam O=sample.CTAGCGCT.aln.sam.sort.bam SORT_ORDER=coordinate CREATE_INDEX=true
#all in a loop
#for sample in *.aln.sam;do java -Xmx128g -XX:ParallelGCThreads=32 -jar picard/build/libs/picard.jar SortSam I=$sample O=sample.sort.bam SORT_ORDER=coordinate CREATE_INDEX=true;done
