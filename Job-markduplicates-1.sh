#!/bin/bash
#SBATCH -J picard           # Job name
#SBATCH -o picard-CTAGCGCT.o%j	# Name of stdout output file
#SBATCH -e picard-CTAGCGCT.e%j	# Name of stderr error file
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 68               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 18:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-user=ftermignoni@fas.harvard.edu
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH -A Keitt-Lab	   # Allocation name (req'd if you have more than 1)

export OMP_NUM_THREADS=3    # 68 total OpenMP threads (1 per KNL core)

module load intel/17.0.4
module load picard/2.11.0

java -Xmx2g -Xms1g -XX:ParallelGCThreads=3 -jar picard/build/libs/picard.jar MarkDuplicates TMP_DIR=tmp I=sample.CTAGCGCT.aln.sam.sort.bam O=sample.CTAGCGCT.aln.sam.sort.bam.dedup.bam METRICS_FILE=sample.CTAGCGCT.aln.sam.sort.bam.dedup.bam.metrics.txt MAX_RECORDS_IN_RAM=1000000 REMOVE_DUPLICATES=true TAGGING_POLICY=All
#all in a loop
#for sample in *.aln.sam;do java -Xmx68g -XX:ParallelGCThreads=3 --MAX_RECORDS_IN_RAM 1000000 -jar picard/build/libs/picard.jar SortSam I=$sample O=sample.sort.bam SORT_ORDER=coordinate CREATE_INDEX=true;done
