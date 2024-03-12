#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=20GB
#SBATCH --account=bbenayou_34
module purge
module load gcc/8.3.0
/project/bbenayou_34/mcgill/subread-2.0.2-source/bin/featureCounts -O -p -M -a /project/bbenayou_34/mcgill/mm39_refGene.gtf -o EMA_microglia_counts.txt E3WT_17a_rep1_STARAligned.sortedByCoord.out.bam E3WT_17a_rep2_STARAligned.sortedByCoord.out.bam E3WT_17a_rep3_STARAligned.sortedByCoord.out.bam E4WT_17a_rep1_STARAligned.sortedByCoord.out.bam E4WT_17a_rep2_STARAligned.sortedByCoord.out.bam E4WT_17a_rep3_STARAligned.sortedByCoord.out.bam E4WT_17a_rep4_STARAligned.sortedByCoord.out.bam E3WT_CTL_rep1_STARAligned.sortedByCoord.out.bam E3WT_CTL_rep2_STARAligned.sortedByCoord.out.bam E3WT_CTL_rep3_STARAligned.sortedByCoord.out.bam E4WT_CTL_rep1_STARAligned.sortedByCoord.out.bam E4WT_CTL_rep2_STARAligned.sortedByCoord.out.bam E4WT_CTL_rep3_STARAligned.sortedByCoord.out.bam E4WT_CTL_rep4_STARAligned.sortedByCoord.out.bam E4WT_CTL_rep5_STARAligned.sortedByCoord.out.bam
