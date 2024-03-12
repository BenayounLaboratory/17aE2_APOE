#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=60GB
#SBATCH --account=bbenayou_34

module purge
module load gcc/8.3.0
module load star/2.7.0e


WorkDir="/project/bbenayou_34/mcgill"

cd $WorkDir;

export GENIN="/project/bbenayou_34/mcgill/index"

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA1_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA1_R2_paired.fq.gz --readFilesCommand zcat --alignEndsProtrude 15 ConcordantPair --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --outFileNamePrefix E3WT_17a_rep1_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA2_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA2_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E3WT_CTL_rep1_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA3_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA3_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_17a_rep1_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA4_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA4_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_CTL_rep1_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA6_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA6_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E3WT_CTL_rep2_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA7_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA7_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_17a_rep2_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA8_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA8_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_CTL_rep2_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA10_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA10_R2_paired.fq.gz --readFilesCommand zcat --runThreadN  2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E3WT_CTL_rep3_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA11_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA11_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_17a_rep3_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA12_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA12_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_CTL_rep3_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA13_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA13_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E3WT_17a_rep2_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA16_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA16_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_CTL_rep4_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA17_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA17_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E3WT_17a_rep3_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA18_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA18_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_17a_rep4_STAR

STAR --genomeDir $GENIN --readFilesIn /project/bbenayou_34/mcgill/FASTQ/EMA19_R1_paired.fq.gz /project/bbenayou_34/mcgill/FASTQ/EMA19_R2_paired.fq.gz --readFilesCommand zcat --runThreadN 2 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3 --outFilterMismatchNmax 2 --outFilterMultimapNmax 50 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM SortedByCoordinate --alignEndsProtrude 15 ConcordantPair --outFileNamePrefix E4WT_CTL_rep5_STAR



