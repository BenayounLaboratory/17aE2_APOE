## README ##

Code for the manuscript:

Cassandra J. McGill, Amy Christensen, Wenjie Qian, Max Thorwald, Jose Godoy Lugo, Sara Namvari, Olivia S. White, Caleb Finch, Bérénice A. Benayoun, and Christian J Pike.
"Protection against APOE4-associated aging phenotypes with the longevity-promoting intervention 17a-estradiol"
(biorXiv preprint at https://www.biorxiv.org/content/10.1101/2024.03.12.584678v1)

The code is arranged by major analysis carried out:

	#1_Lipidomics
              - 1_Limma_analysis_plasma
                   o Run_Limma_Plasma_Lipidomics.R
                   o Plasma_Lipidomics_byClass.R
              - 2_Limma_analysis_cortex
                   o Run_Limma_Cortex_Lipidomics.R
              - 3_LION_bubbleplot
                   o LION_Plasma_Bubbleplots.R
	#2_Microglia_RNAseq
              - 1_STAR_mapping
                   o map_with_star.sh
              - 2_Read_counting
                   o count_reads.sh
              - 3_DESeq2
                   o Run_DESeq2_microglia_RNAseq.R
                   o DESeq2_VolcanoPlot_Microglia_RNAseq.R
              - 4_GSEA
                   o Run_GSEA_Microglia_RNAseq.R
                   o GSEA_Bubbleplot_Microglia_RNAseq.R

All sequencing data was deposited to SRA under accession PRJNA1078754. 
Lipidomics and DNAge data was were deposited to FigShare DOI:10.6084/m9.figshare.25346143 and DOI:10.6084/m9.figshare.25346143.


