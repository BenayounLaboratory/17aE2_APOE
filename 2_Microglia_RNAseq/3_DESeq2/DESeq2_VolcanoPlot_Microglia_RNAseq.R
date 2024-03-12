################################################################################
#################### EMA Microglia GSEA Volcano Plot ###########################
################################################################################

# set working directory
setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/June2023")

# set options
options(stringsAsFactors = F)

# load packages
require(ggplot2)
library("ggrepel")
library("scales")

################################################################################
#################### APOE4 Control vs APOE3 Control ############################
################################################################################

plot.data <- read.table("2023-06-19_EMA_microglia_DESeq2_Analysis_APOE_Genotype_ControlsOnly_all_genes_statistics.txt", header = TRUE)

head(plot.data)
plot.data$newFDR <- -log10(plot.data$padj)
plot.data$gene <- rownames(plot.data)
plot.data$label <- "NS"
plot.data$label[plot.data$log2FoldChange > 0 & plot.data$padj < 0.05] <- "APOE4"
plot.data$label[plot.data$log2FoldChange < 0 & plot.data$padj < 0.05] <- "APOE3"

head(plot.data)

p <- ggplot(plot.data, aes(log2FoldChange, newFDR))
p <- p + geom_point(aes(color = factor(label)), size = 2)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor =
                 element_blank())
p <- p + geom_hline(yintercept = -log10(0.01), color = "gray44",
                    linetype="dashed")
p <- p + scale_colour_manual(values = c("steelblue4","red4", "gray"))
p <- p + xlim(-10,10)
p <- p + ylim(0,60)
p <- p + xlab("\n")
p <- p + ylab("\n")
p <- p + theme(axis.text=element_text(size=15),
               axis.title=element_text(size=15))
p <- p + scale_x_continuous(breaks=c(-10,-5, 0, 5, 10),
                            labels=c("-10","-5", "0", "5", "10"),
                            limits=c(-10,10))
# label the genes
#p <- p+geom_text_repel(data=head(plot.data, 50), aes(label=gene))

print(p)

ggsave(filename = "VolcanoPlotControls.pdf", p, height = 7, width = 10)

################################################################################
#################### APOE4 Treated vs APOE3 Control ############################
################################################################################

plot.data <- read.table("2023-06-19_EMA_microglia_DESeq2_Analysis_APOE4treatvsAPOE3con_all_genes_statistics.txt", header = TRUE)

head(plot.data)
plot.data$newFDR <- -log10(plot.data$padj)
plot.data$gene <- rownames(plot.data)
plot.data$label <- "NS"
plot.data$label[plot.data$log2FoldChange > 0 & plot.data$padj < 0.05] <- "APOE4 Treated"
plot.data$label[plot.data$log2FoldChange < 0 & plot.data$padj < 0.05] <- "APOE3 Control"

head(plot.data)

p <- ggplot(plot.data, aes(log2FoldChange, newFDR))
p <- p + geom_point(aes(color = factor(label)), size = 2)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor =
                 element_blank())
p <- p + geom_hline(yintercept = -log10(0.01), color = "gray44",
                    linetype="dashed")
p <- p + scale_colour_manual(values = c("steelblue4","rosybrown3","gray"))
p <- p + xlim(-10,10)
p <- p + ylim(0,60)
p <- p + xlab("\n")
p <- p + ylab("\n")
p <- p + theme(axis.text=element_text(size=15),
               axis.title=element_text(size=15))
p <- p + scale_x_continuous(breaks=c(-10,-5, 0, 5, 10),
                            labels=c("-10","-5", "0", "5", "10"),
                            limits=c(-10,10))

# label the genes

#p <- p+geom_text_repel(data=head(plot.data, 0), aes(label=gene))

print(p)

ggsave(filename = "VolcanoPlotAPOE417avsAPOE3Con.pdf", p, height = 7, width = 10)

################################################################################
#################### APOE3 Control vs APOE3 Treated ############################
################################################################################

plot.data <- read.table("2023-06-19_EMA_microglia_DESeq2_Analysis_APOE3_all_genes_statistics.txt", header = TRUE)

head(plot.data)
plot.data$newFDR <- -log10(plot.data$padj)
plot.data$gene <- rownames(plot.data)
plot.data$label <- "NS"
plot.data$label[plot.data$log2FoldChange > 0 & plot.data$padj < 0.05] <- "APOE3 Treated"
plot.data$label[plot.data$log2FoldChange < 0 & plot.data$padj < 0.05] <- "APOE3 Control"

head(plot.data)

p <- ggplot(plot.data, aes(log2FoldChange, newFDR))
p <- p + geom_point(aes(color = factor(label)), size = 2)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor =
                 element_blank())
p <- p + geom_hline(yintercept = -log10(0.01), color = "gray44",
                    linetype="dashed")
p <- p + scale_colour_manual(values = c("steelblue4","lightskyblue2","gray"))
p <- p + xlim(-10,10)
p <- p + ylim(0,60)
p <- p + xlab("\n")
p <- p + ylab("\n")
p <- p + theme(axis.text=element_text(size=15),
               axis.title=element_text(size=15))
p <- p + scale_x_continuous(breaks=c(-10,-5, 0, 5, 10),
                            labels=c("-10","-5", "0", "5", "10"),
                            limits=c(-10,10))
# label the genes

#p <- p+geom_text_repel(data=head(plot.data, 500), aes(label=gene))

print(p)

ggsave(filename = "VolcanoPlotAPOE3.pdf", p, height = 7, width = 10)




################################################################################
#################### APOE4 Control vs APOE4 Treated ############################
################################################################################
plot.data <- read.table("2023-06-19_EMA_microglia_DESeq2_Analysis_APOE4_all_genes_statistics.txt", header = TRUE)

head(plot.data)
plot.data$newFDR <- -log10(plot.data$padj)
plot.data$gene <- rownames(plot.data)
plot.data$label <- "NS"
plot.data$label[plot.data$log2FoldChange > 0 & plot.data$padj < 0.05] <- "APOE4 Treated"
plot.data$label[plot.data$log2FoldChange < 0 & plot.data$padj < 0.05] <- "APOE4 Control"

head(plot.data)

p <- ggplot(plot.data, aes(log2FoldChange, newFDR))
p <- p + geom_point(aes(color = factor(label)), size = 2)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor =
                 element_blank())
p <- p + geom_hline(yintercept = -log10(0.01), color = "gray44",
                    linetype="dashed")
p <- p + scale_colour_manual(values = c("red4","rosybrown3","gray"))
p <- p + xlim(-10,10)
p <- p + ylim(0,60)
p <- p + xlab("\n")
p <- p + ylab("\n")
p <- p + theme(axis.text=element_text(size=15),
               axis.title=element_text(size=15))
p <- p + scale_x_continuous(breaks=c(-10,-5, 0, 5, 10),
                            labels=c("-10","-5", "0", "5", "10"),
                            limits=c(-10,10))
# label the genes

#p <- p+geom_text_repel(data=head(plot.data, 500), aes(label=gene))

print(p)

ggsave(filename = "VolcanoPlotAPOE4.pdf", p, height = 7, width = 10)





#######################
sink(paste0(Sys.Date(),"_Volcanoplot_sessionInfo.txt"))
sessionInfo()
sink()
