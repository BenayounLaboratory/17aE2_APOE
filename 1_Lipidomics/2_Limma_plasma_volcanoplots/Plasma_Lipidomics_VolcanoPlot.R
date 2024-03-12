################################################################################
#################### EMA Plasma Lipids  Volcano Plot ###########################
################################################################################

#March 12, 2024

# set working directory
setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/20211027Cassie/Exp3/")

# set options
options(stringsAsFactors = F)

# load packages
require(ggplot2)
library("ggrepel")
library("scales")

################################################################################
#################### APOE4 Control vs APOE3 Control ############################
################################################################################

plot.data <- read.delim("2023-06-16_EMA_Lipidomics_Plasma_Class_VSN_Limma_Plasma_E3vsE4_ALL.txt", header = TRUE)

head(plot.data)
plot.data$newFDR <- -log10(plot.data$adj.P.Val)
plot.data$gene <- rownames(plot.data)
plot.data$label <- "NS"
plot.data$label[plot.data$logFC > 0 & plot.data$adj.P.Val < 0.05] <- "APOE4"
plot.data$label[plot.data$logFC < 0 & plot.data$adj.P.Val < 0.05] <- "APOE3"

head(plot.data)

p <- ggplot(plot.data, aes(logFC, newFDR))
p <- p + geom_point(aes(color = factor(label)), size = 2)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor =
                 element_blank())
p <- p + geom_hline(yintercept = -log10(0.05), color = "gray44",
                    linetype="dashed")
p <- p + scale_colour_manual(values = c("steelblue4","red4", "gray"))
p <- p + xlim(-2,2)
p <- p + ylim(0,5)
p <- p + xlab("\n")
p <- p + ylab("\n")
p <- p + theme(axis.text=element_text(size=15),
               axis.title=element_text(size=15))
p <- p + scale_x_continuous(breaks=c(-2,-1, 0, 1, 2),
                            labels=c("-2","-1", "0", "1", "2"),
                            limits=c(-2,2))
# label the genes
p <- p+geom_text_repel(data=head(plot.data, 10), aes(label=gene))

print(p)
dev.off()
ggsave(filename = "VolcanoPlotControls.pdf", p, height = 5, width = 7)

################################################################################
#################### APOE4 Treated vs APOE3 Control ############################
################################################################################

#clear environment# 

plot.data <- read.delim("2023-06-16_EMA_Lipidomics_Plasma_Class_VSN_Limma_Plasma_E3ConE417a_ALL.txt", header = TRUE)

head(plot.data)
plot.data$newFDR <- -log10(plot.data$adj.P.Val)
plot.data$gene <- rownames(plot.data)
plot.data$label <- "NS"
plot.data$label[plot.data$logFC > 0 & plot.data$adj.P.Val < 0.05] <- "APOE4 Treated"
plot.data$label[plot.data$logFC < 0 & plot.data$adj.P.Val < 0.05] <- "APOE3 Control"

head(plot.data)

p <- ggplot(plot.data, aes(logFC, newFDR))
p <- p + geom_point(aes(color = factor(label)), size = 2)
p <- p + theme_bw()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor =
                 element_blank())
p <- p + geom_hline(yintercept = -log10(0.05), color = "gray44",
                    linetype="dashed")
p <- p + scale_colour_manual(values = c("rosybrown3", "gray"))
p <- p + xlim(-3,3)
p <- p + ylim(0,5)
p <- p + xlab("\n")
p <- p + ylab("\n")
p <- p + theme(axis.text=element_text(size=15),
               axis.title=element_text(size=15))
p <- p + scale_x_continuous(breaks=c(-3,-2,-1, 0, 1, 2,3),
                            labels=c("-3","-2","-1", "0", "1", "2","3"),
                            limits=c(-3,3))
# label the genes
#p <- p+geom_text_repel(data=head(plot.data, 10), aes(label=gene))

print(p)
dev.off()

ggsave(filename = "VolcanoPlotAPOE417avsAPOE3Con.pdf", p, height = 5, width = 7)



#######################
sink(paste0(Sys.Date(),"_Volcanoplot_sessionInfo.txt"))
sessionInfo()
sink()
