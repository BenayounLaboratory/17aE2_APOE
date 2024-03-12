################################################################################
#################### EMA Microglia Bulk RNA-seq  ###############################
################################################################################

setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/June2023")
options(stringsAsFactors = FALSE)

# load libraries for analysis
library(DESeq2)
library(pheatmap)
library('pvclust')
library('bitops')
library('sva')
library('limma')
library(RColorBrewer)
library(fields)

# 2024-03-11

# read in count matrix
my.ema1 <- read.table('EMA_microglia_counts.txt', sep = "\t", header = T) #26264 genes
my.ema <- my.ema1[,-c(2:6)]

colnames(my.ema) <- c("GeneSymbol",                                          
                      "E3WT_17a_rep1",
                      "E3WT_17a_rep2",
                      "E3WT_17a_rep3",
                      "E4WT_17a_rep1",
                      "E4WT_17a_rep2",
                      "E4WT_17a_rep3",
                      "E4WT_17a_rep4",
                      "E3WT_CTL_rep1",
                      "E3WT_CTL_rep2",
                      "E3WT_CTL_rep3",
                      "E4WT_CTL_rep1",
                      "E4WT_CTL_rep2",
                      "E4WT_CTL_rep3",
                      "E4WT_CTL_rep4",
                      "E4WT_CTL_rep5")

my.ema.2 <- my.ema[c("GeneSymbol"     ,
                     "E3WT_CTL_rep1"  ,
                     "E3WT_CTL_rep2"  ,
                     "E3WT_CTL_rep3"  ,
                     "E3WT_17a_rep1"  ,
                     "E3WT_17a_rep2"  ,
                     "E3WT_17a_rep3"  ,
                     "E4WT_CTL_rep1"  ,
                     "E4WT_CTL_rep2"  ,
                     "E4WT_CTL_rep3"  ,
                     "E4WT_CTL_rep4"  ,
                     "E4WT_CTL_rep5"  ,
                     "E4WT_17a_rep1"  ,
                     "E4WT_17a_rep2"  ,
                     "E4WT_17a_rep3"  ,
                     "E4WT_17a_rep4"  ) ]

my.geno       <- factor(c(rep("E3",6),rep("E4",9))) #tell which genotypes 
my.treatment  <- factor(c(rep("CTL",3),rep("17aE2",3),rep("CTL",5),rep("17aE2",4)))  # tell which treatment 
my.group      <- factor( paste0(my.geno,my.treatment)) # for SVA to avoid overcorrection
        
# RIN
# E3WT_17a _rep1	2.5
# E3WT_17a _rep2	7.3
# E3WT_17a _rep3	7.1
# E4WT_17a _rep1	6.6
# E4WT_17a _rep2	2.5
# E4WT_17a _rep3	8.3
# E4WT_17a _rep4	6.1
# E3WT_control _rep1	7.7
# E3WT_control _rep2	7.6
# E3WT_control _rep3	6.2
# E4WT_control _rep1	6.2
# E4WT_control _rep2	5.9
# E4WT_control _rep3	6.6
# E4WT_control _rep4	2.6
# E4WT_control _rep5	5.5

my.rin        <- c(7.7, # E3WT_control _rep1
                   7.6, # E3WT_control _rep2
                   6.2, # E3WT_control _rep3
                   2.5, # E3WT_17a _rep1
                   7.3, # E3WT_17a _rep2
                   7.1, # E3WT_17a _rep3
                   6.2, # E4WT_control _rep1
                   5.9, # E4WT_control _rep2
                   6.6, # E4WT_control _rep3
                   2.6, # E4WT_control _rep4
                   5.5, # E4WT_control _rep5
                   6.6, # E4WT_17a _rep1
                   2.5, # E4WT_17a _rep2
                   8.3, # E4WT_17a _rep3
                   6.1) # E4WT_17a _rep4
                  

# get the genes with no reads out
my.good <- which(apply(my.ema.2[,-1]>0, 1, sum) >= ncol(my.ema.2)/2) # see deseq2 vignette, need to remove too low genes
my.filtered.matrix <- my.ema.2[my.good,-1] # 18545 genes
rownames(my.filtered.matrix) <- my.ema.2$GeneSymbol[my.good]

################################################################################
# 1. Run SVA to remove unwanted variation
# build design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.matrix ), 
                         group     = my.group ,
                         rin       = my.rin
                         )

# Set null and alternative models (ignore batch)
mod1 = model.matrix(~ group + rin, data = dataDesign)
n.sv.be = num.sv(my.filtered.matrix,mod1,method="be") # 1 surrogate variable

# apply SVAseq algortihm
my.svseq = svaseq(as.matrix(my.filtered.matrix), mod1, n.sv=n.sv.be, constant = 0.1)

# remove SV, preserve genotype
my.clean <- removeBatchEffect(log2(my.filtered.matrix + 0.1), 
                              batch=NULL, 
                              covariates=cbind(my.svseq$sv, my.rin),
                              design=mod1[,1:4])

# delog and round data for DEseq2 processing
my.filtered.sva <- round(2^my.clean-0.1)
write.table(my.filtered.sva, file = paste(Sys.Date(),"EMA_microglia_Cleaned_Counts.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
################################################################################

# 2. DESeq2 on cleaned data
my.outprefix <- paste(Sys.Date(),"EMA_microglia_DESeq2_Analysis",sep="_")


################################################################################
############# A. model genotype and treatment together  ########################
################################################################################

# design matrix
dataDesign = data.frame( row.names = colnames( my.filtered.sva ), 
                         geno      = my.geno,
                         treatment = my.treatment,
                         group = my.group)

# get matrix using age as a modeling covariate
dds <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                              colData   = dataDesign,
                              design    = ~ geno + treatment)

# run DESeq normalizations and export results
dds.deseq <- DESeq(dds)

# plot dispersion
my.disp.out <- paste(my.outprefix,"dispersion_plot.pdf",sep="_")

pdf(my.disp.out)
plotDispEsts(dds.deseq)
dev.off()

# normalized expression value
tissue.cts <- getVarianceStabilizedData(dds.deseq)

# color-code 
my.colors <- rep("steelblue4",15)
my.colors[grep("E3WT_17a",colnames(tissue.cts))] <- "lightskyblue2"
my.colors[grep("E4WT_CTL",colnames(tissue.cts))] <- "red4"
my.colors[grep("E4WT_17a",colnames(tissue.cts))] <- "rosybrown3"

# do MDS analysis
mds.result <- cmdscale(1-cor(tissue.cts,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"MDS_plot.pdf",sep="_")
pdf(my.mds.out)
plot(x, y,
     xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling/ EMA Microglia",
     cex=3, pch = 16, col = my.colors,
     xlim = c(-0.06,0.06),
     ylim = c(-0.06,0.06),
     cex.lab = 1.5,
     cex.axis = 1.5)
legend("topright",c("E3 CTL","E3 17aE2","E4 CTL","E4 17aE2"), col = c("steelblue4","lightskyblue2", "red4","rosybrown3"), bty = 'n', pch = 16)
dev.off()

# PCA analysis
my.pos.var <- apply(tissue.cts,1,var) > 0
my.pca <- prcomp(t(tissue.cts[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot.pdf",sep="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, pch = 16, col = my.colors,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5) 
legend("topright",c("E3 CTL","E3 17aE2","E4 CTL","E4 17aE2"), col = c("steelblue4","lightskyblue2", "red4","rosybrown3"), bty = 'n', pch = 16)
dev.off()


# expression range
pdf(paste(my.outprefix,"_Normalized_counts_boxplot.pdf"))
boxplot(tissue.cts,col=my.colors,cex=0.5,ylab="Log2 DESeq2 Normalized counts", las = 2)  
dev.off()

# check microglia marker expression
genes.purity <- read.table("PurityCheck.txt", header =T)
genes.purity <- genes.purity$GeneName

my.num.purity <- length(genes.purity)
my.heatmap.out <- paste(my.outprefix, "PurityCheck.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Cell types, ", my.num.purity, " genes",sep="")
pheatmap(tissue.cts[genes.purity,],
         cluster_cols = F,
         cluster_rows = F,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = T , scale = "column",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

################################################################################
############# A1. Genotype effect with treatment as covariate  #################
################################################################################

res.geno <- results(dds.deseq, contrast = c("geno","E4","E3")) # FC in E4 over E3
res.geno <- res.geno[!is.na(res.geno$padj),]

### get the heatmap of changes at FDR5; exclude NA
genes.geno  <- rownames(res.geno)[res.geno$padj < 0.05]
my.num.geno <- length(genes.geno) # 924 genes

my.heatmap.out <- paste(my.outprefix,"E3_vsE4_treatmentRegressed_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Genotype significant (FDR<5%), ", my.num.geno, " genes",sep="")
pheatmap(tissue.cts[genes.geno,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

################################################################################
############# A2. Treatment effect with genotype as covariate  #################
################################################################################

res.17a  <- results(dds.deseq, contrast = c("treatment","17aE2","CTL")) # FC in treat vs CTL
res.17a    <- res.17a[!is.na(res.17a$padj),]
genes.17a  <- rownames(res.17a)[res.17a$padj < 0.05]
my.num.17a <- length(genes.17a) # 30 genes

### get the heatmap of changes at FDR5; exclude N
my.heatmap.out <- paste(my.outprefix,"treatment_genoRegressed_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Treatment significant (FDR<5%), ", my.num.17a, " genes",sep="")
pheatmap(tissue.cts[genes.17a,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

# output result tables of combined analysis to text files
my.out.stats.geno <- paste0(my.outprefix,"_APOE_Genotype_all_genes_statistics.txt")
my.out.stats.17a  <- paste0(my.outprefix,"_TREATMENT_all_genes_statistics.txt")
write.table(res.geno, file = my.out.stats.geno , sep = "\t" , row.names = T, quote=F)
write.table(res.17a, file = my.out.stats.17a , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.geno <- paste0(my.outprefix,"_APOE_Genotype_FDR5_genes_statistics.txt")
my.out.fdr5.17a  <- paste0(my.outprefix,"_TREATMENT_FDR5_genes_statistics.txt")
write.table(res.geno[genes.geno,], file = my.out.fdr5.geno, sep = "\t" , row.names = T, quote=F)
write.table(res.17a[genes.17a,], file = my.out.fdr5.17a, sep = "\t" , row.names = T, quote=F)

my.rdata.geno <- paste0(my.outprefix,"_APOE_Genotype_DEseq2_object.RData")
my.rdata.17a  <- paste0(my.outprefix,"_TREATMENT_DEseq2_object.RData")
save(res.geno, file = my.rdata.geno)
save(res.17a , file = my.rdata.17a)


################################################################################
############# B. model treatment in each genotype separately  ##################
################################################################################

# get matrix using Treatment as a modeling covariate in E3
dds.E3 <- DESeqDataSetFromMatrix(countData = my.filtered.sva[,my.geno %in% 'E3'],
                                colData    = dataDesign[my.geno %in% 'E3',],
                                design     = ~ treatment)

# run DESeq normalizations and export results
dds.deseq.E3 <- DESeq(dds.E3)
res.E3 <- results(dds.deseq.E3, contrast = c("treatment","17aE2","CTL")) #FC in E3 17a over E3 Control

# get matrix using Treatment as a modeling covariate in E4
dds.E4 <- DESeqDataSetFromMatrix(countData = my.filtered.sva[,my.geno %in% 'E4'],
                                 colData    = dataDesign[my.geno %in% 'E4',],
                                 design     = ~ treatment)

# run DESeq normalizations and export results
dds.deseq.E4 <- DESeq(dds.E4)
res.E4 <- results(dds.deseq.E4, contrast = c("treatment","17aE2","CTL")) #FC in E4 17a over E4 Control


res.E3    <- res.E3[!is.na(res.E3$padj),]
genes.E3  <- rownames(res.E3)[res.E3$padj < 0.05]
my.num.E3<- length(genes.E3) # 1016 genes

### get the heatmap of changes at FDR5; exclude N
my.heatmap.out <- paste(my.outprefix,"E3_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("APOE3 significant (FDR<5%), ", my.num.E3, " genes",sep="")
pheatmap(tissue.cts[genes.E3,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()


res.E4    <- res.E4[!is.na(res.E4$padj),]
genes.E4  <- rownames(res.E4)[res.E4$padj < 0.05]
my.num.E4<- length(genes.E4) # 37 genes

### get the heatmap of changes at FDR5; exclude N
my.heatmap.out <- paste(my.outprefix,"E4_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("APOE4 significant (FDR<5%), ", my.num.E4, " genes",sep="")
pheatmap(tissue.cts[genes.E4,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

my.rdata.E3 <- paste0(my.outprefix,"_TREATMENT_E3_DEseq2_object.RData")
my.rdata.E4 <- paste0(my.outprefix,"_TREATMENT_E4_DEseq2_object.RData")
save(res.E3, file = my.rdata.E3)
save(res.E4, file = my.rdata.E4)

my.out.stats.E3 <- paste0(my.outprefix,"_APOE3_all_genes_statistics.txt")
my.out.stats.E4  <- paste0(my.outprefix,"_APOE4_all_genes_statistics.txt")
write.table(res.E3, file = my.out.stats.E3 , sep = "\t" , row.names = T, quote=F)
write.table(res.E4, file = my.out.stats.E4 , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.E3 <- paste0(my.outprefix,"_APOE3_FDR5_genes_statistics.txt")
my.out.fdr5.E4  <- paste0(my.outprefix,"_APOE4_FDR5_genes_statistics.txt")
write.table(res.E3[genes.E3,], file = my.out.fdr5.E3, sep = "\t" , row.names = T, quote=F)
write.table(res.E4[genes.E4,], file = my.out.fdr5.E4, sep = "\t" , row.names = T, quote=F)


################################################################################
##################### C. Log fold change comparison  ###########################
################################################################################

res.E3 <- as.data.frame(res.E3)
res.E4 <- as.data.frame(res.E4)

colnames(res.E3) <- paste(colnames(res.E3),"E3",sep = "_")
colnames(res.E4) <- paste(colnames(res.E4),"E4",sep = "_")

my.Merged.E34 <- merge(res.E3,res.E4, by = 0 , all = TRUE)
my.Merged.E34 <- my.Merged.E34[!is.na(my.Merged.E34$padj_E3),]
my.Merged.E34 <- my.Merged.E34[!is.na(my.Merged.E34$padj_E4),]

my.spear.cor <- cor.test(my.Merged.E34$log2FoldChange_E3,my.Merged.E34$log2FoldChange_E4, method = 'spearman', exact=FALSE)
my.rho <- signif(my.spear.cor$estimate,3)

#### commonly regulated genes
my.E3_E4.5  <- bitAnd(my.Merged.E34$padj_E3 < 0.05, my.Merged.E34$padj_E4 < 0.05) > 0
my.E3_E4.10 <- bitAnd(my.Merged.E34$padj_E3 < 0.1 , my.Merged.E34$padj_E4 < 0.1) > 0

pdf(paste(Sys.Date(),"Microglia_COMBINED_treatment_E3_vs_E4_FC_scatterplot_FDR5_10.pdf", sep = "_"))
smoothScatter(my.Merged.E34$log2FoldChange_E3,my.Merged.E34$log2FoldChange_E4,
              xlim = c(-8,8), ylim = c(-8,8),
              xlab = "log2(FC) in APOE3 with 17aE2",
              ylab = "log2(FC) in APOE4 with 17aE2",
              main = "Microglia")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-7, 8 , paste("Rho ~ ",my.rho), pos = 4)
text(-7, 7, paste("p ~ ",signif(my.spear.cor$p.value,3)), pos = 4)
points(my.Merged.E34$log2FoldChange_E3[my.E3_E4.10],my.Merged.E34$log2FoldChange_E4[my.E3_E4.10], cex = 0.5, pch = 16, col = "gold1")
points(my.Merged.E34$log2FoldChange_E3[my.E3_E4.5],my.Merged.E34$log2FoldChange_E4[my.E3_E4.5], cex = 0.5, pch = 16, col = "orange3")
dev.off()

write.table(my.Merged.E34[my.E3_E4.10,], file = paste(Sys.Date(),"Microglia_COMBINED_treatment_E3_vs_E4Agreement_FDR10_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)


#### divergently regulated genes
my.E3_notE4 <- bitAnd(my.Merged.E34$padj_E3 < 0.05, my.Merged.E34$padj_E4 > 0.1) > 0
my.E4_notE3 <- bitAnd(my.Merged.E34$padj_E4 < 0.05, my.Merged.E34$padj_E3 > 0.1) > 0

pdf(paste(Sys.Date(),"Microglia_COMBINED_treatment_E3_vs_E4_FC_scatterplot_divergent_FDR5.pdf", sep = "_"))
smoothScatter(my.Merged.E34$log2FoldChange_E3,my.Merged.E34$log2FoldChange_E4,
              xlim = c(-8,8), ylim = c(-8,8),
              xlab = "log2(FC) in APOE3 with 17aE2",
              ylab = "log2(FC) in APOE4 with 17aE2",
              main = "Microglia")
abline(0,1, col = "grey", lty = "dashed")
abline(h = 0, col = "red", lty = "dashed")
abline(v = 0, col = "red", lty = "dashed")
text(-7, 8 , paste("Rho ~ ",my.rho), pos = 4)
text(-7, 7, paste("p ~ ",signif(my.spear.cor$p.value,3)), pos = 4)
points(my.Merged.E34$log2FoldChange_E3[my.E3_notE4],my.Merged.E34$log2FoldChange_E4[my.E3_notE4], cex = 0.5, pch = 16, col = "gray44")
points(my.Merged.E34$log2FoldChange_E3[my.E4_notE3],my.Merged.E34$log2FoldChange_E4[my.E4_notE3], cex = 0.5, pch = 16, col = "dodgerblue3")
legend("topright",
       c(paste("E3_notE4,", sum(my.E3_notE4)),paste("E4_notE3,", sum(my.E4_notE3))),
       col = c("gray44","dodgerblue3"), pch = 16, pt.cex = 0.5,bty = 'n')
dev.off()

write.table(my.Merged.E34[my.E3_notE4,], file = paste(Sys.Date(),"Microglia_COMBINED_treatment_E3notE4_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(my.Merged.E34[my.E4_notE3,], file = paste(Sys.Date(),"Microglia_COMBINED_treatment_E4notE3_FDR5_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

write.table(my.Merged.E34, file = paste(Sys.Date(),"Microglia_COMBINED_treatment_Genotype_Separated_Merged_Table_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(res.E3   , file = paste(Sys.Date(),"Microglia_COMBINED_treatment_E3_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)
write.table(res.E4   , file = paste(Sys.Date(),"Microglia_COMBINED_treatment_E4_ALL_genes_statistics.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

#non-linear relationship, use individual comparisons instead of combined model for analyses 

################################################################################
############# D. model control or treatment separately  ########################
################################################################################

################################################################################
############# D1. Genotype effect between control groups  ######################
################################################################################ 

# get matrix only looking at controls
dds.ctl <- DESeqDataSetFromMatrix(countData = my.filtered.sva[,my.treatment %in% 'CTL'],
                                 colData    = dataDesign[my.treatment %in% 'CTL',],
                                 design     = ~ geno)

# run DESeq normalizations and export results
dds.deseq.ctl <- DESeq(dds.ctl)
res.ctl <- results(dds.deseq.ctl, contrast = c("geno","E4","E3"))
res.ctl    <- res.ctl[!is.na(res.ctl$padj),]
genes.ctl  <- rownames(res.ctl)[res.ctl$padj < 0.05]
my.num.ctl <- length(genes.ctl) # 2819

my.heatmap.out <- paste(my.outprefix,"E3_vsE4_controlsONLY_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Genotype significant (FDR<5%), ", my.num.ctl, " genes",sep="")
pheatmap(tissue.cts[genes.ctl,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

################################################################################
############# D2. Treatment effect between treated groups  #####################
################################################################################ 

# get matrix only looking at treated groups
dds.treat <- DESeqDataSetFromMatrix(countData = my.filtered.sva[,my.treatment %in% '17aE2'],
                                 colData    = dataDesign[my.treatment %in% '17aE2',],
                                 design     = ~ geno)

# run DESeq normalizations and export results
dds.deseq.treat <- DESeq(dds.treat)
res.treat <- results(dds.deseq.treat, contrast = c("geno","E4","E3"))
res.treat    <- res.treat[!is.na(res.treat$padj),]
genes.treat  <- rownames(res.treat)[res.treat$padj < 0.05]
my.num.treat <- length(genes.treat) # 251

my.heatmap.out <- paste(my.outprefix,"E3_vsE4_treatedONLY_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("Treated groups significant (FDR<5%), ", my.num.treat, " genes",sep="")
pheatmap(tissue.cts[genes.treat,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

my.out.stats.ctl <- paste0(my.outprefix,"_APOE_Genotype_ControlsOnly_all_genes_statistics.txt")
my.out.stats.treat  <- paste0(my.outprefix,"_TREATMENT_only_all_genes_statistics.txt")
write.table(res.ctl, file = my.out.stats.ctl , sep = "\t" , row.names = T, quote=F)
write.table(res.treat, file = my.out.stats.treat , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.ctl <- paste0(my.outprefix,"_APOE_Genotype_ControlsOnly_FDR5_genes_statistics.txt")
my.out.fdr5.treat  <- paste0(my.outprefix,"_TREATMENT_only_FDR5_genes_statistics.txt")
write.table(res.ctl[genes.ctl,], file = my.out.fdr5.ctl, sep = "\t" , row.names = T, quote=F)
write.table(res.treat[genes.treat,], file = my.out.fdr5.treat, sep = "\t" , row.names = T, quote=F)

my.rdata.ctl <- paste0(my.outprefix,"_ControlsONLY_DEseq2_object.RData")
my.rdata.treat <- paste0(my.outprefix,"_TREATMENTONLY_DEseq2_object.RData")
save(res.ctl, file = my.rdata.ctl)
save(res.treat, file = my.rdata.treat)

################################################################################
#################### E. APOE4 17a vs APOE3 Control  ############################
################################################################################
# get matrix looking at groups separately
dds.E4treatE3con <- DESeqDataSetFromMatrix(countData = my.filtered.sva,
                              colData   = dataDesign,
                              design    = ~ group)
# run DESeq normalizations and export results
dds.deseq.E4treatE3con <- DESeq(dds.E4treatE3con)
res.E4treatE3con <- results(dds.deseq.E4treatE3con, contrast = c("group","E417aE2","E3CTL"))
res.E4treatE3con    <- res.E4treatE3con[!is.na(res.E4treatE3con$padj),]
genes.E4treatE3con  <- rownames(res.E4treatE3con)[res.E4treatE3con$padj < 0.05] #652
my.num.E4treatE3con <- length(genes.E4treatE3con) # 652 genes

my.heatmap.out <- paste(my.outprefix,"APOE417avsAPOE3Control_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("APOE4 17aE2 vs APOE3 Control significant (FDR<5%), ", my.num.E4treatE3con, " genes",sep="")
pheatmap(tissue.cts[genes.E4treatE3con,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

my.out.stats.E4treatE3con <- paste0(my.outprefix,"_APOE4treatvsAPOE3con_all_genes_statistics.txt")
write.table(res.E4treatE3con, file = my.out.stats.E4treatE3con , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.E4treatE3con <- paste0(my.outprefix,"_APOE4treatvsAPOE3con_FDR5_genes_statistics.txt")
write.table(res.E4treatE3con[genes.E4treatE3con,], file = my.out.fdr5.E4treatE3con, sep = "\t" , row.names = T, quote=F)

my.rdata.E4treatE3con <- paste0(my.outprefix,"_APOE4treatvsAPOE3con_DEseq2_object.RData")
save(res.E4treatE3con, file = my.rdata.E4treatE3con)


################################################################################
#################### F. APOE3 17a vs APOE4 Control  ############################
################################################################################

# run DESeq normalizations and export results

res.E3treatE4con <- results(dds.deseq.E4treatE3con, contrast = c("group","E317aE2","E4CTL"))
res.E3treatE4con    <- res.E3treatE4con[!is.na(res.E3treatE4con$padj),]
genes.E3treatE4con  <- rownames(res.E3treatE4con)[res.E3treatE4con$padj < 0.05] #829 genes
my.num.E3treatE4con <- length(genes.E3treatE4con) 

my.heatmap.out <- paste(my.outprefix,"APOE317avsAPOE4Control_Heatmap_FDR5.pdf", sep = "_")

pdf(my.heatmap.out, onefile = F)
my.heatmap.title <- paste("APOE3 17aE2 vs APOE4 Control significant (FDR<5%), ", my.num.E3treatE4con, " genes",sep="")
pheatmap(tissue.cts[genes.E3treatE4con,],
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         show_rownames = F, scale="row",
         main = my.heatmap.title, 
         cellwidth = 15,
         border = NA)
dev.off()

my.out.stats.new <- paste0(my.outprefix,"_APOE3treatvsAPOE4con_all_genes_statistics.txt")
write.table(res.E3treatE4con, file = my.out.stats.new , sep = "\t" , row.names = T, quote=F)

my.out.fdr5.new <- paste0(my.outprefix,"_APOE3treatvsAPOE4con_FDR5_genes_statistics.txt")
write.table(res.E3treatE4con[genes.E3treatE4con,], file = my.out.fdr5.new, sep = "\t" , row.names = T, quote=F)

my.rdata.E3treatE4con <- paste0(my.outprefix,"_APOE3treatvsAPOE4con_DEseq2_object.RData")
save(res.E3treatE4con, file = my.rdata.E3treatE4con)

#######################
sink(file = paste(Sys.Date(),"Microglia_COMBINED_treatment_RNAseq_analysis_session_Info.txt", sep =""))
sessionInfo()
sink()



