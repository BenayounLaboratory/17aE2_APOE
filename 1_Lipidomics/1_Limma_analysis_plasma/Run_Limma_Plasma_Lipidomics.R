################################################################################
#################### EMA Plasma Lipidomics Species  ############################
################################################################################

# Updated March 7, 2024

#set working directory
setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/20211027Cassie/Exp3")

# set options
options(stringsAsFactors = F)

#load data and transpose # starting with 841 lipids
counts.lipid <- read.delim('20211027_Exp3_EMA_plasma.txt', header = T, check.names = FALSE, row.names=1)
lipid.df <- as.data.frame(t(counts.lipid))


#remove lipids not represented in all samples #565 lipids left
cleaned.lipid <-na.omit(lipid.df) 

#export cleaned and reorganized lipid file for downstream use
write.table(cleaned.lipid, file = paste(Sys.Date(),"Lipidomics_Plasma_Class_cleaned.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

#the next lines are optional if you want to keep in more lipids

#remove lipids that aren't represented in at least 75% of the samples
#lipid.df[is.na(lipid.df)] <- 1
#lipid.df <- log2(lipid.df)
#lipid.df[lipid.df == 0] <- NA
#cleaned.lipid <- lipid.df[rowSums(!is.na(lipid.df)) > 0.75*ncol(lipid.df),]


#load necessary libraries

library(vsn)
library(limma)
library(pheatmap)


my.outprefix <- paste(Sys.Date(),"EMA_Lipidomics_Plasma_Class_VSN",sep="_")


################################################################################
#########################  A. Preprocess data ################################## 
################################################################################

# 1. read raw data and set column order
my.lipid.FA.data.1 <- cleaned.lipid
my.E3.ctl    <- grep("E3WT_ctl", colnames(my.lipid.FA.data.1))
my.E3.17a    <- grep("E3WT_17a", colnames(my.lipid.FA.data.1))
my.E4.ctl    <- grep("E4WT_ctl",colnames(my.lipid.FA.data.1))
my.E4.17a    <- grep("E4WT_17a",colnames(my.lipid.FA.data.1))

# subset and reorder
my.lipid.FA.data <- my.lipid.FA.data.1[,c(my.E3.ctl,my.E3.17a,my.E4.ctl,my.E4.17a)]

# 2. perform variance stabilizing normalization
my.lipid.data.vsn <- as.data.frame(normalizeVSN(as.matrix(my.lipid.FA.data)))
rownames(my.lipid.data.vsn) <- rownames(my.lipid.FA.data)
colnames(my.lipid.data.vsn) <- colnames(my.lipid.FA.data)

#export normalized values
write.table(my.lipid.data.vsn, file = paste(Sys.Date(),"Lipidomics_Plasma_Class_cleaned_VSNnormalized.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# color-code 
my.colors <- rep("steelblue4",16)
my.colors[grep("E3WT_17a",colnames(my.lipid.data.vsn))] <- "lightskyblue2"
my.colors[grep("E4WT_ctl",colnames(my.lipid.data.vsn))] <- "red4"
my.colors[grep("E4WT_17a",colnames(my.lipid.data.vsn))] <- "rosybrown3"

#export to PDF
pdf(paste(Sys.Date(),"Plasma_EMA_Lipidomics_VSN_norm_Class_boxplots.pdf",sep = "_"), width = 8, height = 8)
boxplot(my.lipid.data.vsn, las = 2, col = my.colors)
par(mfrow = c(1,1))
dev.off()

################################################################################
########################  B. Run analysis    ###################################
################################################################################


# 1. do MDS and PCA analysis

# MDS analysis
mds.result <- cmdscale(1-cor(my.lipid.data.vsn,method="spearman"), k = 2, eig = FALSE, add = FALSE, x.ret = FALSE)
x <- mds.result[, 1]
y <- mds.result[, 2]

my.mds.out <- paste(my.outprefix,"MDS_plot.pdf", sep ="_")
pdf(my.mds.out)
plot(x, y, xlab = "MDS dimension 1", ylab = "MDS dimension 2",
     main="Multi-dimensional Scaling",
     cex=3,
     pch = 16, 
     col = my.colors,
     cex.lab = 1.5,
     cex.axis = 1.5)
dev.off()

# PCA analysis
my.pos.var <- apply(my.lipid.data.vsn,1,var) > 0
my.pca <- prcomp(t(my.lipid.data.vsn[my.pos.var,]),scale = TRUE)
x <- my.pca$x[,1]
y <- my.pca$x[,2]

my.summary <- summary(my.pca)

my.pca.out <- paste(my.outprefix,"PCA_plot.pdf", sep ="_")
pdf(my.pca.out)
plot(x,y,
     cex=3, 
     col= my.colors, 
     pch=16,
     xlab = paste('PC1 (', round(100*my.summary$importance[,1][2],1),"%)", sep=""),
     ylab = paste('PC2 (', round(100*my.summary$importance[,2][2],1),"%)", sep=""),
     cex.lab = 1.5,
     cex.axis = 1.5) 
dev.off()



# 2. do limma analysis 
my.genotype <- c(rep("E3WT",8),rep("E4WT",10)) #genotype with treatment regressed out
my.treat <- c(rep("CTL",4),rep("17a",4),rep("CTL",5),rep("17a",5)) #treatment with genotype regressed out
my.con  <- c(rep("E3WT",4),rep("E4WT",5)) #for controls only
my.E3.treat  <- c(rep("CTL",4),rep("17a",4)) #for E3 only 
my.E4.treat <- c(rep("CTL",5),rep("17a",5)) #for E4 only
my.17a <- c(rep("E3WT",4),rep("E4WT",5)) # for treatment only 
my.E3conE417a <- c(rep("E3WT",4), rep("E4WT",5)) # for E3 control vs E4 treated
my.E317aE4con <- c(rep("E3WT",4), rep("E4WT",5)) # for E3 treated vs E4 control

# make individual comparisons

my.lipid.E3.data <- subset(my.lipid.data.vsn, select = -c(9:18))
my.lipid.E4.data <- subset(my.lipid.data.vsn, select = -c(1:8))
my.lipid.E3vsE4.data <- subset(my.lipid.data.vsn, select = -c(5:8,14:18))
my.lipid.E3vsE417a.data <- subset(my.lipid.data.vsn, select = -c(1:4,9:13))
my.lipid.E3conE417a.data <- subset(my.lipid.data.vsn, select =-c(5:13))
my.lipid.E317aE4con.data <- subset(my.lipid.data.vsn, select = -c(1:4,14:18))

# build design matrix
dataDesign = data.frame(row.names = colnames(my.lipid.data.vsn),
                        genotype= my.genotype,
                        treatment = my.treat)
dataDesignCon = data.frame(row.names =colnames(my.lipid.E3vsE4.data),
                           con = my.con)
dataDesignE3 = data.frame(row.names =colnames(my.lipid.E3.data),
                          E3treat = my.E3.treat)
dataDesignE4 = data.frame(row.names =colnames(my.lipid.E4.data),
                          E4treat = my.E4.treat)
dataDesign17a = data.frame(row.names =colnames(my.lipid.E3vsE417a.data),
                           treat = my.17a)
dataDesignE3conE417a = data.frame(row.names =colnames(my.lipid.E3conE417a.data),
                           E3conE417a = my.E3conE417a)
dataDesignE317aE4con = data.frame(row.names =colnames(my.lipid.E317aE4con.data),
                                  E317aE4con = my.E317aE4con)

# Set null and alternative models
my.model = model.matrix(~ genotype + treatment, data = dataDesign)
my.fit.all <- lmFit(my.lipid.FA.data, design=my.model)
my.fit.all <- eBayes(my.fit.all)

my.model.con = model.matrix(~ con, data = dataDesignCon)
my.fit.all.con <- lmFit(my.lipid.E3vsE4.data, design=my.model.con)
my.fit.all.con <- eBayes(my.fit.all.con)

my.model.E3 = model.matrix(~ E3treat, data = dataDesignE3)
my.fit.all.E3 <- lmFit(my.lipid.E3.data, design=my.model.E3)
my.fit.all.E3 <- eBayes(my.fit.all.E3)

my.model.E4 = model.matrix(~ E4treat, data = dataDesignE4)
my.fit.all.E4 <- lmFit(my.lipid.E4.data, design=my.model.E4)
my.fit.all.E4 <- eBayes(my.fit.all.E4)

my.model.17a = model.matrix(~ treat, data = dataDesign17a)
my.fit.all.17a <- lmFit(my.lipid.E3vsE417a.data, design=my.model.17a)
my.fit.all.17a <- eBayes(my.fit.all.17a)

my.model.E3conE417a = model.matrix(~ E3conE417a, data = dataDesignE3conE417a)
my.fit.all.E3conE417a <- lmFit(my.lipid.E3conE417a.data, design=my.model.E3conE417a)
my.fit.all.E3conE417a <- eBayes(my.fit.all.E3conE417a)

my.model.E317aE4con = model.matrix(~ E317aE4con, data = dataDesignE317aE4con)
my.fit.all.E317aE4con <- lmFit(my.lipid.E317aE4con.data, design=my.model.E317aE4con)
my.fit.all.E317aE4con <- eBayes(my.fit.all.E317aE4con)

my.sig.genotype <- topTable(my.fit.all, coef='genotypeE4WT', p.value = 1, number = Inf)
my.sig.treatment <- topTable(my.fit.all, coef='treatmentCTL', p.value = 1, number = Inf)
my.sig.con <- topTable(my.fit.all.con, p.value = 1 , number = Inf)
my.sig.E3 <- topTable(my.fit.all.E3, p.value = 1 , number = Inf)
my.sig.E4 <- topTable(my.fit.all.E4, p.value = 1 , number = Inf)
my.sig.17a <- topTable(my.fit.all.17a, p.value = 1 , number = Inf)
my.sig.E3conE417a <- topTable(my.fit.all.E3conE417a, p.value = 1, number = Inf)
my.sig.E317aE4con <- topTable(my.fit.all.E317aE4con, p.value = 1, number = Inf)


#check for linear relationship between genotype and treatment
my.merged <- merge(my.sig.E3, my.sig.E4, by = "row.names", suffixes = c(".E3",".E4"))

smoothScatter(my.merged$logFC.E3, my.merged$logFC.E4, xlim = c(-2,2), ylim = c(-2,2))

cor.test(my.merged$logFC.E3, my.merged$logFC.E4, method = "spearman")
abline(0, 1, lty = "dashed", col = "grey")
abline(0,-1, lty = "dashed", col = "grey")

#non-linear relationship, use individual comparisons instead of combined model for analyses
my.fdr <- 0.05


################################################################################
########################  C. Plot Diff Lipids    ###############################
################################################################################
# a. E3 Control vs E4 Control heatmap #70 lipids diff

sig.all.names.con <- rownames(my.sig.con[my.sig.con$adj.P.Val < my.fdr,]) #70 lipids diff

# plot significant heatmap
pdf(paste(my.outprefix,"_Heatmap_E3vsE4_FDR", 100* my.fdr, ".pdf", sep =""), onefile = F)
my.heatmap.title <- paste("EMA E3vsE4 CTL significant (FDR<", 100* my.fdr, "%), ", length(sig.all.names.con), " features",sep="")
pheatmap(my.lipid.data.vsn[sig.all.names.con,], 
         scale = 'row',
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         border_color = NA,
         show_rownames = T,
         main = my.heatmap.title,
         cellwidth = 20
)
dev.off()

# b. E3 Control vs E3 17a heatmap # no lipids diff
sig.all.names.E3 <- rownames(my.sig.E3[my.sig.E3$adj.P.Val < my.fdr,]) 

# c. E4 Control vs E4 17a heatmap # no lipids diff
sig.all.names.E4 <- rownames(my.sig.E4[my.sig.E4$adj.P.Val < my.fdr,])

# d. E3 17a vs E4 17a heatmap # 1 lipid diff

sig.all.names.17a <- rownames(my.sig.17a[my.sig.17a$adj.P.Val < my.fdr,]) 

# e. E3 Control vs E4 17a heatmap #3 lipids diff

sig.all.names.E3conE417a <- rownames(my.sig.E3conE417a[my.sig.E3conE417a$adj.P.Val < my.fdr,]) 

# plot significant heatmap
pdf(paste(my.outprefix,"_Heatmap_E3CTLvsE417a_FDR", 100* my.fdr, ".pdf", sep =""), onefile = F)
my.heatmap.title <- paste("EMA E3CTL vs E417a significant (FDR<", 100* my.fdr, "%), ", length(sig.all.names.E3conE417a), " features",sep="")
pheatmap(my.lipid.data.vsn[sig.all.names.E3conE417a,], 
         scale = 'row',
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         border_color = NA,
         show_rownames = T,
         main = my.heatmap.title,
         cellwidth = 20
)
dev.off()

# f. E3 17a vs E4 con heatmap #1 lipid diff

sig.all.names.E317aE4con <- rownames(my.sig.E317aE4con[my.sig.E317aE4con$adj.P.Val < my.fdr,]) 


# output useful info to file
write.table(my.sig.con[my.sig.con$adj.P.Val < my.fdr,], file = paste(my.outprefix,"__Plasma_Limma_E3vsE4_FDR5", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = T)
write.table(my.sig.con, file = paste(my.outprefix,"_Limma_Plasma_E3vsE4_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = T)

write.table(my.sig.17a[my.sig.17a$adj.P.Val < my.fdr,], file = paste(my.outprefix,"__Plasma_Limma_17a_FDR5", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = T)
write.table(my.sig.17a, file = paste(my.outprefix,"_Limma_Plasma_17a_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = T)

write.table(my.sig.E3conE417a[my.sig.E3conE417a$adj.P.Val < my.fdr,], file = paste(my.outprefix,"__Plasma_Limma_E3ConE417a_FDR5", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = T)
write.table(my.sig.E3conE417a, file = paste(my.outprefix,"_Limma_Plasma_E3ConE417a_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = T)

write.table(my.sig.E317aE4con[my.sig.E317aE4con$adj.P.Val < my.fdr,], file = paste(my.outprefix,"__Plasma_Limma_E317aE4con_FDR5", 100* my.fdr, ".txt", sep =""), quote = F, sep = "\t", row.names = T)
write.table(my.sig.E317aE4con, file = paste(my.outprefix,"_Limma_Plasma_E317aE4con_ALL", ".txt", sep =""), quote = F, sep = "\t", row.names = T)


#Use exported files to input lipids into LION for enrichment analysis

################################################################################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()

