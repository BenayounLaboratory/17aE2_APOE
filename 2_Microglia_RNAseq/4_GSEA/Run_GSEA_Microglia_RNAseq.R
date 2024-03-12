################################################################################
#################### EMA Microglia RNA-seq GSEA  ###############################
################################################################################

setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/June2023")
options(stringsAsFactors = FALSE)

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
require(DOSE)

# set the desired organism

organism = library('org.Mm.eg.db', character.only = TRUE)

my.control  <- read.table('2023-06-19_EMA_microglia_DESeq2_Analysis_APOE_Genotype_ControlsOnly_all_genes_statistics.txt', sep = "\t", header = T)
my.treat    <- read.table('2023-06-19_EMA_microglia_DESeq2_Analysis_TREATMENT_only_all_genes_statistics.txt', sep = "\t", header = T)
my.E3       <- read.table('2023-06-19_EMA_microglia_DESeq2_Analysis_APOE3_all_genes_statistics.txt', sep = "\t", header = T)
my.E4       <- read.table('2023-06-19_EMA_microglia_DESeq2_Analysis_APOE4_all_genes_statistics.txt', sep = "\t", header = T)
my.E4treatE3con <- read.table('2023-06-19_EMA_microglia_DESeq2_Analysis_APOE4treatvsAPOE3con_all_genes_statistics.txt', sep = "\t", header = T)
my.E3treatE4con <- read.table('2023-06-19_EMA_microglia_DESeq2_Analysis_APOE3treatvsAPOE4con_all_genes_statistics.txt', sep = "\t", header = T)


# genes lists for enrichment analysis
genes.control   <- rownames(my.control)   
genes.treat      <- rownames(my.treat) 
genes.E3 <- rownames(my.E3)
genes.E4 <- rownames(my.E4)
genes.E4treatE3con <- rownames(my.E4treatE3con)
genes.E3treatE4con <- rownames(my.E3treatE4con)


# keytypes(org.Mm.eg.db)
entrezID.genes.control     <- bitr(genes.control   , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.genes.treat        <- bitr(genes.treat      , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.genes.E3     <- bitr(genes.E3   , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.genes.E4        <- bitr(genes.E4      , fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
entrezID.genes.E4treatE3con <- bitr(genes.E4treatE3con, fromType = "SYMBOL", toType = "ENTREZID", OrgDb ="org.Mm.eg.db")
entrezID.genes.E3treatE4con <- bitr(genes.E3treatE4con, fromType = "SYMBOL", toType = "ENTREZID", OrgDb ="org.Mm.eg.db")

# we want the log2 fold change 
original_gene_list.control <- my.control$log2FoldChange
head(original_gene_list.control)

original_gene_list.treat <- my.treat$log2FoldChange
head(original_gene_list.treat)

original_gene_list.E3 <- my.E3$log2FoldChange
head(original_gene_list.E3)

original_gene_list.E4 <- my.E4$log2FoldChange
head(original_gene_list.E4)

original_gene_list.E4treatE3con <- my.E4treatE3con$log2FoldChange
head(original_gene_list.E4treatE3con)

original_gene_list.E3treatE4con <- my.E3treatE4con$log2FoldChange
head(original_gene_list.E3treatE4con)


# name the vector
names(original_gene_list.control) <- rownames(my.control)
names(original_gene_list.treat) <- rownames(my.treat)
names(original_gene_list.E3) <- rownames(my.E3)
names(original_gene_list.E4) <- rownames(my.E4)
names(original_gene_list.E4treatE3con) <- rownames(my.E4treatE3con)
names(original_gene_list.E3treatE4con) <- rownames(my.E3treatE4con)

# omit any NA values 
gene_list.control<-na.omit(original_gene_list.control)
gene_list.treat<-na.omit(original_gene_list.treat)
gene_list.E3<-na.omit(original_gene_list.E3)
gene_list.E4<-na.omit(original_gene_list.E4)
gene_list.E4treatE3con <- na.omit(original_gene_list.E4treatE3con)
gene_list.E3treatE4con <- na.omit(original_gene_list.E3treatE4con)

# sort the list in decreasing order (required for clusterProfiler)
gene_list.control = sort(gene_list.control, decreasing = TRUE)
head(gene_list.control)
gene_list.treat = sort(gene_list.treat, decreasing = TRUE)
head(gene_list.treat)
gene_list.E3 = sort(gene_list.E3, decreasing = TRUE)
head(gene_list.E3)
gene_list.E4 = sort(gene_list.E4, decreasing = TRUE)
head(gene_list.E4)
gene_list.E4treatE3con = sort(gene_list.E4treatE3con, decreasing = TRUE)
head(gene_list.E4treatE3con)
gene_list.E3treatE4con = sort(gene_list.E3treatE4con, decreasing = TRUE)
head(gene_list.E3treatE4con)

################################################################################
#################### EMA Microglia RNA-seq GSEA ALL  ############################
################################################################################

# control groups
gse.control <- gseGO(geneList=gene_list.control, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "none")

gsea.res.control <- as.data.frame(gse.control)

write.table(gsea.res.control, file = paste(Sys.Date(),"Control_GSEA_ALL_FDR5.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Dot plot
my.control.out <- paste("Control_GSEA_ALL.pdf",sep="_")
pdf(my.control.out, width = 10, height = 10)
dotplot(gse.control, showCategory=20, split=".sign") +facet_grid(.~.sign)
dev.off()

my.control.out <- paste("Control_GSEA_ALL_Top20.pdf",sep="_")
pdf(my.control.out, width = 10, height = 10)
dotplot(gse.control, showCategory=20)
dev.off()

# treated groups
gse.treat <- gseGO(geneList=gene_list.treat, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = org.Mm.eg.db, 
                     pAdjustMethod = "none")


gsea.res.treat <- as.data.frame(gse.treat)

write.table(gsea.res.treat, file = paste(Sys.Date(),"Treat_GSEA_ALL_FDR5.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Dot plot
my.treat.out <- paste("Treat_GSEA_ALL.pdf",sep="_")
pdf(my.treat.out, width = 10, height = 10)
dotplot(gse.treat, showCategory=20, split=".sign") +facet_grid(.~.sign)
dev.off()

my.treat.out <- paste("Treat_GSEA_ALL_Top20.pdf",sep="_")
pdf(my.treat.out, width = 10, height = 10)
dotplot(gse.treat, showCategory=20)
dev.off()

# APOE3
gse.E3 <- gseGO(geneList=gene_list.E3, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = org.Mm.eg.db, 
                     pAdjustMethod = "none")

gsea.res.E3 <- as.data.frame(gse.E3)

write.table(gsea.res.E3, file = paste(Sys.Date(),"APOE3_GSEA_ALL_FDR5.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Dot plot
my.E3.out <- paste("APOE3_GSEA_ALL.pdf",sep="_")
pdf(my.E3.out, width = 10, height = 10)
dotplot(gse.E3, showCategory=20, split=".sign") +facet_grid(.~.sign)
dev.off()

my.E3.out <- paste("APOE3_GSEA_ALL_Top20.pdf",sep="_")
pdf(my.E3.out, width = 10, height = 10)
dotplot(gse.E3, showCategory=20)
dev.off()

# APOE4

gse.E4 <- gseGO(geneList=gene_list.E4, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = org.Mm.eg.db, 
                     pAdjustMethod = "none")

gsea.res.E4 <- as.data.frame(gse.E4)

write.table(gsea.res.E4, file = paste(Sys.Date(),"APOE4_GSEA_ALL_FDR5.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Dot plot
my.E4.out <- paste("APOE4_GSEA_ALL.pdf",sep="_")
pdf(my.E4.out, width = 10, height = 10)
dotplot(gse.E4, showCategory=20, split=".sign") +facet_grid(.~.sign)
dev.off()

my.E4.out <- paste("APOE4_GSEA_ALL_Top20.pdf",sep="_")
pdf(my.E4.out, width = 10, height = 10)
dotplot(gse.E4, showCategory=20)
dev.off()


# APOE4 treated vs APOE3 control

gse.E4treatE3con <- gseGO(geneList=gene_list.E4treatE3con, 
                ont ="ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Mm.eg.db, 
                pAdjustMethod = "none")

gsea.res.E4treatE3con <- as.data.frame(gse.E4treatE3con)

write.table(gsea.res.E4treatE3con, file = paste(Sys.Date(),"APOE4TreatAPOE3Con_GSEA_ALL_FDR5.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Dot plot
my.E4treatE3con.out <- paste("APOE4TreatAPOE3Con_GSEA_ALL.pdf",sep="_")
pdf(my.E4treatE3con.out, width = 10, height = 10)
dotplot(gse.E4treatE3con, showCategory=20, split=".sign") +facet_grid(.~.sign)
dev.off()

my.E4treatE3con.out <- paste("APOE4TreatAPOE3Con_GSEA_ALL_Top20.pdf",sep="_")
pdf(my.E4.out, width = 10, height = 10)
dotplot(gse.E4treatE3con, showCategory=20)
dev.off()

# APOE3 treated vs APOE4 control

gse.E3treatE4con <- gseGO(geneList=gene_list.E3treatE4con, 
                ont ="ALL", 
                keyType = "SYMBOL", 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = org.Mm.eg.db, 
                pAdjustMethod = "none")

gsea.res.E3treatE4con <- as.data.frame(gse.E3treatE4con)

write.table(gsea.res.E3treatE4con, file = paste(Sys.Date(),"APOE3treatAPOE4con_GSEA_ALL_FDR5.txt", sep ="_"), sep = "\t" , row.names = T, quote=F)

# Dot plot
my.E3treatE4con.out <- paste("APOE3TreatAPOE4Con_GSEA_ALL.pdf",sep="_")
pdf(my.E3treatE4con.out, width = 10, height = 10)
dotplot(gse.E3treatE4con, showCategory=20, split=".sign") +facet_grid(.~.sign)
dev.off()

my.E3treatE4con.out <- paste("APOE3TreatAPOE4con_GSEA_ALL_Top20.pdf",sep="_")
pdf(my.E3treatE4con.out, width = 10, height = 10)
dotplot(gse.E3treatE4con, showCategory=20)
dev.off()


#######################
sink(paste0(Sys.Date(),"_Enrichment_analysis_sessionInfo.txt"))
sessionInfo()
sink()
