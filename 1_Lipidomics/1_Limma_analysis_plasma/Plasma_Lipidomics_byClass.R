################################################################################
#################### EMA Plasma Lipidomics Class ###############################
################################################################################

# Updated March 7, 2024

#set working directory
setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/20211027Cassie/Exp3")
options(stringsAsFactors = F)
library(beeswarm)
library(limma)
library(pheatmap)

# add a global data-driven normalization, which is supposed to help with robustness in metabolomics
# VSN : 
#    - https://www.nature.com/articles/srep38881
#    - https://academic.oup.com/bioinformatics/article/30/15/2155/2390365

# use lm to determine effect and significance by class

################################################################################
########################   Class analysis     ##################################
################################################################################

# 1. read raw data and set column order
my.lipid.FA.data.1 <- read.csv('2023-06-16_Lipidomics_Plasma_Class_cleaned.txt', header = T, sep = "\t")
my.E3.ctl    <- grep("E3WT_ctl", colnames(my.lipid.FA.data.1))
my.E3.17a    <- grep("E3WT_17a", colnames(my.lipid.FA.data.1))
my.E4.ctl    <- grep("E4WT_ctl",colnames(my.lipid.FA.data.1))
my.E4.17a    <- grep("E4WT_17a",colnames(my.lipid.FA.data.1))

# 2. subset and reorder
my.lipid.FA.data <- my.lipid.FA.data.1[,c(my.E3.ctl,my.E3.17a,my.E4.ctl,my.E4.17a)]

################################################################################
# 3. perform variance stabilizing normalization
my.lipid.data.vsn <- as.data.frame(normalizeVSN(as.matrix(my.lipid.FA.data)))
rownames(my.lipid.data.vsn) <- rownames(my.lipid.FA.data)
colnames(my.lipid.data.vsn) <- colnames(my.lipid.FA.data)

################################################################################
# 4. retrieve annotations
my.lipid.data.vsn <- cbind(rownames(my.lipid.data.vsn),my.lipid.data.vsn)
colnames(my.lipid.data.vsn)[1] <- "Lipid_ID"



library(dplyr)
library(tidyr)

#warning message ok here, it's separating the row names into two columns
my.lipid.data.vsn <- my.lipid.data.vsn %>% separate(Lipid_ID, c("Lipid_Family", "Lipid_ID"))

################################################################################
# 5. sum by lipid type
my.lipid.data.type <- aggregate(my.lipid.data.vsn[,-c(1:2)], 
                                by = list(my.lipid.data.vsn$Lipid_Family), 
                                FUN = sum)
rownames(my.lipid.data.type) <- my.lipid.data.type$Group.1

# 6. export for downstream analysis
write.table(my.lipid.data.type,file = paste0(Sys.Date(), "_Plasma_Lipid_by_Class.txt"), sep = "\t", quote = F)

#heatmap of different classes
pdf(paste(Sys.Date(),"Plasma_Lipidomics_by_Lipid_Family_BCA_VSN_heatmap_SORTED_BY_CLASS.pdf",sep = "_"), width = 8, height = 5, onefile = F)
pheatmap(my.lipid.data.type[c("CE","Cer","DG","FA","HexCER","LPC","LPE","PA","PC","PG","PI","SM","TG"),-1], 
         scale = 'row', 
         cluster_cols = F,
         cluster_rows = T,
         colorRampPalette(rev(c("#CC3333","#FF9999","#FFCCCC","white","#CCCCFF","#9999FF","#333399")))(50),
         cellwidth = 10)
dev.off()


#######################
sink(file = paste(my.outprefix,"_session_Info.txt", sep =""))
sessionInfo()
sink()
