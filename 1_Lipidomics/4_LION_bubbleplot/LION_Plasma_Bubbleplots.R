################################################################################
########################## Make LION Bubbleplots ###############################
################################################################################

# Updated March 7, 2024

setwd("~/Dropbox/2023_EMA_Manuscript/Analyses/20211027Cassie/Exp3")
options(stringsAsFactors = F)
library(ggplot2) 
library(scales) 
theme_set(theme_bw())

# input is data from LION

########################################################
# read and combine
my.APOE4.lion   <- read.csv('LION-enrichment-job2-APOE4.csv')
my.APOE3.lion <- read.csv('LION-enrichment-job3-APOE3.csv')

my.APOE4.lion$Enrichment    <- - my.APOE4.lion$Significant/my.APOE4.lion$Expected
my.APOE3.lion$Enrichment  <- my.APOE3.lion$Significant/my.APOE3.lion$Expected

my.lion.data <- rbind(my.APOE4.lion,my.APOE3.lion)

# filter significant and write to file for table
my.lion.data <- my.lion.data[my.lion.data$FDR.q.value < 0.05,]

my.lion.data$geno <- ifelse(my.lion.data$Enrichment < 0, "APOE4", "APOE3")  # APOE4/APOE3 avg flag

write.table(my.lion.data, file = paste0(Sys.Date(),"_LION_Combined_FDR5_Table.txt"), row.names = F, quote = F, sep = "\t")
########################################################

########################################################
# select top each direction and plot
my.APOE4.lion   <- my.APOE4.lion[my.APOE4.lion$FDR.q.value < 0.05,]
my.APOE3.lion <- my.APOE3.lion[my.APOE3.lion$FDR.q.value < 0.05,]

# remove anything that is significant in both (if any)

my.pos.sort <- sort(my.APOE3.lion$FDR.q.value, index.return = T, decreasing = F) # most significant
my.neg.sort <- sort(my.APOE4.lion$FDR.q.value, index.return = T, decreasing = F)   # most significant

my.lion.data.2 <- rbind(my.APOE3.lion[my.pos.sort$ix[1:5],],
                        my.APOE4.lion[my.neg.sort$ix[1:3],])

# create -log10 FDR for plotting
my.lion.data.2$minlog10fdr  <- -log10(my.lion.data.2$FDR.q.value)
colnames(my.lion.data.2)[2] <- "Description"

# create and preserve wanted display order
my.lion.data.2$geno <- ifelse(my.lion.data.2$Enrichment < 0, "APOE4", "APOE3")  # APOE4/APOE3 avg flag
my.lion.geno.sorted <- my.lion.data.2[order(my.lion.data.2$geno), ]

my.lion.geno.sorted$Tissue <- factor(rep("Plasma",length(  my.lion.geno.sorted$Description)))
my.lion.geno.sorted$Description <- factor(my.lion.geno.sorted$Description, levels = rev(unique(my.lion.geno.sorted$Description)))

# APOE3/APOE4 color scale
my.max <- max(my.lion.geno.sorted$Enrichment)
my.min <- min(my.lion.geno.sorted$Enrichment)

# make bubble plot
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
#my.color.vector.geno <- c("steelblue4","skyblue3","lightskyblue1","lightcyan","white","rosybrown1","indianred1","red","red4")
my.color.vector.geno <- c("red4","red","indianred1","rosybrown1","white","lightcyan","lightskyblue1","skyblue3","steelblue4")
my.plot <- ggplot(my.lion.geno.sorted,aes(x=Tissue,y=Description,colour=Enrichment,size=minlog10fdr))+ theme_bw()+ geom_point(shape = 16)
my.plot <- my.plot + ggtitle("LION Analysis") + labs(x = "-log10(pvalue)", y = "")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector.geno, na.value = "grey50", guide = "colourbar", values = my.scaled)
my.plot

my.pdfname <- paste(Sys.Date(),"LION_geno_BALLOON_plot_significant_pathways.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 6, width=6)
print(my.plot)
dev.off()  

