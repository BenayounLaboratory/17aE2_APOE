################################################################################
#################### EMA Microglia GSEA Bubbleplot #############################
################################################################################

#March 12, 2024

# load libraries

library('pheatmap') 
library(ggplot2) 
library(scales) 

# set working directory
setwd("C:/Users/livis/Dropbox/2023_EMA_Manuscript_Cassie/Code/Input Files/RNA-seq")

# set options
options(stringsAsFactors = F)

####################
##### FUNCTION #####
####################
# This function makes bubble plots from APOE3 and APOE4 GSEA statistics (FDR 5)
# INPUT: my.data.name = output file name describing the gene set name and FDR
#        FDR = false discovery rate (0.05)
#        my.threshold = if threshold is 2, the bubble plot will display pathways
#                 that are present in at least 2 or more columns
#        my.colnames = names of columns for the bubble plot

get_enrich_balloons <- function(my.data.name, FDR, my.threshold = 2,
                                my.colnames = c("APOE3", "APOE4") ) {
  
  if (FDR == 5) {
    my.enrich.sets <- list.files(".", pattern = "*FDR5.txt")
  }
  
  # reorder files based on colnames
  my.columns <<- c() 
  
  for (i in 1:length(my.colnames)) {
    my.columns <<- c(my.columns,
                     grep(paste("_", my.colnames[i], "_", sep = ""), my.enrich.sets))
  }
  
  # get data from significant FDR
  my.tissues <- vector(length=length(my.colnames), mode="list")
  names(my.tissues) <- my.colnames
  my.pathways <- c()
  
  for ( i in 1:length(my.columns)) {
    my.file <- my.enrich.sets[my.columns[i]]
    print(my.file);
    my.tissues[[i]] <- read.csv(paste("./", my.file, sep=""), sep="\t", header=T) 
    row.names(my.tissues[[i]]) <- my.tissues[[i]]$Description
    my.pathways <- unique(c(my.pathways, rownames(my.tissues[[i]])))
    
  }
  print(my.pathways)
  
  ## prepapre output data
  # p-val matrix
  my.matrix <<- matrix(0,length(my.pathways),length(my.colnames)) # default: -log10(1) pval == 0 no enrichment
  
  # Enrichment matrix
  my.matrix2 <<- matrix(0,length(my.pathways),length(my.colnames)) # initialize with Enrichment = 0 if no enrich
  
  # matrix with record of significance
  my.matrix3 <<- matrix(0,length(my.pathways),length(my.colnames)) # to get sigificant pathways
  
  colnames(my.matrix)  <<- my.colnames
  colnames(my.matrix2) <<- my.colnames
  rownames(my.matrix)  <<- my.pathways
  rownames(my.matrix2) <<- my.pathways
  colnames(my.matrix3) <<- my.colnames
  rownames(my.matrix3) <<- my.pathways
  
  # collect data from files
  for (i in 1:length(my.pathways)) {
    #print(my.pathways[i])
    
    for (j in 1:length(my.colnames)) { # tissues 
      
      my.id <- which(rownames(my.tissues[[j]]) %in% my.pathways[i])
      #print(paste(j, length(my.id)))
      if(length(my.id) == 1) { # if was significant in this geno (and not on both tail ends, which would be 2)
      #print('IN')
        
        my.matrix[i,j] <<- -log10(my.tissues[[j]]$p.adjust[my.id]+1e-10) # log(0) is undefined
        
        my.matrix2[i,j] <<- my.tissues[[j]]$NES[my.id]
        
        my.matrix3[i,j] <<- 1
        
      }
    }
  }
  
  # find pathways significant in threshold (minimum) or more tissues (2 potential from our experiments)
  my.sigs <- apply(my.matrix3[,c(1:2)],1,sum) >= my.threshold    
  my.res.enrich <<- data.frame(my.matrix2[my.sigs,])
  my.pval.enrich <- data.frame(my.matrix[my.sigs,])
  
  # sort by average change
  my.average <- apply(my.res.enrich[,c(1,2)],1,mean)
  my.sorted <- sort(my.average,index.return=T,decreasing=T)
  my.res.enrich2 <<- my.res.enrich[my.sorted$ix,]
  
  my.pval.enrich2 <- data.frame(my.pval.enrich[my.sorted$ix,])
  
  my.txtname <- paste(Sys.Date(),"Microglia_APOE3_APOE4_Enrichment_table",my.data.name,"pathways_significant_in",my.threshold,"or_more_pvalenrich.txt", sep="_")
  
  write.table(my.res.enrich2,file=my.txtname,sep="\t",quote=F)
  
  my.res.enrich2$Pathnames <<- rownames(my.res.enrich2)
  
  # format for ggplot
  my.res.enrich2$Pathnames <<- rownames(my.res.enrich2)
  my.res.enrich3 <<- cbind(my.res.enrich2[,c('Pathnames',my.colnames[1])],rep(my.colnames[1],dim(my.res.enrich2)[1]),my.pval.enrich2[,my.colnames[1]])
  colnames(my.res.enrich3) <<- c('Pathnames','NES','condition','minusLog10Pval')
  for ( h in 2:length(my.colnames)) {
    my.new <- cbind(my.res.enrich2[,c('Pathnames',my.colnames[h])],rep(my.colnames[h],dim(my.res.enrich2)[1]),my.pval.enrich2[,my.colnames[h]])
    colnames(my.new) <- colnames(my.res.enrich3)
    my.res.enrich3 <<- rbind(my.res.enrich3, 
                             my.new)
  }
  
  my.max <<- max(my.res.enrich3$NES)
  my.min <<- min(my.res.enrich3$NES)
  my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
  my.scaled <- rescale(my.values, to = c(0, 1))
  my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
  
  # to preserve the wanted order
  my.res.enrich3$condition <- factor(my.res.enrich3$condition, levels = unique(my.res.enrich3$condition))
  my.res.enrich3$Pathnames <- factor(my.res.enrich3$Pathnames, levels = rev(unique(my.res.enrich3$Pathnames)))
  
  my.pdfname <- paste(Sys.Date(),"Microglia_APOE3_APOE4_Enrichment_BALLOON_plot",my.data.name,"pathways_significant_in",my.threshold,"or_more.pdf", sep="_")
  
  pdf(my.pdfname, onefile=T, height = max(5, sum(my.sigs)/3), width=24)
  
  my.plot <- ggplot(my.res.enrich3,aes(x = condition, y = reorder(minusLog10Pval, order) ,colour = NES, size = minusLog10Pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16) + theme_bw()
  my.plot <- ggplot(my.res.enrich3,aes(x=condition,y=Pathnames,colour=NES,size=minusLog10Pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16)
  
  my.plot <- my.plot + ggtitle("GSEA") + labs(x = "Condition", y = "Pathways")
  my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
  print(my.plot)
  dev.off()  
}


# source code
get_enrich_balloons("GSEA_ALL_FDR5_pathways", 5)

#####To get the most significant pathways####
# Order the data by 'minusLog10Pval' in descending order
#ordered_data <- my.res.enrich3[order(my.res.enrich3$minusLog10Pval, decreasing = TRUE), ]

# Select the top 30 rows
#top_30_rows <- head(ordered_data, 30)

# Remove duplicates from the "Pathnames" column
#unique_pathnames <- unique(top_30_rows$Pathnames)

# Order the unique data by 'minusLog10Pval' again
#ordered_unique_data <- top_30_rows[order(top_30_rows$minusLog10Pval, decreasing = TRUE), ]

# Select the top 20 rows
#top_20_final <- head(ordered_unique_data, 20)

# Extract the "Pathnames" column
#selected_pathnames <- top_20_final$Pathnames

# Print or use the selected_pathnames vector as needed
#print(selected_pathnames)

# used my.res.enrich3 and pulled the top and bottom 5 NES for both genotypes

# Create a function to get the indices of the 5 highest and 5 lowest values within each group
get_indices <- function(x, condition) {
  indices <- c(head(order(x), 5), tail(order(x), 5))
  condition_labels <- rep(condition, length(indices))
  return(data.frame(indices, condition_labels))
}

# Initialize an empty data frame to store indices and condition labels
all_indices <- data.frame(indices = integer(0), condition_labels = character(0), stringsAsFactors = FALSE)

# Apply the function within each group for each condition
unique_conditions <- unique(my.res.enrich3$condition)
for (condition_value in unique_conditions) {
  indices_df <- get_indices(my.res.enrich3$NES[my.res.enrich3$condition == condition_value], condition_value)
  all_indices <- rbind(all_indices, indices_df)
}

# Select the corresponding rows
filtered_data <- my.res.enrich3[all_indices$indices, ]
filtered_data$condition <- all_indices$condition_labels


# for specific pathways of interest
my.selected.pathways<- filtered_data$Pathnames
my.selected.pathways <- unique(my.selected.pathways)
my.selected.pathways <- c(my.selected.pathways, paste0(my.selected.pathways, "1"))
my.res.enrich4 <- my.res.enrich3[my.selected.pathways,]
my.res.enrich4 <- na.omit(my.res.enrich4)
my.res.enrich4$order <- c(38:1)

my.pdfname <- paste(Sys.Date(),"Microglia_APOE3_APOE4_Enrichment_BALLOON_plot_TOP5updownNES_eachgeno.pdf", sep="_")

pdf(my.pdfname, onefile=T, height = 6, width=12)
my.plot <- ggplot(my.res.enrich4, aes(x = condition, y = reorder(Pathnames, order) ,colour = NES, size = minusLog10Pval))+ theme(text = element_text(size=16))+ geom_point(shape = 16) + theme_bw()
my.color.vector <- c("darkblue","dodgerblue4","dodgerblue3","dodgerblue1","white","lightcoral","brown1","firebrick2","firebrick4")
my.max <<- max(my.res.enrich4$NES)
my.min <<- min(my.res.enrich4$NES)
my.values <- c(my.min,0.75*my.min,0.5*my.min,0.25*my.min,0,0.25*my.max,0.5*my.max,0.75*my.max,my.max)
my.scaled <- rescale(my.values, to = c(0, 1))
my.plot <- my.plot + ggtitle("GSEA") + labs(x = "Condition", y = "Pathways")
my.plot <- my.plot + scale_colour_gradientn(colours = my.color.vector,space = "Lab", na.value = "grey50", guide = "colourbar", values = my.scaled)
print(my.plot)
dev.off()  

#######################
sink(paste0(Sys.Date(),"_Bubbleplot_sessionInfo.txt"))
sessionInfo()
sink()
