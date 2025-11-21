#Loading packages
library(tidyverse)
library(dplyr)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(Rcpp)
library(EnhancedVolcano)
library(statmod)
library(limma)
library(gplots)
library(RColorBrewer)
library(janitor)

#Reading data into R
# List of barcodes
samples <- c("sample1", "sample2", "sample3", "sample4", "sample5", "sample6",
             "sample7", "sample8", "sample9", "sample10", "sample11", "sample12",
             "sample13", "sample14", "sample15", "sample16", "sample17", "sample18",
             "sample19", "sample20", "sample21", "sample22", "sample23", "sample24",
             "sample25", "sample26", "sample27", "sample28", "sample29", "sample30")

samples <- c("sample31", "sample32", "sample33", "sample34", "sample35", "sample36",
             "sample37", "sample38", "sample39", "sample40", "sample41", "sample42",
             "sample43", "sample44", "sample45", "sample46", "sample47", "sample48",
             "sample49", "sample50", "sample51", "sample52", "sample53", "sample54",
             "sample55", "sample56", "sample57", "sample58", "sample59", "sample60")

samples <- c("sample61", "sample62", "sample63", "sample64", "sample65", "sample66",
             "sample67", "sample68", "sample69", "sample70", "sample71", "sample72",
             "sample73", "sample74", "sample75", "sample76", "sample77", "sample78",
             "sample79", "sample80")

# Initialize an empty list to store data frames
sample_data <- list()

# Loop through each barcode and read the corresponding file
for (sample in samples) {
  file_name <- paste0(sample, "_counts.tsv")  # Construct the file name
  sample_data[[sample]] <- read.csv(file_name, sep = '\t', header = TRUE)  # Read and store in list
}

##Accessing data from each specific barcode
for (i in 1:length(samples)) {
  assign(paste0("sample", i), sample_data[[samples[i]]])
}

#creating a variable with matching file names
samples <- paste0("sample", 1:30)

##picking only the necessary collumns
for (i in 1:length(samples)) {
  current_sample <- get(samples[i])
  current_sample <- current_sample[, c(1, 2)]
  assign(samples[i], current_sample)
}

#merging the different samples into one dataframe
all_transcripts <- Reduce(function(x,y) merge(x,y,by="Gene_id",all=TRUE) ,list(sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12,
                                                                               sample13, sample14, sample15, sample16, sample17, sample18, sample19, sample20, sample21, sample22, sample23, sample24,
                                                                               sample25, sample26, sample27, sample28, sample29, sample30))

all_transcripts <- Reduce(function(x,y) merge(x,y,by="Gene_id",all=TRUE) ,list(sample1, sample2, sample3, sample4, sample5, sample6, sample7, sample8, sample9, sample10, sample11, sample12,
                                                                               sample13, sample14, sample15, sample16, sample17, sample18, sample19, sample20))
#Replacing NAs with 0
all_transcripts[is.na(all_transcripts)] <- 0

#filtering riboosomal RNA
#filtered_data <- all_transcripts[!all_transcripts$Gene_id %in% ribosomal$Gene_id, ]

#Summing to identify the most abundant transcript
cumulation <- all_transcripts %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

#Reading the sample metadata into R
sample_metadata <- read.csv("sample_metadata.csv", sep = ',', header = TRUE)

#converting the transcript name collumn into row names
read_count <- all_transcripts %>% remove_rownames %>% column_to_rownames(var="Gene_id")

#making a DGE list
dgeObj <- DGEList(read_count, group = sample_metadata$MB_infection_status)


# Filter lowly expressed genes
keep <- filterByExpr(dgeObj)
dgeObj <- dgeObj[keep, , keep.lib.sizes = FALSE]

# Step 4: Normalize for library size using TMM
dgeObj <- calcNormFactors(dgeObj, method = "TMM")

# Check normalization factors
dgeObj$samples

#Normalize the data
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

# Log-transform CPM (adding a small offset to avoid log of zero)
logCPM <- cpm(dgeObj, log = TRUE, prior.count = 1)
pca <- prcomp(t(logCPM), scale. = TRUE)

summary(pca)

#Extract the PCA component
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = sample_metadata$MB_infection_status)

#plot the PCA 
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "24 hours PBM guts", x = "PC1", y = "PC2") +
  theme_minimal()

#####For differential expression expression analysis#######
# Apply sample grouping based on sample_type from which the sample was derived
design <- model.matrix(~0+sample_metadata$MB_infection_status)
colnames(design) <- levels(as.factor(sample_metadata$MB_infection_status))

# Estimate dispersions for tags
filtered.counts.dge <- estimateDisp(dgeObj, design, robust = TRUE)

# Fit a generalized likelihood model to the DGELIST using sample grouping
fit <- glmFit(filtered.counts.dge,design)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(as.factor(sample_metadata$MB_infection_status)), 2))

comparisons <- list()
for (i in 1:nrow(condition_pairs)) {
  comparisons[[i]] <- as.character(condition_pairs[i,])
}

# vector to store deferentially expressed genes
sig_genes <- c()

# iterate over the contrasts, and perform a differential expression test for each pair
for (conds in comparisons) {
  # generate string contrast formula
  contrast_formula <- paste(conds, collapse=' - ') 
  
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=design)
  contrast_lrt <- glmLRT(fit, contrast=contrast_mat)
  topGenes <- topTags(contrast_lrt, n=Inf, p.value=0.05, adjust.method = "BH")
  
  # Grab highly ranked genes
  sig_genes <- union(sig_genes, rownames(topGenes$table))
}
sig_genes

# Filter out genes which were not differentially expressed for any contrast
de.genes <- filtered.counts.dge[rownames(filtered.counts.dge) %in% sig_genes,]

dim(de.genes$counts)
DE_data <- de.genes$counts

# Mbpos compared to Mbneg
MBNEG_vs_MBNEG_PLplus_lrt <- glmLRT(fit, contrast=c(-1,1,0,0))
MBNEG_vs_MbposPLneg_lrt <- glmLRT(fit, contrast=c(-1,0,1,0))
MBNEG_PLplus_vs_Mbpos_PLplus_lrt <- glmLRT(fit, contrast=c(0,-1,0,1))
Mbpos_PLneg_vs_Mbpos_PLplus_lrt <- glmLRT(fit, contrast=c(0,0,-1,1))

#MBNEG compared to MBNEGPLplus
topGenes_MBNEG_vs_MBNEG_PLplus<- topTags(MBNEG_vs_MBNEG_PLplus_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_MBNEG_vs_MBNEG_PLplus)
DE_topGenes_MBNEG_vs_MBNEG_PLplus <- topGenes_MBNEG_vs_MBNEG_PLplus$table
write.csv(DE_topGenes_MBNEG_vs_MBNEG_PLplus, "External_vs_ExternalPLplus.csv")

#MBNEG compared to MbposPLneg
topGenes_MBNEG_vs_MbposPLneg <- topTags(MBNEG_vs_MbposPLneg_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_MBNEG_vs_MbposPLneg)
DE_topGenes_MBNEG_vs_MbposPLneg <- topGenes_MBNEG_vs_MbposPLneg$table
write.csv(DE_topGenes_MBNEG_vs_MbposPLneg, "External_vs_MbposPLneg.csv")

#MBNEG_PLplus compared to Mbpos_PLplus
topGenes_MBNEG_PLplus_vs_Mbpos_PLplus <- topTags(MBNEG_PLplus_vs_Mbpos_PLplus_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_MBNEG_PLplus_vs_Mbpos_PLplus)
DE_topGenes_MBNEG_PLplus_vs_Mbpos_PLplus <- topGenes_MBNEG_PLplus_vs_Mbpos_PLplus$table
write.csv(DE_topGenes_MBNEG_PLplus_vs_Mbpos_PLplus, "ExternalPLplus_vs_MbposPLplus.csv")

#MBpos_PLneg compared to Mbpos_PLplus
topGenes_Mbpos_PLneg_vs_Mbpos_PLplus <- topTags(Mbpos_PLneg_vs_Mbpos_PLplus_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_Mbpos_PLneg_vs_Mbpos_PLplus)
DE_topGenes_Mbpos_PLneg_vs_Mbpos_PLplus <- topGenes_Mbpos_PLneg_vs_Mbpos_PLplus$table
write.csv(DE_topGenes_Mbpos_PLneg_vs_Mbpos_PLplus, "MbposPLneg_vs_MbposPLplus.csv")

#Gene descriptions
description <- read.csv("description.csv", sep = ',', header = TRUE)
#converting the rowname to collumn name
DE_topGenes_Mbpos_PLneg_vs_Mbpos_PLplus <- rownames_to_column(DE_topGenes_Mbpos_PLneg_vs_Mbpos_PLplus, "Input_ID")
#merging the dataset to the gene description
dataset <- merge(DE_topGenes_Mbpos_PLneg_vs_Mbpos_PLplus, description, by = "Input_ID", all = FALSE)

