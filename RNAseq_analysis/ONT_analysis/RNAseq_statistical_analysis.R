##RNAseq statistical analysis script inclusive of PCA plotting, differential expression analysis and visualization, enrichment analysis

getwd()
setwd("/Users/jwahura/Desktop/RNAseq_analysis/sugar_fed_guts/all_corrected_reads/")

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
library(Glimma)
library(org.Ag.eg.db)
library(gplots)
library(RColorBrewer)

#Reading data into R
# List of barcodes
barcodes <- c("barcode01", "barcode02", "barcode03", "barcode04", "barcode05", "barcode06",
              "barcode07", "barcode08", "barcode09", "barcode10", "barcode11", "barcode12")

# Initialize an empty list to store data frames
barcode_data <- list()

# Loop through each barcode and read the corresponding file
for (barcode in barcodes) {
  file_name <- paste0(barcode, "_counts.tsv")  # Construct the file name
  barcode_data[[barcode]] <- read.csv(file_name, sep = '\t', header = TRUE)  # Read and store in list
}

##Accessing data from each specific barcode
for (i in 1:length(barcodes)) {
  assign(paste0("barcode", i), barcode_data[[barcodes[i]]])
}

#creating a variable with matching file names
barcodes <- paste0("barcode", 1:12)

##picking only the necessary collumns
for (i in 1:length(barcodes)) {
  current_barcode <- get(barcodes[i])
  current_barcode <- current_barcode[, c(1, 3)]
  assign(barcodes[i], current_barcode)
}

#merging the different samples into one dataframe
all_transcripts <- Reduce(function(x,y) merge(x,y,by="transcript_name",all=TRUE) ,list(barcode1, barcode2, barcode3, barcode4, barcode5, barcode6, barcode7, barcode8, barcode9, barcode10, barcode11, barcode12))

#Replacing NAs with 0
all_transcripts[is.na(all_transcripts)] <- 0

#Summing to identify the most abundant transcript
cumulation <- all_transcripts %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

#Reading the sample metadata into R
sample_metadata <- read.csv("sample_metadata.csv", sep = ',', header = TRUE)

#converting the transcript name collumn into row names
read_count <- all_transcripts %>% remove_rownames %>% column_to_rownames(var="transcript_name")

#Obtain CPMs
myCPM <- cpm(read_count)
#Have a look at the output
head(myCPM)

##Filtering lowly abundant transcripts/genes
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

# Summary of how many TRUEs there are in each row
table(rowSums(thresh))

# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- read_count[keep,]
summary(keep)
dim(counts.keep)

#making a DGE list
dgeObj <- DGEList(counts.keep, group = sample_metadata$sample_type)

#listing the dge object
dgeObj
# See what slots are stored in dgeObj
names(dgeObj)
# Library size information is stored in the samples slot
dgeObj$samples

#Normalize the data
dgeObj <- calcNormFactors(dgeObj)
dgeObj$samples

# Log-transform CPM (adding a small offset to avoid log of zero)
logCPM <- cpm(dgeObj, log = TRUE, prior.count = 1)

#Run PCA
pca <- prcomp(t(logCPM))

#Extract the PCA component
pca_data <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], Group = sample_metadata$sample_type)

#plot the PCA 
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4) +
  labs(title = "BF72_fat body", x = "PC1", y = "PC2") +
  theme_minimal()

#####For differential expression expression analysis#######
# Apply sample grouping based on sample_type from which the sample was derived
design <- model.matrix(~0+sample_metadata$sample_type)
colnames(design) <- levels(as.factor(sample_metadata$sample_type))

# Estimate dispersions for tags
filtered.counts.dge <- estimateDisp(dgeObj, design, robust = TRUE)

# Fit a generalized likelihood model to the DGELIST using sample grouping
fit <- glmFit(filtered.counts.dge,design)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(as.factor(sample_metadata$sample_type)), 2))

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
MB_NEG_vs_Mbneg_lrt <- glmLRT(fit, contrast=c(-1,1,0))
MB_NEG_vs_Mpos_lrt <- glmLRT(fit, contrast=c(-1,0,1))
Mbneg_vs_Mpos_lrt <- glmLRT(fit, contrast=c(0,-1,1 ))

# MB_NEG compared to Mbneg
topGenes_MB_NEG_Mbneg <- topTags(MB_NEG_vs_Mbneg_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_MB_NEG_Mbneg)
DE_MB_NEG_Mbneg <- topGenes_MB_NEG_Mbneg$table

# MB_NEG compared to Mpos
topGenes_MB_NEG_Mpos<- topTags(MB_NEG_vs_Mpos_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_MB_NEG_Mpos)
DE_topGenes_MB_NEG_Mpos <- topGenes_MB_NEG_Mpos$table

# Mbneg compared to Mpos
topGenes_Mbneg_Mpos <- topTags(Mbneg_vs_Mpos_lrt, adjust.method = "BH", p.value = 1, n=Inf)
dim(topGenes_Mbneg_Mpos)
DE_topGenes_Mbneg_Mpos <- topGenes_Mbneg_Mpos$table

##For visualization
#  Mbpos compared to Mbneg
MB_NEG_Mbneg_de.genes <- decideTests(MB_NEG_vs_Mbneg_lrt, adjust.method = "BH", p.value = 0.05)
MB_NEG_Mpos_de.genes <- decideTests(MB_NEG_vs_Mpos_lrt, adjust.method = "BH", p.value = 0.05)
Mbneg_Mpos_de.genes <- decideTests(Mbneg_vs_Mpos_lrt, adjust.method = "BH", p.value = 0.05)

# summary
summary(MB_NEG_Mbneg_de.genes)
summary(MB_NEG_Mpos_de.genes)
summary(Mbneg_Mpos_de.genes)

# create a dataframe with all the data on differential gene expression
MB_NEG_Mbneg_data <- topGenes_MB_NEG_Mbneg$table
MB_NEG_Mpos_data <- topGenes_MB_NEG_Mpos$table
Mbneg_Mpos_data <- topGenes_Mbneg_Mpos$table

#plotting a volcano plot
MB_NEG_Mbneg_data$Expression_level <- "not significant"
MB_NEG_Mbneg_data$label <- NA
MB_NEG_Mpos_data$Expression_level <- "not significant"
MB_NEG_Mpos_data$label <- NA
Mbneg_Mpos_data$Expression_level <- "not significant"
Mbneg_Mpos_data$label <- NA


# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
MB_NEG_Mbneg$Expression_level[MB_NEG_Mbneg$logFC > 2 & MB_NEG_Mbneg$FDR < 0.05 & MB_NEG_Mbneg$PValue < 0.05] <- "up-regulated"
MB_NEG_Mpos_data$Expression_level[MB_NEG_Mpos_data$logFC > 2 & MB_NEG_Mpos_data$FDR < 0.05 & MB_NEG_Mpos_data$PValue < 0.05] <- "up-regulated"
Mbneg_Mpos$Expression_level[Mbneg_Mpos$logFC > 2 & Mbneg_Mpos$FDR < 0.05 & Mbneg_Mpos$PValue < 0.05] <- "up-regulated"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
MB_NEG_Mbneg_data$Expression_level[MB_NEG_Mbneg_data$logFC < -2 & MB_NEG_Mbneg_data$FDR < 0.05 & MB_NEG_Mbneg_data$PValue < 0.05] <- "down-regulated"
MB_NEG_Mbneg_data <- rownames_to_column(MB_NEG_Mbneg_data, "gene_symbol")
MB_NEG_Mpos_data$Expression_level[MB_NEG_Mpos_data$logFC < -2 & MB_NEG_Mpos_data$FDR < 0.05 & MB_NEG_Mpos_data$PValue < 0.05] <- "down-regulated"
MB_NEG_Mpos_data <- rownames_to_column(MB_NEG_Mpos_data, "gene_symbol")
Mbneg_Mpos_data$Expression_level[Mbneg_Mpos_data$logFC < -2 & Mbneg_Mpos_data$FDR < 0.05 & Mbneg_Mpos_data$PValue < 0.05] <- "down-regulated"
Mbneg_Mpos_data <- rownames_to_column(Mbneg_Mpos_data, "gene_symbol")

#adding labels
MB_NEG_Mbneg_data$label[MB_NEG_Mbneg_data$Expression_level != "not significant"] <- MB_NEG_Mbneg_data$gene_symbol[MB_NEG_Mbneg_data$Expression_level != "not significant"]
MB_NEG_Mpos_data$label[MB_NEG_Mpos_data$Expression_level != "not significant"] <- MB_NEG_Mpos_data$gene_symbol[MB_NEG_Mpos_data$Expression_level != "not significant"]
Mbneg_Mpos_data$label[Mbneg_Mpos_data$Expression_level != "not significant"] <- Mbneg_Mpos_data$gene_symbol[Mbneg_Mpos_data$Expression_level != "not significant"]

##plotting a volcano plot
ggplot(MB_NEG_Mbneg_data, aes(x = logFC, y = -log10(PValue),color = Expression_level, label=label),alpha = 0.7) +
  geom_point()  + theme_minimal() + geom_text() +
  # Scatter plot with color by significance
  scale_color_manual(values = c("blue", "black", "red")) +  # Colors for categories
  theme_minimal() +  # Use a clean theme
  labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 p-value"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # P-value threshold line
  geom_vline(xintercept = c(-2,2), linetype = "dashed") +  # Fold change threshold lines
  theme(legend.position = "top")


#Enrichment analysis
enriched_terms <- read.table("hiddenGoEnrichmentResult.tsv", sep = "\t", header = TRUE)
bar <- ggplot(data = enriched_terms, aes(x = factor(Name), y = Result.count, 
                                         fill=P.value))+ 
  coord_flip()+
  geom_bar(stat = 'identity') +
  labs(x="counts", y="Enriched terms")
bar  

library(STRINGdb)
#gene ontology
read_count1 <- rownames_to_column(read_count, "transcript_id")
enrichment <- string_db$get_enrichment( read_count1)
head(enrichment, n=20)


merged_data <- merge(DE_mbnegmbpos, transcript_ids, by = "Gene_ID", all = TRUE)
merged_data[order(merged_data$FDR, decreasing = TRUE),]  
Mbpos_data <- tibble::rownames_to_column(Mbpos_data, var = "transcript_id")
Mbpos_data <- rownames_to_column(Mbpos_data, "gene_symbol")



