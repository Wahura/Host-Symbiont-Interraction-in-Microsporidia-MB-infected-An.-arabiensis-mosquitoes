getwd()
setwd("/Users/jwahura/Desktop/Proteomics_Mosquito ovaries/anarabiensis")

#loading libraries
library(DEqMS)
library(gplots)
library(ggpubr)
library(ggrepel)

#reading raw data into R
proteomics_data <- read.csv("proteomics_dataset.csv", sep = ',', header = TRUE)
sample_metadata <- read.csv("metadata_V2.csv",sep = ',', header = TRUE)

#getting samples with a qvalue less than 0.001
dat <- proteomics_data %>% filter(Q.value<0.01)
dataset <- dat[,c(1,12:19)]

#converting the accession numbers into rownames
dataset <- dataset %>% remove_rownames %>% column_to_rownames(var="Protein.IDs")

#log transforming the data
dat.log = log2(dataset)
dat.log[sapply(dat.log, is.infinite)] <- NA
dat.log = na.omit(dat.log)

#checking if the data is median centred
boxplot(dat.log,las=2,main="Ovary TMT8plex data")

#normalizing to equal median representation
dat.log = equalMedianNormalization(dat.log)

# Apply sample grouping based on sample_type from which the sample was derived
design <- model.matrix(~0+sample_metadata$condition)
colnames(design) <- levels(as.factor(sample_metadata$condition))

#fitting a linear model
fit1 <- lmFit(dat.log, design)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(as.factor(sample_metadata$condition)), 2))

comparisons <- list()
for (i in 1:nrow(condition_pairs)) {
  comparisons[[i]] <- as.character(condition_pairs[i,])
}

# iterate over the contrasts, and perform a differential expression test for each pair
for (conds in comparisons) {
  # generate string contrast formula
  contrast_formula <- paste(conds, collapse=' - ') 
  
  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=design)
} 

#fitting the model
fit2 <- contrasts.fit(fit1,contrasts = c(-1,1))
fit3 <- eBayes(fit2)

#getting a summary of the expression levels of proteins
summary(decideTests(fit3))

#performing differential abundance analysis
results <- topTable(fit3, adjust = "fdr", sort.by = "P", number = Inf)
results <- rownames_to_column(results, "Protein_ID")

#identifying the differentially expressed proteins
sig_results <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1, ]

#writing the significant results for annotation
write.csv(sig_results, "significant_results.csv")

#Reading the annotated datasets back to R
significant_data <- read.csv("significant_proteins.csv", sep = ',', header = TRUE)

#merging the annotated dataframe with the results dataframe
merged_data <- merge(results, significant_data, by = "Protein_ID", all = TRUE)

#writing the merged table for further annotation
write.csv(merged_data, "merged_dataset.csv")

#filtering proteins with NA values
merged_data <- merged_data %>% filter(!is.na(logFC))

#creating a new significance collumn
merged_data$Significance <- ifelse(abs(merged_data$logFC) > 1 & merged_data$adj.P.Val < 0.05,
                               ifelse(merged_data$logFC > 0, "Upregulated", "Downregulated"),
                               "Non-significant")

# Create volcano plot
volcano_plot <- ggplot(merged_data, aes(x = logFC, y = -log10(adj.P.Val),color = Significance)) +
  geom_point(size = 1, aes(color = abs(logFC) > 1 & adj.P.Val < 0.05)) +
  geom_text(data = subset(merged_data, abs(logFC) > 1 & adj.P.Val < 0.05), 
            aes(label = Protein_name),  
            size = 3, nudge_y = 0.1) +  # Add text labels for significant genes
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted p-value)") +
  theme_minimal() +  
  theme(panel.grid = element_blank(),
axis.line = element_line(color = "black"))+
  scale_color_manual(values = c("red", "blue", "green", "black"),
                     labels = c("Upregulated", "Downregulated", "Non-significant"))

print(volcano_plot)




