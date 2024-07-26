getwd()
setwd("/Users/jwahura/Documents/16s_qpcr/microbiome_analysis/statistical_analysis/")

#loading the required packages
library(tidyverse)
library(naniar)
library(dplyr)
library(seqRFLP)
library(ape)
library(phangorn)
library(seqinr)
library(msa)
library(Biostrings)
library(phyloseq)
library(vegan)
library(microbiome)
library(janitor)
library(metagMisc)
library(gplots)
library(ggpubr)
library(DECIPHER)
library(taxonomizr)
library(edgeR)
library(picante)
library(psych)
library(janitor)
library(schoolmath)
library(igraph)

#Reading datasets into R
taxonomy <- read.csv("taxonomy.tsv", sep = '\t', header = TRUE)
sequence <- read.csv("sequence.tsv", sep = '\t', header = TRUE)
feature_table <- read.csv("feature_table.tsv", sep = '\t', header = TRUE)

#Filtering the unwanted taxons based on taxonomy
filtered_tax <- taxonomy[ !grepl("Mitochondria", taxonomy$Taxon) , ]
filtered_tax <- filtered_tax[ !grepl("Chloroplast", filtered_tax$Taxon) , ]
filtered_tax <- filtered_tax[ !grepl("Unassigned", filtered_tax$Taxon) , ]
filtered_tax <- filtered_tax[ !grepl("Archaea", filtered_tax$Taxon) , ]
filtered_tax <- filtered_tax[ !grepl("Cyanobacteria", filtered_tax$Taxon) , ] 
filtered_tax <- filtered_tax[ !grepl("Chloroflexi", filtered_tax$Taxon) , ] 
filtered_tax <- filtered_tax[ !grepl("Eukaryota", filtered_tax$Taxon) , ]

#separating the taxonomy table to the various taxonomic levels
my_taxonomy <- filtered_tax %>% as.data.frame() %>%
  separate(col = 2,  into = c("Kingdom", "Phylum", "Class", "Order", "Family", "genus","Species"), sep = ";")

#Removing the strings added by silva as identifiers of the various taxonomic levels
my_taxonomy$Kingdom <- sub("d__", "", my_taxonomy$Kingdom)
my_taxonomy$Phylum <- sub("p__", "", my_taxonomy$Phylum)
my_taxonomy$Class <- sub(" c__", "", my_taxonomy$Class)
my_taxonomy$Order <- sub("o__", "", my_taxonomy$Order)
my_taxonomy$Family <- sub("f__", "", my_taxonomy$Family)
my_taxonomy$genus <- sub(" g__", "", my_taxonomy$genus)
my_taxonomy$Species <- sub(" s__", "", my_taxonomy$Species)

#retrieving all genera not NA
silva_classified <- my_taxonomy %>% filter(!is.na(my_taxonomy$genus)) #all OTUs well classified by silva at the genus level
silva_classified <- silva_classified[,-9]
to_blast <- my_taxonomy %>% filter(is.na(my_taxonomy$genus))#all OTUs uncharacterized by silva at the genus level that need to e identified by blast

#merging the OTUs to blast with their corresponding sequences so as to extract the sequences to blast
merged_data <- merge(data, my_taxonomy, by = "Feature_ID", all = FALSE)
#Exracting the sequences for blasting
my_sequences <- merged_data[, c(1,2)]
my_blast_sequences <- dataframe2fas(my_sequences, file = "my_blast_sequences.fasta")

#Reading blast results to R
blast_results <- read.csv("my_blast_results.csv", sep = ',', header = TRUE)

#merging the silva classified taxonomy with the blast classified ones
silva_blast <- as.data.frame(bind_rows(silva_classified, blast_results))

net <- read.csv("ASVs_for_blasting_network_analysis.csv", sep = ',', header = TRUE)
#subsetting the values present for second analysis
subset <- silva_blast[!silva_blast$Feature_ID %in% net$Feature_ID, ]

#preparing a fasta file for blasting of the second batch
network_sequences <- merge(net, sequence, by = "Feature_ID", all = FALSE)
network_blast_sequences <- dataframe2fas(network_sequences, file = "network_blast_sequences.fasta")
#none could be annotated using blast, so excluded from downstream analysis

#merging the taxonomy table with the abundance table
tax_feature <- merge(subset, feature_table, by = "Feature_ID", all = FALSE)

#Extracting the abundance table
feature_table <- tax_feature[,c(-2:-8)]

#Creating a matrix of the taxonomy and feature tables for creation of a phyloseq object
#taxonomy table
my_taxonomy <- silva_blast %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
my_taxonomy <- as.matrix(my_taxonomy)
TAX = tax_table(my_taxonomy)

#feature table
my_feature_table <- feature_table
my_feature_table <- my_feature_table %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
my_feature_table <- as.matrix(my_feature_table)
OTU = otu_table(my_feature_table, taxa_are_rows = TRUE)

#Reading the sample meta data into R
sample_metadata <- read.csv("sample_metadata.csv", sep = ',', header = TRUE)
sdata <- sample_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
samdata = sample_data(sdata)

#creating a phyloseq object
physeq = phyloseq(OTU, TAX, samdata)
physeq

#filtering OTUs with an aubdace of less than five
filtered_physeq <- prune_taxa(taxa_sums(physeq) > 5, physeq)
filtered_physeq

#Extracting the filtered taxonomy and feature tables for barplot plotting
tax_table <- phyloseq_to_df(filtered_physeq, addtax = T, addtot = T, addmaxrank = F)

#ordering data for barplot generation at the genus level
genus_level <- tax_table[,c(7,9:98)]

#grouping the data at the genus level
group <- genus_level %>%
  group_by(genus)%>%
  summarise_each(funs(sum), NBF_10G,PN_BF24_10G,PN_BF48_10G,PN_NBF_11G,NBF_12G,PN_BF24_12G,PN_BF48_12G,BF48_13G,NBF_13G, 
                 PN_BF24_13G,PN_NBF_14G,PN_BF24_15G,NBF_16G,BF48_17G,PN_BF24_17G,BF48_18G,NBF_18G,BF48_19G,PN_NBF_20G,
                 BF48_21G,PN_NBF_21G,BF24_22G,PN_NBF_23G,BF24_25G,PN_BF48_28G,BF48_29G,NBF_2G,PN_BF24_2G,BF48_30G,   
                 PN_NBF_30G,BF24_31G,BF24_32G,BF24_33G,BF48_33G,BF48_34G,PN_NBF_34G,PN_NBF_35G,PN_NBF_36G,BF24_3G, 
                 BF48_3G,NBF_3G,NBF_41G,PN_BF48_41G,PN_NBF_41G,PN_BF48_42G,PN_BF48_43G,NBF_44G,PN_BF48_44G,
                 PN_BF48_45G,BF48_49G,PN_BF48_49G,BF48_4G,NBF_4G,PN_BF24_4G,BF48_51G,BF48_52G,NBF_53G,BF48_54G,BF48_55G,NBF_55G,  
                 BF48_57G,NBF_59G,NBF_5G,PN_BF24_5G,BF24_61G,BF48_61G,BF24_62G,BF48_62G,BF24_63G,BF24_64G,   
                 BF24_66G,NBF_66G,BF24_69G,BF24_6G,NBF_6G,PN_BF24_6G,NBF_70G,BF24_71G,BF24_72G,BF24_74G,BF24_75G,
                 NBF_75G,BF24_79G,BF48_7G,NBF_7G,BF24_83G,BF24_84G,NBF_8G,PN_BF24_9G,PN_BF48_9G,)


#Summing to identify the most abundant genera
cumulation <- group %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

#Extracting the various groups
NBF_internal <- group[,c(1,2,6,10,14,18,42,48,54,64,76)]
NBF_PN <- group[,c(1,5,12,20,22,24,31,37,38,39,45)]
NBF_positive <- group[,c(1,28,43,58,61,63,73,78,83,86,89)]
BF24_positive <- group[,c(1,25,70,71,74,79,80,81,82,87,88)]
BF24_internal <- group[,c(1,23,32,33,34,40,66,68,72,75,84)] 
BF24_PN <- group[,c(1,3,7,11,13,16,29,55,65,77,90)]
BF48_positive <- group[,c(1,27,30,36,51,53,56,57,59,60,62)]
BF48_internal <- group[,c(1,9,15,17,19,21,35,41,67,69,85)] 
BF48_PN <- group[,c(1,4,8,26,44,46,47,49,50,52,91)]

#calculating the average
NBF_internal_total <- NBF_internal %>% adorn_totals(c("col"))
NBF_internal_total <- mutate(NBF_internal_total, NBF_internal=rowSums(NBF_internal_total[12])/10)
NBF_internal_total <- NBF_internal_total[,c(1,13)]

NBF_PN_total <- NBF_PN %>% adorn_totals(c("col"))
NBF_PN_total <- mutate(NBF_PN_total, NBF_PN=rowSums(NBF_PN_total[12])/10)
NBF_PN_total <- NBF_PN_total[,c(1,13)]

NBF_positive_total <- NBF_positive %>% adorn_totals(c("col"))
NBF_positive_total <- mutate(NBF_positive_total, NBF_positive=rowSums(NBF_positive_total[12])/10)
NBF_positive_total <- NBF_positive_total[,c(1,13)]

BF24_positive_total <-BF24_positive %>% adorn_totals(c("col"))
BF24_positive_total <- mutate(BF24_positive_total, BF24_positive=rowSums(BF24_positive_total[12])/10)
BF24_positive_total <- BF24_positive_total[,c(1,13)]

BF24_internal_total <-BF24_internal %>% adorn_totals(c("col"))
BF24_internal_total <- mutate(BF24_internal_total, BF24_internal=rowSums(BF24_internal_total[12])/10)
BF24_internal_total <- BF24_internal_total[,c(1,13)]

BF24_PN_total <-BF24_PN %>% adorn_totals(c("col"))
BF24_PN_total <- mutate(BF24_PN_total, BF24_PN=rowSums(BF24_PN_total[12])/10)
BF24_PN_total <- BF24_PN_total[,c(1,13)]

BF48_positive_total <-BF48_positive %>% adorn_totals(c("col"))
BF48_positive_total <- mutate(BF48_positive_total, BF48_positive=rowSums(BF48_positive_total[12])/10)
BF48_positive_total <- BF48_positive_total[,c(1,13)]

BF48_internal_total <-BF48_internal %>% adorn_totals(c("col"))
BF48_internal_total <- mutate(BF48_internal_total, BF48_internal=rowSums(BF48_internal_total[12])/10)
BF48_internal_total <- BF48_internal_total[,c(1,13)]

BF48_PN_total <-BF48_PN %>% adorn_totals(c("col"))
BF48_PN_total <- mutate(BF48_PN_total, BF48_PN=rowSums(BF48_PN_total[12])/10)
BF48_PN_total <- BF48_PN_total[,c(1,13)]

#merging the various dataframes
merged <- Reduce(function(x,y) merge(x,y,by="genus",all=TRUE) ,list(NBF_positive_total,NBF_internal_total,NBF_PN_total,BF24_positive_total,BF24_internal_total,BF24_PN_total,BF48_positive_total,BF48_internal_total,BF48_PN_total))

#chosing the genera to represent as the most abundant
to_represent <- c("Rhodococcus","Pseudomonas","Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Serratia", "Achromobacter","Elizabethkingia","Brevundimonas","Ralstonia", "Acidovorax", "Variovorax", "Phyllobacterium","Mesorhizobium","Aeromonas", "Asaia", "Acinetobacter", "Rahnella1", "Chryseobacterium", "Enterococcus", "Delftia")

#aggregating the rest of the phyla as others
grouped_data <- aggregate(merged[-1], list(genus = replace(merged$genus,!(merged$genus %in% to_represent), "Others")), sum)
View(grouped_data) 

#converting the abudance into percentage
bar <- adorn_percentages(grouped_data, denominator = "col", na.rm = TRUE)

#Gathering the data
bar <- bar %>%
  gather(value = "abundance", key = "sample_names", -genus)

#ordering the data on abundance basis for plotting
bar$genus <- reorder(bar$genus, bar$abundance)
bar$genus <- factor(bar$genus, levels=rev(levels(bar$genus)))

bar$genus <- reorder(bar$genus, bar$abundance)
bar$genus <- factor(bar$genus, levels=rev(levels(bar$genus)))
bar$genus <- factor(bar$genus, levels=c("Rhodococcus","Pseudomonas","Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Serratia", "Achromobacter","Elizabethkingia","Brevundimonas","Ralstonia", "Acidovorax", "Variovorax", "Phyllobacterium","Mesorhizobium","Aeromonas", "Asaia", "Acinetobacter", "Rahnella1", "Chryseobacterium", "Enterococcus", "Delftia","Others"))


#specifying the genera to be italicized
lbs = brk = levels(bar$genus)
lbs[match("Rhodococcus", brk)] = expression(italic("Rhodococcus"))
lbs[match("Pseudomonas", brk)] = expression(italic("Pseudomonas"))
lbs[match("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", brk)] = expression(italic("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"))
lbs[match("Serratia", brk)] = expression(italic("Serratia"))
lbs[match("Achromobacter", brk)] = expression(italic("Achromobacter"))
lbs[match("Elizabethkingia", brk)] = expression(italic("Elizabethkingia"))
lbs[match("Brevundimonas", brk)] = expression(italic("Brevundimonas"))
lbs[match("Ralstonia", brk)] = expression(italic("Ralstonia"))
lbs[match("Acidovorax", brk)] = expression(italic("Acidovorax"))
lbs[match("Variovorax", brk)] = expression(italic("Variovorax"))
lbs[match("Phyllobacterium", brk)] = expression(italic("Phyllobacterium"))
lbs[match("Mesorhizobium", brk)] = expression(italic("Mesorhizobium"))
lbs[match("Aeromonas", brk)] = expression(italic("Aeromonas"))
lbs[match("Asaia", brk)] = expression(italic("Asaia"))
lbs[match("Acinetobacter", brk)] = expression(italic("Acinetobacter"))
lbs[match("Rahnella1", brk)] = expression(italic("Rahnella1"))
lbs[match("Chryseobacterium", brk)] = expression(italic("Chryseobacterium"))
lbs[match("Enterococcus", brk)] = expression(italic("Enterococcus"))
lbs[match("Delftia", brk)] = expression(italic("Delftia"))
lbs[match("Others", brk)] = expression(plain("Others"))

#pallet
myPalette <- c('#89C5DA', "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#plotting the barplot
p <- ggplot(bar,aes(x = fct_inorder(sample_names), y = abundance), labs(fill= genus), group=row.names(bar))+ xlab("Gut Tissue source")+ ylab("Abundance") + geom_col(aes(fill = genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_fill_manual(values = myPalette, labels = lbs, breaks = brk)+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10)+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = "grey"))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90))+
          theme(axis.title.x = element_text(size = 10, angle = 0)))
p + ylab("Abundance") + scale_y_continuous(labels = scales::percent_format(accuracy = 1))+ labs(tag = "A", plot.tag.position = c(0.2, -0.1))

#plotting venn diagrams
#Removing all genera with zero hits in the dataframes
NBFinternal <- NBF_internal_total[!(NBF_internal_total$NBF_internal == 0),]
NBFPN <- NBF_PN_total[!(NBF_PN_total$NBF_PN == 0),]
NBFpositive <- NBF_positive_total[!(NBF_positive_total$NBF_positive == 0),]
BF24positive <- BF24_positive_total[!(BF24_positive_total$BF24_positive == 0),]
BF24internal <- BF24_internal_total[!(BF24_internal_total$BF24_internal == 0),]
BF24PN <- BF24_PN_total[!(BF24_PN_total$BF24_PN == 0),]
BF48positive <- BF48_positive_total[!(BF48_positive_total$BF48_positive == 0),]
BF48internal <- BF48_internal_total[!(BF48_internal_total$BF48_internal == 0),]
BF48PN <- BF48_PN_total[!(BF48_PN_total$BF48_PN == 0),]

NBFinternal <- NBFinternal %>% column_to_rownames(var = "genus")
NBFPN <- NBFPN %>% column_to_rownames(var = "genus")
NBFpositive <- NBFpositive %>% column_to_rownames(var = "genus")
BF24positive <- BF24positive %>% column_to_rownames(var = "genus")
BF24internal <- BF24internal %>% column_to_rownames(var = "genus")
BF24PN <- BF24PN %>% column_to_rownames(var = "genus")
BF48positive <- BF48positive %>% column_to_rownames(var = "genus")
BF48internal <- BF48internal %>% column_to_rownames(var = "genus")
BF48PN <- BF48PN %>% column_to_rownames(var = "genus")

# get the Genus vector which is in the rownames
NBFinternal = rownames(NBFinternal)
NBFPN = rownames(NBFPN)
NBFpositive = rownames(NBFpositive)
BF24positive = rownames(BF24positive)
BF24internal = rownames(BF24internal)
BF24PN = rownames(BF24PN)
BF48positive = rownames(BF48positive)
BF48internal = rownames(BF48internal)
BF48PN = rownames(BF48PN)

vd <- list(BF48positive, BF48internal, BF48PN)
names(vd) = c("BF48 positive", "BF48 internal", "BF 48PN")
venn(vd)

require(VennDiagram)


vp <- venn.diagram(vd, 
                   fill = 1:3, alpha = 0.3, filename = NULL);

grid.draw(vp)


#function for identifying the common elements
Intersect <- function (vd) {  
  # Multiple set version of intersect
  # vd is a list
  if (length(vd) == 1) {
    unlist(vd)
  } else if (length(vd) == 2) {
    intersect(vd[[1]], vd[[2]])
  } else if (length(vd) > 2){
    intersect(vd[[1]], Intersect(vd[-1]))
  }
}

Intersect(vd)

#Beta diversity estimation
#Extracting sequences to be included in the study for plotting phylogenetic trees
seq <- merge(sequence, tax_table, by.x = "Feature_ID", by.y = "OTU", all = FALSE)
seq <- seq[, c(1,2)]

#converting the filtered sequences to fasta format ad writing the fasta file to the working directory
my_sequences <- dataframe2fas(seq, file = "my_mosquito_tree_sequences.fasta")

#Reading the sequences back to R
my_sequences <- readDNAStringSet("my_mosquito_tree_sequences.fasta")
names(my_sequences)

#Running multiple sequence alignment
alignment <- AlignSeqs(my_sequences, anchor = NA)

#constructng the phylogenetic tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)
fit = pml(treeNJ, data=phang.align)

#fitting the tree using the GTR model
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) 

saveRDS(fitGTR, "mosquito_phangorn_tree.RDS")
phangorn <- readRDS("mosquito_phangorn_tree.RDS")

#Extracting the tree from the GTR model
phylo_tree <- phangorn$tree
phylogenetic_tree <- phy_tree(phylo_tree)
phylogenetic_tree <- root(phy_tree(phylo_tree), sample(taxa_names(phylo_tree), 1), resolve.root = TRUE)
is.rooted(phylogenetic_tree)

#separating the samples per the treatment conditions
NBF <- feature_table[,c(1,2,5,6,10,12,14,18,20,22,24,28,31,37,38,39,42,43,45,48,54,58,61,63,64,73,76,78,83,86,89)]
BF24 <- feature_table[,c(1,3,7,11,13,16,23,25,29,32,33,34,40,55,65,66,68,70,71,72,74,75,77,79,80,81,82,84,87,88,90)]
BF48 <- feature_table[,c(1,4,8,9,15,17,19,21,26,27,30,35,36,41,44,46,47,49,50,51,52,53,56,57,59,60,62,67,69,85,91)]
MBpos <- feature_table[,c(1,28,86,89,43,58,61,73,78,83,63,53,27,30,51,56,62,59,36,60,57,25,70,71,74,79,87,88,81,82,80)]
MB_internal <- feature_table[,c(1,42,54,64,76,2,6,10,14,18,48,23,33,34,66,68,72,84,40,75,32,41,19,21,35,9,17,67,85,15,69)]
PN_controls <- feature_table[,c(1,5,12,20,22,24,31,37,38,39,45,3,7,11,13,16,29,55,65,77,90,4,8,26,44,46,47,49,50,52,91)]
  
  
sdata_NBF <- sample_metadata[c(1:30),]
sdata_BF24 <- sample_metadata[c(31:60),]
sdata_BF48 <- sample_metadata[c(61:90),]
sdata_MBpos <- sample_metadata[c(1:10,31:40,61:70),]
sdata_MBinternal <- sample_metadata[c(11:20,41:50,71:80),]
sdata_PN <- sample_metadata[c(21:30,51:60,81:90),]

#creating the corresponding phyloseq objects
#feature table
NBF_feature_table <- NBF %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
NBF_feature_table <- as.matrix(NBF_feature_table)
NBF_OTU = otu_table(NBF_feature_table, taxa_are_rows = TRUE)

BF24_feature_table <- BF24 %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
BF24_feature_table <- as.matrix(BF24_feature_table)
BF24_OTU = otu_table(BF24_feature_table, taxa_are_rows = TRUE)

BF48_feature_table <- BF48 %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
BF48_feature_table <- as.matrix(BF48_feature_table)
BF48_OTU = otu_table(BF48_feature_table, taxa_are_rows = TRUE)

MBpos_feature_table <- MBpos %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
MBpos_feature_table <- as.matrix(MBpos_feature_table)
MBpos_OTU = otu_table(MBpos_feature_table, taxa_are_rows = TRUE)

MB_internal_feature_table <- MB_internal %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
MB_internal_feature_table <- as.matrix(MB_internal_feature_table)
MB_internal_OTU = otu_table(MB_internal_feature_table, taxa_are_rows = TRUE)

PN_controls_feature_table <- PN_controls %>% remove_rownames %>% column_to_rownames(var="Feature_ID")
PN_controls_feature_table <- as.matrix(PN_controls_feature_table)
PN_controls_OTU = otu_table(PN_controls_feature_table, taxa_are_rows = TRUE)
#sample meta data 
sdata_NBF <- sdata_NBF %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
NBF_samdata = sample_data(sdata_NBF)

sdata_BF24 <- sdata_BF24 %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
BF24_samdata = sample_data(sdata_BF24)

sdata_BF48 <- sdata_BF48 %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
BF48_samdata = sample_data(sdata_BF48)

sdata_MBpos <- sdata_MBpos %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
MBpos_samdata = sample_data(sdata_MBpos)

sdata_MBinternal <- sdata_MBinternal %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
MBinternal_samdata = sample_data(sdata_MBinternal)

sdata_PN <- sdata_PN %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
PN_samdata = sample_data(sdata_PN)

#creating a phyloseq object
NBF_physeq = phyloseq(NBF_OTU, TAX, NBF_samdata,phylogenetic_tree)
NBF_physeq

BF24_physeq = phyloseq(BF24_OTU, TAX, BF24_samdata,phylogenetic_tree)
BF24_physeq

BF48_physeq = phyloseq(BF48_OTU, TAX, BF48_samdata,phylogenetic_tree)
BF48_physeq

MBpos_physeq = phyloseq(MBpos_OTU, TAX, MBpos_samdata,phylogenetic_tree)
MBpos_physeq

MBinternal_physeq = phyloseq(MB_internal_OTU, TAX, MBinternal_samdata,phylogenetic_tree)
MBinternal_physeq

PN_physeq = phyloseq(PN_controls_OTU, TAX, PN_samdata,phylogenetic_tree)
PN_physeq


mypallet =c("#F52100", "#0071BC","#FBAE17", "#65BD62","#3F4921",'#89C5DA', "#DA5724", "#74D944","#CE50CA","#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E","#508578")
#Weighted unifrac
ordu = ordinate(PN_physeq, "PCoA", "unifrac", weighted = TRUE)
p <- plot_ordination(PN_physeq, ordu, color="Treatment")+ geom_point(size=2) +
  scale_color_manual(values = mypallet) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + stat_ellipse()+labs(tag = "A", plot.tag.position = c(0.2, -0.1))

#unweighted unifrac
ordu = ordinate(PN_physeq, "PCoA", "unifrac", weighted = FALSE)
p <- plot_ordination(PN_physeq, ordu, color="Treatment")+ geom_point(size=2) +
  scale_color_manual(values = mypallet) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + stat_ellipse()+labs(tag = "B", plot.tag.position = c(0.2, -0.1))

#ordinating the phyloseq object
ordu = ordinate(PN_physeq, "PCoA", "bray")
p <- plot_ordination(PN_physeq, ordu, color="Treatment")+ geom_point(size=2) +
  scale_color_manual(values = mypallet) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + stat_ellipse() + labs(tag = "C", plot.tag.position = c(0.2, -0.1))


#alpha diversity estimation
#checking out the total read counts in the samples
reads <- sample_sums(filtered_physeq)
reads

summary(sample_sums(filtered_physeq))

#Extracting the otu table from the phyloseq object and plotting the rarefaction curve
otu_tab <- t(abundances(filtered_physeq))
raremax <- 30000


#Extracting the otu table from the phyloseq object and plotting the rarefaction curve
otu_tab <- t(abundances(filtered_physeq))
head(otu_tab)

otu_df <- as.data.frame(otu_tab)

#plotting the rarefaction curve
r <- rarecurve(otu_df, step=10000, lwd=2, ylab="OTU",  label=F)

Nmax <- sapply(r, function(x) max(attr(x, "Subsample")))
Smax <- sapply(r, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "number of reads/sample",
     ylab = "number of OTUs", type = "n")
abline(v = raremax)
for (i in seq_along(r)) {
  N <- attr(r[[i]], "Subsample")
  lines(N, r[[i]], col = "black")
}

set.seed(9242)
#calculatin an even sampling depth for all the samples
rarefied <- rarefy_even_depth(filtered_physeq, sample.size = 30000)
rarefied

#calculating the alpha diversity
diversity <- alpha(rarefied, index = "all")
diversity <- rownames_to_column(diversity, "Sample_ID")
diversity <- merge(sample_metadata, diversity, by= "Sample_ID", all= TRUE)
diversity <- read.csv("diversity.csv", sep = ',', header = TRUE)
NBF_diversity <- diversity[c(1:30),]
BF24_diversity <- diversity[c(31:60),]
BF48_diversity <- diversity[c(61:90),]
MBpos <- diversity[c(1:10,31:40,61:70),]
MB_Internal <- diversity[c(11:20,41:50,71:80),]
PN <- diversity[c(21:30,51:60,81:90),]



#shannon diversity
shannon <- PN[, c(1:3,8)]
P <- ggboxplot(shannon, "MB_status","diversity_shannon",
               color = "MB_status",
               add = "jitter", linetype = "solid", Family = "Palatino Linotype", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10)+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = "grey"))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90))+
          theme(axis.title.x = element_text(size = 10, angle = 0))) + stat_compare_means()
P   + theme(legend.position = "none") + xlab("Gut Tissue Source") + ylab("Shannon index") + labs(tag = "A", plot.tag.position = c(0.2, -0.1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))

#chao1 diversity
chao1 <- PN[,c(1:3,5)]
P <- ggboxplot(chao1, "MB_status","chao1",
               color = "MB_status", 
               add = "jitter", linetype = "solid", Family = "Palatino Linotype", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10)+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = "grey"))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90))+
          theme(axis.title.x = element_text(size = 10, angle = 0))) + stat_compare_means()
P + theme(legend.position = "none") + xlab("Gut Tissue Source") + ylab("chao1 index")+ labs(tag = "B", plot.tag.position = c(0.2, -0.1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))

#Evenness
Evenness <- PN[,c(1:3,13)]
P <- ggboxplot(Evenness, "MB_status","evenness_simpson",
               color = "MB_status",
               add = "jitter", linetype = "solid", Family = "Palatino Linotype", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10)+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = "grey"))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90))+
          theme(axis.title.x = element_text(size = 10, angle = 0))) + stat_compare_means()
P + theme(legend.position = "none") + xlab("Gut Tissue Source") + ylab("Evenness")+ labs(tag = "C", plot.tag.position = c(0.2, -0.1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))

#Faith's phylogeetic diversity
FD_OTU <- as.data.frame(rarefied@otu_table)
df.pd <- pd(t(FD_OTU), phylogenetic_tree,include.root=T)
df.pd <- rownames_to_column(df.pd, "Sample_ID")
df.pd_edited <- merge(sample_metadata, df.pd, by = "Sample_ID", all = FALSE)

df.pd_edited <- read.csv("df.pd.edited.csv", sep = ',', header= TRUE)

P <- ggboxplot(df.pd_edited, "MB_status","PD",
               color = "MB_status",
               add = "jitter", linetype = "solid", Family = "Palatino Linotype", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
  theme(legend.text = element_text(size = 10, colour = "black"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10)+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = "grey"))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90))+
          theme(axis.title.x = element_text(size = 10, angle = 0))) + stat_compare_means()
P + theme(legend.position = "none") + xlab("Gut Tissue Source") + ylab("Faith's phylogenetic diversity (PD)")+ labs(tag = "D", plot.tag.position = c(0.2, -0.1))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, family = "Arial"))

#Prevalence analysis
# make the genus column rownames
prevalence_estimation <- group %>% column_to_rownames("genus")

# Count of present bacteria 
apply(prevalence_estimation,1, function(x) sum(x != 0))

# Count of present bacteria expressed as percentage
apply(prevalence_estimation,1, function(x) sum(x != 0)/ncol(prevalence_estimation)*100)

# extract bacteria with a certain percentage of prevalence
prevalence_value <- 90
test <- apply(prevalence_estimation,1, function(x) sum(x != 0)/ncol(prevalence_estimation)*100)
rownames(prevalence_estimation[test > prevalence_value,])

#creating a dataframe with the prevalence values for all the genus members
prevalence <- as.data.frame(apply(prevalence_estimation,1, function(x) sum(x != 0)/ncol(prevalence_estimation)*100))

#Permanova pairwise analysis
#Calculate bray curtis distance matrix
physeq_bray <- phyloseq::distance(filtered_physeq, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(filtered_physeq))

# Adonis test
permanova <- adonis(physeq_bray ~ MB_status, data = sampledf)
permanova
# P-value
print(as.data.frame(permanova$aov.tab)["MB_status", "Pr(>F)"])

source("parwise.adonis.r")
# convert OTU ID and sample id to rownames
permanova_feature <- group %>% column_to_rownames(var= "genus")
#sample_metadata <- sample_metadata %>% column_to_rownames(var= "sample.id")

# transpose the feature table
permanova_feature_transposed <- t(permanova_feature)

# create sample.ids column from rownames
permanova_feature_transposed <- permanova_feature_transposed %>% as.data.frame() %>% rownames_to_column(var = "Sample_ID")

# merge the sample metadata with the feature table using the sample ids
permanova_feature_transposed <- merge(sample_metadata, permanova_feature_transposed, by = "Sample_ID", all = FALSE)

# perform permanova analysis 
permanova <- pairwise.adonis(permanova_feature_transposed[,!names(permanova_feature_transposed) %in% c("Sample_ID", "Treatment","MB_status")],
                             factors = permanova_feature_transposed$MB_status,
                             sim.method = 'bray',
                             p.adjust.m ='bonferroni')
permanova

#Differential abundance testing
# Convert phyloseq object to edgeR object
otu_table <- group
otu_table <- otu_table %>% remove_rownames %>% column_to_rownames(var="genus")
otu_table <- as.matrix(otu_table)

#grouping factor
sdat <- factor(sample_metadata$MB_status)

# Ensure OTU table has non-zero dimensions
if (dim(otu_table)[1] > 0 && dim(otu_table)[2] > 0) {
  dge <- DGEList(counts = otu_table, group = sdat)
} else {
  stop("OTU table has zero dimensions. Please check your input data.")
}

# Inspect the DGEList object
str(dge)

# Calculate normalization factors
dge <- calcNormFactors(dge)

# Apply sample grouping based on sample_type from which the sample was derived
design <- model.matrix(~0+sample_metadata$MB_status)
colnames(design) <- levels(as.factor(sample_metadata$MB_status))

# Estimate dispersion
dge <- estimateDisp(dge,design, robust = TRUE)

# Perform differential abundance testing
fit <- glmFit(dge, design)

# generate a list of all possible pairwise contrasts
condition_pairs <- t(combn(levels(as.factor(sample_metadata$MB_status)), 2))

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
de.genes <- dge[rownames(dge) %in% sig_genes,]

dim(de.genes$counts)
DE_data <- de.genes$counts

order <- c(BF24_internal, BF24_PN,BF24_positive, BF48_negative,BF48_PN,BF48_positive,NBF_internal,NBF_PN,NBF_positive)

# comparison of treatments across time points
BF24_internal_PN <- glmLRT(fit, contrast=c(1,-1,0,0,0,0,0,0,0))
BF24_internal_positive <- glmLRT(fit, contrast=c(-1,0,1,0,0,0,0,0,0))
BF24_PN_positive <- glmLRT(fit, contrast=c(0,-1,1,0,0,0,0,0,0))

BF48_internal_PN <- glmLRT(fit, contrast=c(0,0,0,1,-1,0,0,0,0))
BF48_internal_positive <- glmLRT(fit, contrast=c(0,0,0,-1,0,1,0,0,0))
BF48_PN_positive <- glmLRT(fit, contrast=c(0,0,0,0,-1,1,0,0,0))

NBF_internal_PN <- glmLRT(fit, contrast=c(0,0,0,0,0,0,1,-1,0))
NBF_internal_positive <- glmLRT(fit, contrast=c(0,0,0,0,0,0,-1,0,1))
NBF_PN_positive <- glmLRT(fit, contrast=c(0,0,0,0,0,0,0,-1,1))

MBpos_24_48 <- glmLRT(fit, contrast=c(0,0,-1,0,0,1,0,0,0))
MBpos_24_NBF <- glmLRT(fit, contrast=c(0,0,1,0,0,0,0,0,-1))
MBpos_48_NBF <- glmLRT(fit, contrast=c(0,0,0,0,0,1,0,0,-1))

internal_24_48 <- glmLRT(fit, contrast=c(-1,0,0,1,0,0,0,0,0))
internal_24_NBF <- glmLRT(fit, contrast=c(1,0,0,0,0,0,-1,0,0))
internal_48_NBF <- glmLRT(fit, contrast=c(0,0,0,1,0,0,-1,0,0))

PN_24_48 <- glmLRT(fit, contrast=c(0,-1,0,0,1,0,0,0,0))
PN_24_NBF <- glmLRT(fit, contrast=c(0,1,0,0,0,0,0,-1,0))
PN_48_NBF <- glmLRT(fit, contrast=c(0,0,0,0,1,0,0,-1,0))


# Mbpos compared to Mbneg
BF24_internal_PN_topGenes <- topTags(BF24_internal_PN, adjust.method = "BH", p.value = 1, n=Inf)
dim(BF24_internal_PN_topGenes)
DE_BF24_internal_PN <- BF24_internal_PN_topGenes$table
write.csv(DE_BF24_internal_PN, "BF24_internal_PN.csv")

BF24_internal_positive_topGenes <- topTags(BF24_internal_positive, adjust.method = "BH", p.value = 1, n=Inf)
dim(BF24_internal_positive_topGenes)
DE_BF24_internal_positive <- BF24_internal_positive_topGenes$table
write.csv(DE_BF24_internal_positive, "BF24_internal_positive.csv")

BF24_PN_positive_topGenes <- topTags(BF24_PN_positive, adjust.method = "BH", p.value = 1, n=Inf)
dim(BF24_PN_positive_topGenes)
DE_BF24_PN_positive <- BF24_PN_positive_topGenes$table
write.csv(DE_BF24_PN_positive, "DE_BF24_PN_positive.csv")

BF48_internal_PN_topGenes <- topTags(BF48_internal_PN, adjust.method = "BH", p.value = 1, n=Inf)
dim(BF48_internal_PN_topGenes)
DE_BF48_internal_PN <- BF48_internal_PN_topGenes$table
write.csv(DE_BF48_internal_PN, "BF48_internal_PN.csv")

BF48_internal_positive_topGenes <- topTags(BF48_internal_positive, adjust.method = "BH", p.value = 1, n=Inf)
dim(BF48_internal_positive_topGenes)
DE_BF48_internal_positive <- BF48_internal_positive_topGenes$table
write.csv(DE_BF48_internal_positive, "BF48_internal_positive.csv")

BF48_PN_positive_topGenes <- topTags(BF48_PN_positive, adjust.method = "BH", p.value = 1, n=Inf)
dim(BF48_PN_positive_topGenes)
DE_BF48_PN_positive <- BF48_PN_positive_topGenes$table
write.csv(DE_BF48_PN_positive, "DE_BF48_PN_positive.csv")

NBF_internal_PN_topGenes <- topTags(NBF_internal_PN, adjust.method = "BH", p.value = 1, n=Inf)
dim(NBF_internal_PN_topGenes)
DE_NBF_internal_PN <- NBF_internal_PN_topGenes$table
write.csv(DE_NBF_internal_PN, "NBF_internal_PN.csv")

NBF_internal_positive_topGenes <- topTags(NBF_internal_positive, adjust.method = "BH", p.value = 1, n=Inf)
dim(NBF_internal_positive_topGenes)
DE_NBF_internal_positive <- NBF_internal_positive_topGenes$table
write.csv(DE_NBF_internal_positive, "NBF_internal_positive.csv")

NBF_PN_positive_topGenes <- topTags(NBF_PN_positive, adjust.method = "BH", p.value = 1, n=Inf)
dim(NBF_PN_positive_topGenes)
DE_NBF_PN_positive <- NBF_PN_positive_topGenes$table
write.csv(DE_NBF_PN_positive, "DE_NBF_PN_positive.csv")

MBpos_24_48_topGenes <- topTags(MBpos_24_48, adjust.method = "BH", p.value = 1, n=Inf)
dim(MBpos_24_48_topGenes)
DE_MBpos_24_48 <- MBpos_24_48_topGenes$table
write.csv(DE_MBpos_24_48, "DE_MBpos_24_48_positive.csv")

MBpos_24_NBF_topGenes <- topTags(MBpos_24_NBF, adjust.method = "BH", p.value = 1, n=Inf)
dim(MBpos_24_NBF_topGenes)
DE_MBpos_24_NBF <- MBpos_24_48_topGenes$table
write.csv(DE_MBpos_24_NBF, "DE_MBpos_24_NBF.csv")

MBpos_48_NBF_topGenes <- topTags(MBpos_48_NBF, adjust.method = "BH", p.value = 1, n=Inf)
dim(MBpos_48_NBF_topGenes)
DE_MBpos_48_NBF <- MBpos_24_48_topGenes$table
write.csv(DE_MBpos_48_NBF, "DE_MBpos_48_NBF.csv")

internal_24_48_topGenes <- topTags(internal_24_48, adjust.method = "BH", p.value = 1, n=Inf)
dim(internal_24_48_topGenes)
DE_internal_24_48 <- internal_24_48_topGenes$table
write.csv(DE_internal_24_48, "DE_internal_24_48.csv")

internal_24_NBF_topGenes <- topTags(internal_24_NBF, adjust.method = "BH", p.value = 1, n=Inf)
dim(internal_24_NBF_topGenes)
DE_internal_24_NBF <- internal_24_NBF_topGenes$table
write.csv(DE_internal_24_NBF, "DE_internal_24_NBF.csv")

internal_48_NBF_topGenes <- topTags(internal_48_NBF, adjust.method = "BH", p.value = 1, n=Inf)
dim(internal_48_NBF_topGenes)
DE_internal_48_NBF <- internal_48_NBF_topGenes$table
write.csv(DE_internal_48_NBF, "DE_internal_48_NBF.csv")

PN_24_48_topGenes <- topTags(PN_24_48, adjust.method = "BH", p.value = 1, n=Inf)
dim(PN_24_48_topGenes)
DE_PN_24_48 <- PN_24_48_topGenes$table
write.csv(DE_PN_24_48, "DE_PN_24_48.csv")

PN_24_NBF_topGenes <- topTags(PN_24_NBF, adjust.method = "BH", p.value = 1, n=Inf)
dim(PN_24_NBF_topGenes)
DE_PN_24_NBF <- PN_24_NBF_topGenes$table
write.csv(DE_PN_24_NBF, "DE_PN_24_NBF.csv")

PN_48_NBF_topGenes <- topTags(PN_48_NBF, adjust.method = "BH", p.value = 1, n=Inf)
dim(PN_48_NBF_topGenes)
DE_PN_48_NBF <- PN_48_NBF_topGenes$table
write.csv(DE_PN_48_NBF, "DE_PN_48_NBF.csv")

#generating a correlation network
#specifying the dataset to work with
dataframe <- merged

#Transposing the dataframe
data_frame <- as.data.frame(t(dataframe))

#converting the first row into column names
data_frame <- data_frame %>%
  row_to_names(row_number = 1)

#making the dataframe numeric
data_frame <- data_frame %>% mutate_all(as.numeric) %>% as_tibble()

#Reading data into R
dataframe <- read.csv("network_data.csv", sep = ",", header = TRUE)

#Removing the rownames
data_frame <- dataframe %>%
  column_to_rownames(var = "samples")


#data normality test
shapiro <- lapply(as.vector(data_frame), shapiro.test)
shapiro_test <- sapply(shapiro, `[`, c("statistic","p.value"))
shapiro_test

#calculating pearson's correlation
correlation_matrix=corr.test(data_frame, method="pearson", adjust="fdr")

# assign symmetric matrix of r-values
r_matrix=correlation_matrix$r

# replace all values in the bottom triangle of the r-matrix to zeros
r_matrix[row(r_matrix)>col(r_matrix)]=0

# replace all values in the diagonal of the r-matrix to zeros
diag(r_matrix)=0

# assign asymmetric matrix of p-values 
p_matrix=correlation_matrix$p

# replace all values in the bottom triangle of the p-matrix to one
p_matrix[row(p_matrix)>col(p_matrix)]=1

# replace all values in the diagonal of the r-matrix to one
diag(p_matrix)=1

# load "reshape2" package
library(reshape2)

# transform r-matrix to r-table with three columns
r_table=melt(r_matrix)

# transform p-matrix to p-table with three columns
p_table=melt(p_matrix)

# save obtained r and p tables to the corresponded file
write.csv(r_table, "network_r_table.csv", row.names = FALSE) 
write.csv(p_table, "network_p_table.csv", row.names = FALSE)

# rename the "value" column in r_table
colnames(r_table) <- c("Var1", "Var2", "correlation")

# add the p-value column to the r_table from the p_table
r_table$pvalue <- p_table$value

# filter the r_table for rows with absolute correlation > 0.65 and pvalue < 0.05
r_table_filtered <- r_table %>% filter(abs(correlation) > 0.65 & pvalue < 0.05)

# save the filtered r_table to csv
write.csv(r_table_filtered, "network_r_table_filtered.csv", row.names = FALSE) 

# Node list
nodes <- c(as.vector(r_table_filtered$Var1), as.vector(r_table_filtered$Var2))
nodes <- unique(nodes)
nodes <- as.data.frame(nodes)

# create a unique id for each node name
nodes <- nodes %>% rowid_to_column("id")

# rename the "nodes" column to "label"
colnames(nodes) <- c("id", "label")

# Edges List
edges <- r_table_filtered %>% 
  left_join(nodes, by = c("Var1" = "label"))  

colnames(edges) <- c("Var1", "Var2", "correlation", "pvalue", "from")

edges <- edges %>% 
  left_join(nodes, by = c("Var2" = "label")) 

colnames(edges) <- c("Var1", "Var2", "correlation", "pvalue", "from", "to")  

# select the columns needed to construct the network
edges <- select(edges, from, to, correlation)

# create an igraph object
routes_igraph <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

# plot the network
png(filename = "MB_network.png", type = "cairo", res =700, units = 'in',
    width = 15, height = 20, pointsize = 9)

pdf(filename = "MB_network.pdf", type = "cairo", res =700, units = 'in',
    width = 15, height = 20, pointsize = 9)
plot(routes_igraph, 
     edge.arrow.size = 0.8,
     edge.width=1,
     vertex.label.dist=2,
     vertex.label.degree=15,
     vertex.label.cex=1.05, 
     edge.color = "black",
     edge.lty=ifelse(is.positive(r_table_filtered$correlation), "solid","dashed"))

dev.off()

pdf(file = "bact_network.pdf",
    width = 7, height = 6, pointsize = 9)

plot(routes_igraph, 
     edge.arrow.size = 0.8,
     edge.width=1,
     vertex.label.dist=2,
     vertex.label.degree=15,
     vertex.label.cex=1.05, 
     edge.color = "black",
     edge.lty=ifelse(is.positive(r_table_filtered$correlation), "solid","dashed"))

dev.off() %>%
  row_to_names(row_number = 1)









