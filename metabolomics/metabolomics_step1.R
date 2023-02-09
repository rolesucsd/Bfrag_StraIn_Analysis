#################
# Title: metabolomics_step1.R
# Author- Renee Oles
# Purpose: Download and separate inititial metabolomics data
# Output: full metabolomics raw matrix and PCA plots after normalization
# Date- 1/12/2023
# TODO: ADD IN METADATA ONCE IT's FIXED
################

# Libraries
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(plotly)
library(Hotelling)
library(KEGGREST)

# IMPORT AND DATA-WRANGLING
#####
# Input
anno <-read.csv("resources/bfrag_quant.csv")
hits <- read_tsv("resources/top_hits.tsv")
# TODO: UPDATE METADATA - LUKE 
# Phylo is division/other information
phylo <- read.table("resources/meta_update.txt",sep="\t", header=T)
# Metadata is general sample and assignment information
metadata <- read.table("resources/metadata.txt", sep="\t", header=T, fill=T)

# Connect to KEGG pathways
pathways <- unique(hits$KEGG)
kegg_pathways <- keggLink("pathway", pathways)
kegg_pathways <- data.frame("KEGG"=names(kegg_pathways), "PATH"=kegg_pathways)
kegg_pathways$KEGG <- gsub("cpd:","",kegg_pathways$KEGG)
kegg_pathways$PATH <- gsub("path:","",kegg_pathways$PATH)
kegg_pathways <- as.data.frame(kegg_pathways)
write.table(kegg_pathways, "full_results/kegg_pathways.txt", sep="\t", quote=F, row.names = F)

# Write hits pathway
hits <- left_join(hits,kegg_pathways)
hits_sub <- hits[,c(1,6,3,45,36:42)]
colnames(hits_sub)[1] <- "Scan"
write.table(hits_sub, "full_results/hits_pathways.txt", sep="\t", quote=F, row.names = F)

# Make matrix with anno dataframe
colnames(anno)[1] <- "Scan"
anno <- right_join(hits_sub,anno)
anno_sub <- anno[!duplicated(anno$Scan),]
anno_sub <- as.data.frame(t(anno_sub[,-c(2:9)]))
colnames(anno_sub) <- anno_sub[1,]
anno_sub$Sample <- rownames(anno_sub)
anno_sub <- left_join(metadata,anno_sub, by="Sample")
colnames(phylo)[4] <- "ID"
phylo$ID <- as.factor(phylo$ID)
colnames(anno_sub)[1] <- "meta_ID"

# Subset nhp dataset
nhp <- anno_sub[anno_sub$Cohort %in% "NHP",]

anno_sub <- right_join(phylo,anno_sub[,c(1,2,3,6:ncol(anno_sub))], by="ID")
write.table(hits_sub, "full_results/sample_metabolite_matrix_raw_full.txt", sep="\t", quote=F, row.names = F)
write.table(anno_sub, "full_results/sample_meta.txt", sep="\t", quote=F, row.names = F)
#####


# SEPARATE BY DATA TYPE
#####
control <- anno_sub[anno_sub$Strain %in% "BHI-S media control",]
control1 <- control[,c(8:ncol(control))]
control1 <- sapply(control1, as.numeric)
control1 <- as.data.frame(t(colMeans(control1)))
control <- cbind(control[1,c(1:8)],control1)
control$Sample <- "BHI"

# B fragilis
bfrag <- anno_sub[anno_sub$Strain %in% "B. fragilis" | anno_sub$Strain %in% "BHI-S media control",]
bfrag[,c(19:ncol(bfrag))] <- sapply(bfrag[,c(19:ncol(bfrag))],as.numeric)
bfrag[,c(19:ncol(bfrag))] <- bfrag[,c(19:ncol(bfrag))] %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) < 1))
bfrag_samples <- bfrag$Sample
bfrag_matrix <- as.data.frame(t(bfrag[,c(19:ncol(bfrag))]))

#####


# PCA
#####
library(zCompositions)
library(vegan)
PCA <- as.data.frame(t(read.csv("resources/bfrag_quant.csv", row.names = 1)))
PCA <- na.omit(PCA[-c(1:12),])
PCA_clr <- PCA %>% decostand(method = "rclr")

# ALTERNATIVE OPTION FOR NORMALIZATION THE OTHER OPTION IS DECONSTAND ABOVE
# Remove columns that have 2 or less positive values (or else the 0 correction will not work)
PCA2 <- colSums(PCA==0)
PCA2 <- PCA2 < 253
PCA2 <- PCA2[PCA2 == FALSE]
PCA2 <- names(PCA2)
PCA <- PCA[,-c(as.numeric(PCA2))]
# Apply a correction for 0s
PCA2 <- cmultRepl(PCA, output = 'p-counts')
# Find the geometric mean of the row
gm <- apply(PCA2,1,function(x) exp(mean(log(x))))
PCA3 <- PCA2
# Center transform the data
for(i in 1:nrow(PCA3)){  PCA3[i,] <- PCA3[i,]/gm[i]}
# Take the natural log of the data
PCA4 <- log(PCA3)

# Add the metadata to the table
PCA_clr$Sample <- row.names(PCA_clr)
PCA_clr <- right_join(metadata,PCA_clr)
PCA_clr <- right_join(phylo,PCA_clr)

# Caculate the PCA of the data
#pca <- prcomp(PCA4[,c(11:ncol(PCA4))], center = TRUE, scale. = TRUE)
pca <- prcomp(PCA_clr[,c(11:ncol(PCA4))], center = TRUE, scale. = TRUE)

summary(pca)
# Plot the PCA per whatever metadata column
autoplot(pca, data=PCA_clr, colour="Media") + theme_classic()
ggsave("full_results/pca_media.png")
autoplot(pca, data=PCA_clr, colour="Assigned") + theme_classic()
ggsave("full_results/pca_assigned.png")

#####

