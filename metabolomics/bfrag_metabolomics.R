#################
# Title: bfrag_metabolomics.R
# Author- Renee Oles
# Purpose: 
# Output:  
################

# Load libraries
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(plotly)
library(Hotelling)
library(reshape2)
library(stats)
library(ggforce)
#library(KEGGgraph)
library(KEGGREST)
library(zCompositions)
library(vegan)

# Requirement: run metabolomics_step1.R to get bfrag matrix and dataframe

# NORMALIZE THE DATA # 
bfrag1 <- cbind(bfrag[,c(1:18)],bfrag[,c(19:ncol(bfrag))] %>% decostand(method = "rclr"))

# PCA PLOT # 
#####
pca_clr <- bfrag %>% 
  select_if(negate(function(col) is.numeric(col) && length(unique(col))==1))
pca <- prcomp(pca_clr[,c(14:ncol(pca_clr))], center = TRUE, scale. = TRUE)
summary(pca)
# Plot the PCA per whatever metadata column
pca_plot <- cbind(pca_clr[,c(1:12)],pca[["x"]])
ggplot(pca_plot, aes(x=PC1, y=PC2,colour=Division)) +
  geom_text(label=pca_plot$Sample, size = 3)+
  theme_classic()
#####

# CORRECT FOR CONTROLS # 
#####
#Now we need to divide it up to subtract the correct control from each of the samples
# subtract the hmo control from each of the samples
control <- bfrag1[!(bfrag1$Assigned %in% "Hiu") & bfrag1$Strain %in% "BHI-S media control",]
control_rest <- colMeans(control[,c(19:ncol(control))])
outliers <- c("Sample37", "Sample72", "Sample97", "Sample213", "Sample108", "Sample62")
bfrag_rest <- bfrag1[!(bfrag1$Assigned %in% "Hiu") & !(bfrag1$Strain %in% "BHI-S media control"),]
bfrag_rest <- bfrag_rest[!(bfrag_rest$meta_ID %in% outliers),]
bfrag_rest <- bfrag_rest[!duplicated(bfrag_rest$Sample),]
bfrag_minus_ctrl <- sweep(bfrag_rest[c(19:ncol(bfrag_rest))], 2, control_rest)
bfrag_minus_ctrl <- cbind(bfrag_rest[,c(1:19)], bfrag_minus_ctrl)
d1 <- bfrag_rest[bfrag_rest$Division %in% 1,]
d2 <- bfrag_rest[bfrag_rest$Division %in% 2,]
d1_sum <- colMeans(d1[,c(19:ncol(d1))])
d2_sum <- colMeans(d2[,c(19:ncol(d2))])
bfrag1 <- rbind(control, d1)
bfrag_rest <- d1_sum
#####

# PCA WITH CONTROL CORRECTION # 
#####
# PCA plot after correcting for media control
pca_clr <- bfrag_minus_ctrl[,c(19:ncol(bfrag_minus_ctrl))] %>% 
  select_if(negate(function(col) is.numeric(col) && length(unique(col))==1))
pca <- prcomp(pca_clr, center = TRUE, scale. = TRUE)
summary(pca)
# Plot the PCA per whatever metadata column
pca_plot <- cbind(bfrag_minus_ctrl[,c(1:18)],pca[["x"]])
ggplot(pca_plot, aes(x=PC1, y=PC2,colour=as.factor(Group))) +
  #geom_text(label=pca_plot$Sample, size = 3)+
  labs(y= "PC2 13.1%", x = "PC1 15.8%") + 
#  labs(y= "PC3 11.9%", x = "PC2 13.1%") + 
#  labs(y= "PC3 11.9%", x = "PC4 11.2%") + 
  geom_point(size=2.5)+
  scale_color_manual(values=c("Healthy"="#0E131F", "Infection" = "#62688D", "IBD" = "#3C6E71", "Unknown" = "#BBBFC2"))+
  theme_classic()
ggsave("bfrag/pca_division_group.png", dpi=300, height=3.2, width=5, units="in")
#####

# METABOLOMICS MATRIX # 
#####
# create a matrix from a dataframe that has each sample as the column and each of the metabolites as a row
metabolomics_matrix <- as.data.frame(t(bfrag_minus_ctrl[,c(19:ncol(bfrag_minus_ctrl))]))
colnames(metabolomics_matrix) <- bfrag_minus_ctrl$Sample
metabolomics <-  melt(bfrag_minus_ctrl[,c(1,2,6,7,19:ncol(bfrag_minus_ctrl))], id.vars = c("Sample","Division", "Group","Site"))
colnames(metabolomics)[5] <- "Scan"
metabolomics$Scan <- as.numeric(as.character(metabolomics$Scan))
metabolomics <- left_join(metabolomics, hits_sub)
write.table(metabolomics_matrix, "bfrag/bfrag_metabolomics_matrix_norm_full.txt", quote=F, sep="\t")
write.table(t(metabolomics_matrix), "bfrag/bfrag_metabolomics_matrix_norm_full_t.txt", quote=F, sep="\t")
#####

# PLOTS # 
#####
# See the distribution of the metabolomics data
# Violin_plot
ggplot(metabolomics, aes(x=ID,y=value))+
  geom_violin()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90))+
  ylab(label="ln Intensity")

# calculate the spread between all the metabolites
spread <- as.data.frame(t(apply(metabolomics_matrix,1,function(x) summary(x))))
spread1 <- as.data.frame((apply(metabolomics_matrix,1,function(x) sd(x))))
spread <- cbind(rownames(metabolomics_matrix), spread, spread1)
colnames(spread)[1] <- "Scan"
colnames(spread)[8] <- "SD"
spread$Scan <- as.numeric(spread$Scan)
metabolomics <- left_join(metabolomics,spread)
metabolomics$q_spread <- abs(metabolomics$`1st Qu.`- metabolomics$`3rd Qu.`)
metabolomics$r_spread <- abs(metabolomics$Min.- metabolomics$Max.)

metabolomics$Scan <- as.character(metabolomics$Scan)
metabolomics <- metabolomics %>%  mutate(Compound_Name = coalesce(Compound_Name,Scan))

# Plot annotated metabolites
ggplot(metabolomics[!is.na(metabolomics$Scan) & metabolomics$q_spread <= 2 & metabolomics$r_spread >= 5,], aes(x=Compound_Name,y=value))+
  geom_boxplot()+
  theme_classic()+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, size=9))+
  #  facet_grid(~subclass, scales = "free",space="free")+
  ylab(label="normalized abundance")
ggsave("bfrag/divergent_boxplot.png",dpi=300)

classes <- metabolomics[!duplicated(metabolomics$class),c(8:14)]
write.table(classes, "classes.txt", row.names=F, sep="\t",quote=F)
colors <- read.delim("class_color.txt")
metabolomics <- left_join(metabolomics, colors)
col <- as.character(metabolomics$color)
names(col) <- as.character(metabolomics$class)

#known_metabolites <- metabolomics[!is.na(metabolomics$subclass) & metabolomics$Min. < 0 & metabolomics$Max. > 0 & metabolomics$SD >= 0.0 & metabolomics$SD < 0.2,]
known_metabolites <- metabolomics[!is.na(metabolomics$class) & metabolomics$class != "N/A",]
#known_metabolites <- metabolomics[metabolomics$superclass %in% "Phenylpropanoids and polyketides",]
#known_metabolites <- metabolomics[!is.na(metabolomics$KEGG),]

core <- known_metabolites[!is.na(known_metabolites$Scan) & known_metabolites$SD <= 0.25 & known_metabolites$Mean != 0,]
col_current <- col[names(col) %in% core$class]
ggplot(core, aes(x=Compound_Name,y=value, fill=class, color = class))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(size=0.3, color="black")+
  theme_classic()+
  theme(axis.text = element_text(size=12))+
  scale_fill_manual(values=col_current)+
  scale_color_manual(values=col_current)+
  coord_flip()+
  theme(axis.text.x = element_text(angle = 90, size=9))+
  #  facet_grid(~subclass, scales = "free",space="free")+
  ylab(label="Media corrected normalized abundance") 
ggsave("violin_plots/known_metabolites_sdle025_box.png",dpi=300, width = 14, height=8, units="in")

variable <- known_metabolites[known_metabolites$SD >= 0.75,]
col_current <- col[names(col) %in% variable$class]
ggplot(variable, aes(x=Compound_Name,y=value, fill = class, color = class))+
  geom_boxplot(outlier.size = 0.5)+
  geom_point(size=0.3, color="black")+
  theme_classic()+
  theme(axis.text = element_text(size=12))+
#  facet_wrap(~PATH, scales = "free")+
  scale_fill_manual(values=col_current)+
  scale_color_manual(values=col_current)+
  guides(fill=guide_legend(ncol =1))+
  coord_flip()+
  ylab(label="Media corrected normalized abundance") 
# key is g = greater than l = less than and e = equal
ggsave("violin_plots/known_metabolites_sdge075_box.png",dpi=300, width = 14, height=10, units="in")
#ggsave("known_metabolites_sdge01_sdl02_maxg0_minl0.png",dpi=300, width = 15, height=20, units="in")
#####

# HEATMAP #
#####
# Make a heatmap of annotated metabolites 
library(pheatmap)
library(scales)
metabolomics_heatmap <- metabolomics[metabolomics$q_spread >= 1 & metabolomics$r_spread >= 3,]
metabolomics_heatmap <- reshape(metabolomics_heatmap[,c(1,7,6)], timevar = "Sample", idvar = "Compound_Name", direction = "wide")
rownames(metabolomics_heatmap) <- metabolomics_heatmap[,1]
colnames(metabolomics_heatmap) <- gsub("value.","",colnames(metabolomics_heatmap))

# Set the class color 
class_color <- data.frame("Compound_Name"= metabolomics_heatmap[,c(1)])
class_color <- left_join(class_color, unique(metabolomics[,c(7,10,26)]))
class_color[is.na(class_color$class),2] <- "N/A"
class_color[class_color$class == "N/A",3] <- "white"
class_color_sub <- unique(class_color[,c(2,3)])
class_color_list <- class_color_sub$color
names(class_color_list) <- unique(class_color_sub$class)
classdf <- data.frame(row.names = class_color$Compound_Name, 
                      category = class_color$class)
# Format heatmap
metabolomics_heatmap <- metabolomics_heatmap[,-1]

# Set the sample color
sample_color <- data.frame("Sample"= colnames(metabolomics_heatmap))
sample_color <- left_join(sample_color, bfrag[,c(1,6)])
sample_color <- unique(sample_color)
sample_color[is.na(sample_color$Group),2] <- "Unknown"
sample_color_sub <- unique(sample_color[,c(2)])
sample_color <- data.frame(row.names = sample_color$Sample, 
                      group = sample_color$Group)

main_list <- list(category=class_color_list, group = c("Healthy"="#0E131F", "Infection" = "#4A4E69", "IBD" = "#3C6E71", "Unknown" = "#E7EBEF"))

# Set the color sclae
scale <- colorRampPalette(c("#440154", "white", "#2A788E"))
# Save the plot
png(filename="heatmap.png",width=12, height=17.5, units="in", res=300)
pheatmap(metabolomics_heatmap, 
        column_title = "Strains", row_title = "Metabolites", annotation_row = classdf, annotation_col = sample_color, annotation_colors = main_list, show_colnames = F, show_rownames = F, color = scale(30))
dev.off()
#####


# DIFFERENTIAL ABUNDANCE TESTING #
#####
calc_ttest <- function(full_matrix, gr1, gr2, maxAdjP, minLog2FC, c1, c2) {
  #  df <- df[df$group == gr1 | df$group == gr2,]
  df_ttest <- metabolomics_summary[,c(gr1,gr2)]
  rownames(df_ttest) <- rownames(metabolomics_summary)
  #Log2 fold change group2 - group1
  df_ttest$Log2FC <- df_ttest[,2] - df_ttest[,1]
  df_ttest$pval <- apply(full_matrix, 1, function(x) t.test(x[c(c1)], x[c(c2)])$p.value)
  #Benjamini-Hochberg correction for multiple testing
  df_ttest$adjPval <- p.adjust(df_ttest$pval, method = "BH")
  df_ttest$Log10adjPval <- -1*log10(df_ttest$adjPval)
  #Add the categorical column for easier visualization
  df_ttest <- df_ttest[df_ttest$pval != "NaN",]
  df_ttest$Diff_Abund <- apply(
    df_ttest, 1, function(x) {
      if(x[["pval"]] != "NaN"){
        if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] >= minLog2FC) {
          return( paste(gr2, "Produced") )
        } else if (x[["adjPval"]] <= maxAdjP & x[["Log2FC"]] <= -1*minLog2FC) {
          return( paste(gr2, "Depleted") )
        } else {
          return('Non-significant')
        }
      }
    }
  )
  df_ttest
}

# Find differential metabolites
gr1 <- "BHI-S media control"
gr2 <- "B. fragilis"
rownames(bfrag) = seq(length=nrow(bfrag))
c1 <- which(grepl("BHI-S media control", bfrag$Strain))
c2 <- which(grepl("B. fragilis", bfrag$Strain))
minLog2FC = 0
maxAdjP = 0.05
metabolomics_summary <- data.frame("BHI-S media control"=control_rest, "B. fragilis"=bfrag_rest, check.names = FALSE)
bfrag_matrix <- t(bfrag[,-c(1:18)])
df_ttest <- calc_ttest(bfrag_matrix, gr1, gr2, maxAdjP, minLog2FC, c1, c2)
df_ttest$Scan <- rownames(df_ttest)
#write.table(df_ttest,paste(name,"tp1_vs_tp2.txt",sep=""),quote=FALSE,sep="\t")
df_ttest$Scan <- as.numeric(as.character(df_ttest$Scan))
df_ttest <- left_join(df_ttest, hits_sub)
df_ttest$Scan <- as.factor(df_ttest$Scan)
df_ttest <- df_ttest %>%  mutate(Compound_Name = coalesce(Compound_Name,Scan))
df_ttest_kegg <- df_ttest[!is.na(df_ttest$KEGG),c(9,3,10,11)]
write.table(df_ttest_kegg, "df_ttest_kegg_compounds.txt", quote=F, sep="\t")
write.table(df_ttest, "df_ttest_bfrag_v_media.txt", quote=F, sep="\t")


# Volcano Plot
volc <- ggplot(df_ttest[!duplicated(df_ttest$Scan),],aes(x = Log2FC, y = Log10adjPval, colour = as.factor(Diff_Abund))) +
  geom_point(shape=19, size=1, alpha = 0.6)+
  geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65", linetype = 2) +
#  geom_vline(xintercept = 0, colour = "gray65") +
  #geom_vline(xintercept = -1*minLog2FC, colour = "gray65", linetype = 2) +
  geom_vline(xintercept = minLog2FC, colour = "gray65", linetype = 2) +
  ggtitle(paste("T-test ", gr1, " vs ", gr2," Adjusted P-value<=", maxAdjP, " Log2 FC>=", minLog2FC,
                sep="")) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=10),
        plot.title = element_text(size=10)) +
  scale_color_manual(values=c("#440154","#2A788E","grey"))+
#  geom_text(data = subset(df_ttest, (Log2FC >=minLog2FC | Log2FC <= -minLog2FC) & adjPval < 0.05),
#            aes( Log2FC, Log10adjPval, label = Compound_Name),
#            alpha = 0.6, hjust = 0.5, vjust = -0.6, size=2)+
  labs(x = paste("Log2 FC", gr2, "-", gr1), y = "-Log10 Adj. P-value" )
volc
ggsave("volcano_plot.png", dpi=300, height=5, width=5, units="in")

# Plot of metabolite groups
df_ttest[is.na(df_ttest$subclass),17] <- "N/A"
df_ttest[is.na(df_ttest$class),12] <- "N/A"
df_ttest[is.na(df_ttest$superclass),12] <- "N/A"
ggplot(df_ttest[!duplicated(df_ttest$Scan),], aes(x = class, fill= Diff_Abund))+
  geom_bar(position="fill")+
  geom_text(
    aes(label=..count..),
    color = "white",
    stat="count",
    position=position_fill(vjust=0.5)) +
  scale_fill_manual(values=c("#440154","#2A788E","grey"))+
  coord_flip()+
  theme_classic()+
  theme(axis.text=element_text(size=15))
ggsave("bfrag/prop_plot_class.png", dpi=300, height=10, width=9, units="in")



# Find differential metabolites between D1 and D2
gr1 <- "1"
gr2 <- "2"
c1 <- which(grepl(1, bfrag$Division))
c2 <- which(grepl(2, bfrag$Division))
minLog2FC = 0
maxAdjP = 0.05
metabolomics_summary <- data.frame("1" = d1, "2" = d2, check.names = FALSE)
matrix_bfrag_minus_ctrl <- as.data.frame(t(bfrag_minus_ctrl[,c(14:ncol(bfrag_minus_ctrl))]))
# for division test sub t.test for wilcox.test
df_ttest <- calc_ttest(matrix_bfrag_minus_ctrl, gr1, gr2, maxAdjP, minLog2FC, c1, c2)
df_ttest$Scan <- rownames(df_ttest)
#write.table(df_ttest,paste(name,"tp1_vs_tp2.txt",sep=""),quote=FALSE,sep="\t")
df_ttest$Scan <- as.numeric(as.character(df_ttest$Scan))
df_ttest <- left_join(df_ttest, anno[,c(1:3)])
df_ttest$Scan <- as.factor(df_ttest$Scan)
df_ttest <- df_ttest %>%  mutate(Compound_Name = coalesce(Compound_Name,Scan))

# Volcano Plot
volc <- ggplot(df_ttest,aes(x = Log2FC, y = Log10adjPval, colour = Diff_Abund )) +
  geom_point(shape=19, size=1, alpha = 0.6)+
  geom_hline(yintercept = -1*log10(maxAdjP), colour = "gray65", linetype = 2) +
  #  geom_vline(xintercept = 0, colour = "gray65") +
  geom_vline(xintercept = -1*minLog2FC, colour = "gray65", linetype = 2) +
  geom_vline(xintercept = minLog2FC, colour = "gray65", linetype = 2) +
  ggtitle(paste("T-test ", gr1, " vs ", gr2," Adjusted P-value<=", maxAdjP, " Log2 FC>=", minLog2FC,
                sep="")) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.text = element_text(size=8),
        plot.title = element_text(size=10)) +
  scale_color_manual(values=c("#78a2f0","#f09c78","grey"))+
  #  geom_text(data = subset(df_ttest, (Log2FC >=minLog2FC | Log2FC <= -minLog2FC) & adjPval < 0.05),
  #            aes( Log2FC, Log10adjPval, label = Compound_Name),
  #            alpha = 0.6, hjust = 0.5, vjust = -0.6, size=2)+
  labs(x = paste("Log2 FC", gr2, "-", gr1), y = "-Log10 Adj. P-value" )
volc
ggsave("volcano_plot.png", dpi=300, height=5, width=6, units="in")

# Plot of metabolite groups
ggplot(df_ttest[!is.na(df_ttest$subclass),], aes(x = subclass, fill= Diff_Abund))+
  geom_bar(position="fill")+
  scale_fill_manual(values=c("#78a2f0","#f09c78","grey"))+
  coord_flip()+
  theme_classic()
ggsave("prop_plot.png", dpi=300, height=10, width=6, units="in")


# Kegg annotation
df_ttest$Scan <- as.numeric(df_ttest$Scan)
df_ttest_kegg <- left_join(df_ttest, hits_sub, by="Scan")
df_ttest_kegg <- df_ttest_kegg[!is.na(df_ttest_kegg$KEGG),]
df_ttest_kegg <- df_ttest_kegg[,c()]
#####

