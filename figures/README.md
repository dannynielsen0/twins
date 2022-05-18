
rm(list=ls())
setwd("~/Desktop/twins/data")


#libraries
library(qiime2R)
library(phyloseq)
library(vegan)
library(data.table)
library(ggplot2)

physeq <- qza_to_phyloseq(features="deblur_table_final.qza", tree = "asvs-tree.qza",
                taxonomy = "classification.qza", metadata = "metadata2.txt")

#convert staph culture to factor
physeq@sam_data$Staph.culture <- as.factor(physeq@sam_data$Staph.culture)


#plot sampling depth and rarefaction 

#sampling depth
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
summary(sdt)



##Rarefy
set.seed(666)
physeq_rare <- rarefy_even_depth(physeq,sample.size = min(sample_sums(physeq)))


###Alpha diversity by location and between staph carriage
plot_richness(physeq, x="Location", measures="Chao1") + facet_grid(~Staph.culture)


#beta diversity
#convert to RRA

physeq_RRA <- transform_sample_counts(physeq, function(x) x/sum(x)) #converts to relative abundance

#create ordinations
PCoA_bray <- ordinate(physeq_rare, method = "PCoA", distance = "bray", scale=TRUE, center=TRUE) #ordination using bray-curtis distances
PCoA_unifrac <- ordinate(physeq_rare, method = "PCoA", distance = "wunifrac") #ordination using bray-curtis distances

#plot ordinations
PCoA_plot1 <- plot_ordination(physeq_rare, PCoA_unifrac, color = "Location", shape="Staph.culture", axes = 1:2)

PCoA_plot1 <- PCoA_plot1 + geom_point(size = 3) + 
  scale_color_manual(values= c("forestgreen", "maroon", "blue")) +
  scale_shape_manual(values=c(15,3)) + theme_bw() + 
  ggtitle("Microbiome Composition") + #scale_x_reverse() + #scale_y_reverse() +
  theme(plot.title = element_text(size = 18, face= "bold")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title=element_text(size=16), 
        legend.text=element_text(size=14)) +
  theme(plot.title = element_text(size = 18, face= "bold"))

PCoA_plot1


###PERMANOVA
#create metadata for permanova
perm_meta <- as(physeq_RRA@sam_data, "data.frame")

#make bray distance matrices for 16S composition
perm_dist_bray <- phyloseq::distance(physeq_RRA, method="bray", scale=TRUE, center=TRUE)
perm_dist_unifrac <- phyloseq::distance(physeq_RRA, method="wunifrac")

#run adonis for PERMANOVA
bray_perm <- adonis(data=perm_meta, perm_dist_unifrac~Location + Staph.culture)

bray_perm

