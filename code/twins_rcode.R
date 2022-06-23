
rm(list=ls())
setwd("~/Desktop/twins/data")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("phyloseq")
library(phyloseq)


install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
library(qiime2R)

BiocManager::install("lefser")
BiocManager::install("DESeq2")


#libraries
library(qiime2R)
library(phyloseq)
library(vegan)
library(data.table)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(factoextra)
library(lefser)
library(DESeq2)
library(SummarizedExperiment)
library(ggpubr)
library(plotly)


physeq <- qza_to_phyloseq(features="deblur_table_final.qza", tree = "asvs-tree.qza",
                          taxonomy = "classification.qza", metadata = "metadata2.txt")

saveRDS(physeq, "physeq_twins.rds")
physeq <- readRDS("physeq_twins.rds")



#rename staph to Staph positive and Staph negative
levels(physeq@sam_data$Staph.culture) <- as.factor(c("staph-negative", "staph-positive"))
#change other variable to factor/character
physeq@sam_data$Location <- as.factor(physeq@sam_data$Location)
physeq@sam_data$Description <- as.character(physeq@sam_data$Description)

#plot sampling depth and rarefaction 

#sampling depth
sdt = data.table(as(physeq@sam_data, "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
summary(sdt) 
# TotalReads   
# Min.   : 2527  
# 1st Qu.: 5362  
# Median : 8530  
# Mean   :10589  
# 3rd Qu.:13815  
# Max.   :46685  
summary(sample_sums(physeq_rare))



##Rarefy
set.seed(666)
physeq_rare <- rarefy_even_depth(physeq,sample.size = min(sample_sums(physeq)))
summary(sample_sums(physeq_rare)) #using above rarefy command, we have this resulting: 2527

#save rarefied object
saveRDS(physeq_rare, "physeq_rare_twins.rds")


###Alpha diversity by location and between staph carriage (*** 0.001, ** 0.01, * 0.05)

#shannon diversity
alpha_div_plot_shannon <- plot_richness(physeq_rare, x="Location", measures="Shannon") + geom_boxplot() + 
  facet_grid(~Staph.culture) + theme_bw() + ylab("Shannon Diversity") +
  geom_signif(test="wilcox.test",  map_signif_level = TRUE,
              comparisons = combn(levels(physeq_rare@sam_data$Location),2, simplify = F)[-2],
              step_increase = 0.2)

#Chao1 diversity
alpha_div_plot_chao1 <- plot_richness(physeq_rare, x="Location", measures="Chao1") + geom_boxplot() + 
  facet_grid(~Staph.culture) + theme_bw() + ylab("Chao1 Diversity") +
  geom_signif(test="wilcox.test",  map_signif_level = TRUE,
              comparisons = combn(levels(physeq_rare@sam_data$Location),2, simplify = F)[-2],
              step_increase = 0.2)

#Simpson Diversity
alpha_div_plot_simpson <- plot_richness(physeq_rare, x="Location", measures="Simpson") + geom_boxplot() + 
  facet_grid(~Staph.culture) + theme_bw() + ylab("Simpson Diversity") +
  geom_signif(test="wilcox.test",  map_signif_level = TRUE,
              comparisons = combn(levels(physeq_rare@sam_data$Location),2, simplify = F)[-2],
              step_increase = 0.2)


ggsave(plot=alpha_div_plot_shannon, "../figures/shannon_div.pdf", width = 10, height =10 , device='pdf', dpi=500)
ggsave(plot=alpha_div_plot_chao1, "../figures/chao1_div.pdf", width = 10, height =10 , device='pdf', dpi=500)
ggsave(plot=alpha_div_plot_simpson, "../figures/simpson_div.pdf", width = 10, height =10 , device='pdf', dpi=500)


#beta diversity
#convert to RRA, if wanted
physeq_RRA <- transform_sample_counts(physeq, function(x) x/sum(x)) #converts to relative abundance

#create ordinations
PCoA_bray <- ordinate(physeq_rare, method = "PCoA", distance = "bray", scale=TRUE, center=TRUE) #ordination using bray-curtis distances
PCoA_unifrac <- ordinate(physeq_rare, method = "PCoA", distance = "unifrac") #ordination using unifrac distances
PCoA_wunifrac <- ordinate(physeq_rare, method = "PCoA", distance = "wunifrac") #ordination using weighted unifrac distances


###Plot Ordinations###

#bray
df_bray <- plot_ordination(physeq_rare, PCoA_bray, axes = c(1:3), justDF = T)
#get % variation explained for axes 1:3
PCoA_bray$values$Relative_eig[1:3]
# 0.21384710 0.10168689 0.07773238
 
#plot 3d plot with plot_ly

#bray
bray_3d <- plot_ly(
  df_bray, 
  x = ~Axis.1, 
  y = ~Axis.2, 
  z = ~Axis.3, 
  color = ~Location,
  colors = c("forestgreen", "maroon", "blue"),
  symbol = ~Staph.culture,
  symbols=c(1,15),
  mode = "markers",
  type = "scatter3d")
bray_3d <- bray_3d %>% layout(title="Bray-Curtis", scene = list(xaxis = list(title = 'Axis 1-21.4%'),
                                 yaxis = list(title = 'Axis 2-10.2%'),
                                 zaxis = list(title = 'Axis 3-7.8%')))


#unifrac
df_unifrac <- plot_ordination(physeq_rare, PCoA_unifrac, axes = c(1:3), justDF = T)
#get % variation explained for axes 1:3
PCoA_unifrac$values$Relative_eig[1:3]
# 0.14867343 0.14310856 0.05954748

#plot 3d plot with plot_ly
unifrac_3d <- plot_ly(
  df_unifrac, 
  x = ~Axis.1, 
  y = ~Axis.2, 
  z = ~Axis.3, 
  color = ~Location,
  colors = c("forestgreen", "maroon", "blue"),
  symbol = ~Staph.culture,
  symbols=c(1,15),
  mode = "markers",
  type = "scatter3d")
unifrac_3d <- unifrac_3d %>% layout(title="Unifrac", scene = list(xaxis = list(title = 'Axis 1-14.9%'),
                                           yaxis = list(title = 'Axis 2-14.3%'),
                                           zaxis = list(title = 'Axis 3-6.0%')))


#weighted unifrac
df_wunifrac <- plot_ordination(physeq_rare, PCoA_wunifrac, axes = c(1:3), justDF = T)
#get % variation explained for axes 1:3
PCoA_wunifrac$values$Relative_eig[1:3]
# 0.5395968 0.4217380 0.1327447

#plot 3d plot with plot_ly
wunifrac_3d <- plot_ly(
  df_wunifrac, 
  x = ~Axis.1, 
  y = ~Axis.2, 
  z = ~Axis.3, 
  color = ~Location,
  colors = c("forestgreen", "maroon", "blue"),
  symbol = ~Staph.culture,
  symbols=c(1,15),
  mode = "markers",
  type = "scatter3d")
wunifrac_3d <- wunifrac_3d %>% layout(title="Weighted Unifrac", scene = list(xaxis = list(title = 'Axis 1-54.0%'),
                                                                  yaxis = list(title = 'Axis 2-42.2%'),
                                                                  zaxis = list(title = 'Axis 3-13.3%')))



###PERMANOVA
#create metadata for permanova
perm_meta <- as(physeq_rare@sam_data, "data.frame")

#make bray distance matrices for 16S composition
perm_dist_bray <- phyloseq::distance(physeq_rare, method="bray")
perm_dist_unifrac <- phyloseq::distance(physeq_rare, method="unifrac")
perm_dist_wunifrac <- phyloseq::distance(physeq_rare, method="wunifrac")


#run adonis for PERMANOVA
bray_perm <- adonis2(perm_dist_bray~Location + Staph.culture, data=perm_meta)
write.table(bray_perm, "../figures/bray_permanova.txt", sep= "\t")

unifrac_perm <- adonis2(perm_dist_unifrac~Location + Staph.culture, data=perm_meta)
write.table(unifrac_perm, "../figures/unifrac_permanova.txt", sep= "\t")

wunifrac_perm <- adonis2(perm_dist_wunifrac~Location + Staph.culture, data=perm_meta)
write.table(wunifrac_perm, "../figures/wunifrac_permanova.txt", sep= "\t")



