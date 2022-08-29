
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

#For alpha diversity plotting, make location_staph variable
physeq_rare@sam_data$loc_staph <- as.factor(paste(physeq_rare@sam_data$Location, physeq_rare@sam_data$Staph.culture, sep = "_"))

###Alpha diversity by location and between staph carriage (*** 0.001, ** 0.01, * 0.05)

div_df <- data.frame(physeq_rare@sam_data)

my_comparisons <- list(c("Hand_staph-negative","Nose_staph-negative"),
                        c("Hand_staph-negative","Nose_staph-positive"),
                        c("Hand_staph-negative","Throat_staph-negative"),
                        c("Hand_staph-negative","Throat_staph-positive"),
                       c("Hand_staph-positive","Nose_staph-negative"),
                       c("Hand_staph-positive","Nose_staph-positive"),
                       c("Nose_staph-negative","Throat_staph-negative"),
                       c("Nose_staph-negative","Throat_staph-positive"),
                       c("Nose_staph-positive","Throat_staph-negative"),
                       c("Nose_staph-positive","Throat_staph-positive"))
                       
                       
#shannon diversity
alpha_div_plot_shannon <- plot_richness(physeq_rare, x="loc_staph", measures="Shannon") + geom_boxplot() + 
  theme_bw() + ylab("Shannon Diversity") + xlab("") +
  stat_compare_means(comparisons = my_comparisons,
                     label="p.signif", method="wilcox.test", hide.ns = TRUE, vjust = 0.5,
                     label.y = c(3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3,5.5,5.7)) +
  #stat_compare_means(comparisons = combn(levels(physeq_rare@sam_data$loc_staph),2, simplify = F)[-15],
                    #label="p.signif", method="wilcox.test", hide.ns = TRUE) +
  scale_y_continuous(breaks=seq(1,4,by=1)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size=12, family = "sans")) +
  theme(text=element_text(size=12,family="sans"))


my_comparisons_chao1 <- list(c("Hand_staph-negative","Nose_staph-negative"),
                       c("Hand_staph-negative","Nose_staph-positive"),
                       c("Hand_staph-negative","Throat_staph-negative"),
                       c("Hand_staph-negative","Throat_staph-positive"),
                       c("Hand_staph-positive","Nose_staph-negative"),
                       c("Hand_staph-positive","Nose_staph-positive"),
                       c("Hand_staph-positive","Throat_staph-negative"),
                       c("Hand_staph-positive","Throat_staph-positive"),
                       c("Nose_staph-negative","Throat_staph-negative"),
                       c("Nose_staph-negative","Throat_staph-positive"),
                       c("Nose_staph-positive","Throat_staph-negative"),
                       c("Nose_staph-positive","Throat_staph-positive"),
                       c("Throat_staph-negative","Throat_staph-positive"))


#Chao1 diversity
alpha_div_plot_chao1 <- plot_richness(physeq_rare, x="loc_staph", measures="Chao1") + geom_boxplot() + 
  theme_bw() + ylab("Chao1 Diversity") + xlab("") +
  stat_compare_means(comparisons = my_comparisons_chao1,
                     label="p.signif", method="wilcox.test", hide.ns = TRUE, vjust = 0.5,
                     label.y = c(275,295,315,335,355,375,395,415,435,455,475, 495, 515)) +
  scale_y_continuous(breaks=seq(0,300,by=100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size=12, family = "sans")) +
  theme(text=element_text(size=12,family="sans"))


my_comparisons_simpson <- list(c("Hand_staph-negative","Nose_staph-negative"),
                             c("Hand_staph-negative","Nose_staph-positive"),
                             c("Hand_staph-positive","Nose_staph-negative"),
                             c("Hand_staph-positive","Nose_staph-positive"),
                             c("Nose_staph-negative","Throat_staph-negative"),
                             c("Nose_staph-negative","Throat_staph-positive"),
                             c("Nose_staph-positive","Throat_staph-negative"),
                             c("Nose_staph-positive","Throat_staph-positive"))



#Simpson Diversity
alpha_div_plot_simpson <- plot_richness(physeq_rare, x="loc_staph", measures="Simpson") + geom_boxplot() + 
  theme_bw() + ylab("Simpson Diversity") + xlab("") +
  stat_compare_means(comparisons = my_comparisons_simpson,
                     label="p.signif", method="wilcox.test", hide.ns = TRUE, vjust = 0.5,
                     label.y = c(1,1.1,1.2,1.3,1.4,1.5,1.6,1.7)) +
  scale_y_continuous(breaks=seq(0,1.5,by=.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size=12, family = "sans")) +
  theme(text=element_text(size=12,family="sans"))
  

ggsave(plot=alpha_div_plot_shannon, "../figures/shannon_div.jpg", width = 10, height =8, device='jpg', dpi=500)
ggsave(plot=alpha_div_plot_chao1, "../figures/chao1_div.jpg", width = 10, height =8, device='jpg', dpi=500)
ggsave(plot=alpha_div_plot_simpson, "../figures/simpson_div.jpg", width = 10, height =8, device='jpg', dpi=500)


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



