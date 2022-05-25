
rm(list=ls())
setwd("~/Desktop/twins/data")


#libraries
library(qiime2R)
library(phyloseq)
library(vegan)
library(data.table)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(factoextra)

physeq <- qza_to_phyloseq(features="deblur_table_final.qza", tree = "asvs-tree.qza",
                          taxonomy = "classification.qza", metadata = "metadata2.txt")


#convert staph culture to factor
physeq@sam_data$Staph.culture <- as.factor(physeq@sam_data$Staph.culture)


#plot sampling depth and rarefaction 

#sampling depth
sdt = data.table(as(sample_data(physeq), "data.frame"),
                 TotalReads = sample_sums(physeq), keep.rownames = TRUE)
summary(sdt)
summary(sample_sums(physeq_rare))



##Rarefy
set.seed(666)
physeq_rare <- rarefy_even_depth(physeq,sample.size = min(sample_sums(physeq)))
summary(sample_sums(physeq_rare)) #using above rarefy command, we have this resulting: 2527



###Alpha diversity by location and between staph carriage
alpha_div_plot_shannon <- plot_richness(physeq_rare, x="Location", measures="Shannon") + facet_grid(~Staph.culture)
alpha_div_plot_chao1 <- plot_richness(physeq_rare, x="Location", measures= "Chao1") + facet_grid(~Staph.culture)

ggsave(plot=alpha_div_plot_shannon, "../figures/shannon_div.pdf", width = 10, height =10 , device='pdf', dpi=500)
ggsave(plot=alpha_div_plot_chao1, "../figures/chao1_div.pdf", width = 10, height =10 , device='pdf', dpi=500)


#beta diversity
#convert to RRA, if wanted

physeq_RRA <- transform_sample_counts(physeq, function(x) x/sum(x)) #converts to relative abundance

#create ordinations
PCoA_bray <- ordinate(physeq_rare, method = "PCoA", distance = "bray", scale=TRUE, center=TRUE) #ordination using bray-curtis distances
PCoA_unifrac <- ordinate(physeq_rare, method = "PCoA", distance = "unifrac") #ordination using bray-curtis distances

#plot ordinations
PCoA_plot1 <- plot_ordination(physeq_rare, PCoA_bray, color = "Location", shape="Staph.culture", axes = 1:2)

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

ggsave(plot=PCoA_plot1, "../figures/U_Unifrac_pcoa.pdf", width = 10, height =10 , device='pdf', dpi=500)



###PERMANOVA
#create metadata for permanova
perm_meta <- as(physeq_rare@sam_data, "data.frame")

#make bray distance matrices for 16S composition
perm_dist_bray <- phyloseq::distance(physeq_rare, method="bray", scale=TRUE, center=TRUE)
perm_dist_unifrac <- phyloseq::distance(physeq_RRA, method="wunifrac")

#run adonis for PERMANOVA
bray_perm <- adonis(data=perm_meta, perm_dist_bray~Location + Staph.culture)



###Picrust 

#Use KO pathways
#change the file to .csv and read in this way, the other format was creating some issues with weird downstream problems 
#probably due to weird characters throwing off the df

KO_crust <- read.csv("../picrust/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat_descrip.csv", header=TRUE, stringsAsFactors = FALSE)


#set up df
KO_crust <- KO_crust[,c(2:444)] #remove first column
crust_descript <- t(KO_crust) #transform so that samples are rows

#make first line the column header
colnames(crust_descript) <- crust_descript[1,]
crust_descript <- crust_descript[-1, ] 

#cbind the picrust data to the physeq metadata
crust_dat <- cbind(physeq@sam_data, crust_descript)
crust_dat[,5:6303] <- as.numeric(unlist(crust_dat[,5:6303])) #make path abundances numeric


#plot PCA of functional pathways

pca <- prcomp(sqrt(crust_dat[,5:6303]), center=T, scale=T)

# visualize
biplot <- fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                          select.var = list(contrib = 25), #draw top 25 arrows
                          #select.var = list(name = c("Q375E", "Q375P")),  #alternative to draw specific substitution loadings
                          addEllipses = TRUE,
                          habillage = crust_dat$Location,
                          col.ind = crust_dat$Staph.culture,
                          palette = c("forestgreen", "maroon", "blue"),
                          ellipse.level=0.95,
                          geom=c("point"), pointsize = 3.5,   #change to geom=c("point","text") for sample ID
                          ind.shape = crust_dat$Location,
                          ind.fill = crust_dat$Staph.culture,
                          invisible = c( "quali"), #remove enlarged symbol for group mean
                          title = "Picrust KO Pathways")
biplot

#PERMANOVA with adonis

dist <- vegdist(crust_dat[,5:6303], method="bray")
vegan::adonis(dist ~ Location + Staph.culture, data=crust_dat)

# Permutation: free
# Number of permutations: 999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Location        2     8.856  4.4278  38.913 0.14997  0.001 ***
#   Staph.culture   1     0.356  0.3557   3.126 0.00602  0.017 *  
#   Residuals     438    49.839  0.1138         0.84401           
# Total         441    59.050                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1




