
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

physeq <- qza_to_phyloseq(features="deblur_table_final.qza", tree = "asvs-tree.qza",
                          taxonomy = "classification.qza", metadata = "metadata2.txt")

saveRDS(physeq, "physeq_twins.rds")


#convert staph culture to factor
physeq@sam_data$Staph.culture <- as.factor(physeq@sam_data$Staph.culture)
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

KO_crust <- read.csv("KEGG_cat2_data.csv", header=TRUE, stringsAsFactors = FALSE) #for the csv
#KO_crust <- read.table("KEGG_cat2.tsv", sep="\t", header = FALSE, row.names = 1) #for txt
KO_crust$X <- gsub(" ", "_", KO_crust$X) #add _ to the spaces between parts of names here

#set up df
crust_descript <- t(KO_crust) #transform so that samples are rows

#make first line the column header
colnames(crust_descript) <- crust_descript[1,]#only if reading in the csv
crust_descript <- crust_descript[-1, ] #only for csv

#cbind the picrust data to the physeq metadata
crust_dat <- cbind(physeq@sam_data, crust_descript)
crust_dat[,5:43] <- as.numeric(unlist(crust_dat[,5:43])) #make path abundances numeric


#plot PCA of functional pathways

pca <- prcomp(sqrt(crust_dat[,5:43]), center=T, scale=T)

# visualize
biplot <- fviz_pca_biplot(pca, repel = TRUE, axes = c(1,2),
                          select.var = list(contrib = 10), #draw top 10 arrows
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
                          title = "KEGG Level 2 Pathways")
biplot

#PERMANOVA with adonis

dist <- vegdist(sqrt(crust_dat[,5:43]), method="bray")
vegan::adonis2(dist ~ Location + Staph.culture, data=crust_dat)

# Permutation: free
# Number of permutations: 999
# 
# vegan::adonis2(formula = dist ~ Location + Staph.culture, data = crust_dat)
# Df SumOfSqs      R2       F Pr(>F)    
# Location        2   0.8481 0.06462 15.2420  0.001 ***
#   Staph.culture   1   0.0907 0.00691  3.2617  0.056 .  
# Residual      438  12.1856 0.92847                   
# Total         441  13.1244 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#LEFSER analysis

#only 2 contrasts at a time,
#let's look at staph pos and neg in each location separately
#subset datasets for each location
throat_dat <- subset(crust_dat, Location == "Throat")
throat_dat$Location <- levels(droplevels(throat_dat$Location))
throat_dat$Location <- factor(throat_dat$Location)

nose_dat <- subset(crust_dat, Location == "Nose")
hand_dat <- subset(crust_dat, Location == "Hand")

#look at table for # of samples fitting each Staph condition in the throat, nose, hand
table(throat_dat$Location, throat_dat$Staph.culture)
table(nose_dat$Location, nose_dat$Staph.culture)
table(hand_dat$Location, hand_dat$Staph.culture)

#set up data for the summarized experiment
throat_counts <- data.frame(throat_dat[,5:43])
throat_colData <- throat_dat[,1:4]

#construct the summarized experiment with throat data

throat_exp <- SummarizedExperiment(assays=list(counts=t(throat_counts)),
                     colData=throat_colData)

#run lefser for throat data
throat_res <- lefser(throat_exp, groupCol = "Staph.culture", blockCol = NULL)



#plot results
throat_plot <- lefserPlot(throat_res)
throat_plot

