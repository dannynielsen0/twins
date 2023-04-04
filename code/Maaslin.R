rm(list=ls())
setwd("~/Desktop/twins/data")


#load libraries

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiomeMarker")

library(microbiomeMarker)
library(phyloseq)
library(lefser)
library(vegan)
library(Maaslin2)
library(ggpubr)
library(rstatix)
library(data.table)

#load example data
zeller_data <- zeller14@assays@data

###Load the Picrust data

crust_dat <- readRDS("crust_dat_LEFSER.rds")

#remove hand samples

throat_dat <- subset(crust_dat, crust_dat$Location=="Hand")
throat_dat <- crust_dat



###Maaslin 

meta <- throat_dat[,1:4] #set meta data
data <- throat_dat[,5:43] #set input data

#run Maaslin
fit_data = Maaslin2(
  input_data = data, 
  input_metadata = meta, 
  output = "Maaslin_Hand_byStaph",
  fixed_effects = c("Staph.culture"))
  #random_effects = c("Location"))
  #reference = "Location,Throat")

sig_res <- read.table("Maaslin_location/significant_results.tsv", sep="\t", header=TRUE)


#melt to long format the data with location and staph as the variables

crust_dat_long <- melt(data=crust_dat[,c(1,4:43)])

#change staph code 0,1 to meaningful name
levels(crust_dat_long$Staph.culture) <-c("Staph negative", "Staph positive")


#plot the levels

KEGG_2_location <- ggplot(data=crust_dat_long, aes(x=variable, y=sqrt(value), fill=Staph.culture)) +
  geom_boxplot(outlier.shape = NA, lwd=0.75) + coord_flip() + 
  theme_bw() + ylim(0,2000) +
  facet_wrap(~Location) +
  scale_fill_manual(values= c("white", "grey")) +
  theme(text=element_text(size=12,  family="sans")) + xlab("KEGG Level 2 Pathway") +
  ylab("Square-root Abundance") +
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
  guides(fill = guide_legend(reverse = TRUE)) 



ggsave(plot=KEGG_2_location, "../figures/KEGG_2_location.jpg", width = 12, height =8 , device='jpg', dpi=600)

  




