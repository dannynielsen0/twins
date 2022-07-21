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

#load example data
zeller_data <- zeller14@assays@data

###Load the Picrust data

crust_dat <- readRDS("crust_dat_LEFSER.rds")

#remove hand samples

throat_dat <- subset(crust_dat, crust_dat$Location=="Nose")
throat_dat <- crust_dat



###Maaslin 

meta <- throat_dat[,1:4] #set meta data
data <- throat_dat[,5:43] #set input data

#run Maaslin
fit_data = Maaslin2(
  input_data = data, 
  input_metadata = meta, 
  output = "Maaslin_location",
  fixed_effects = c("Staph.culture"),
  random_effects = c("Location"))
  #reference = "Location,Throat")

sig_res <- read.table("Maaslin_location/significant_results.tsv", sep="\t", header=TRUE)




