
#Heatmap plotting for Twins

rm(list=ls())
setwd("~/Desktop/twins/data")

#read in rarefied phyloseq object
physeq1 <- readRDS("physeq_rare_twins.rds")

#make the location_staph category
physeq1@sam_data$location_staph <- paste(physeq1@sam_data$Location, physeq1@sam_data$Staph.culture, sep="_")

#need to trim to 50 most abundant species, or OTUS?

top_50 <- prune_taxa(names(sort(taxa_sums(physeq1),TRUE)[1:50]), physeq1)
top_50_glom <- tax_glom(top_50, taxrank = 'Genus')

#plot the heatmap
plot_heatmap(top_50_glom, "PCoA", "bray", "Staph.culture", "Genus", weighted=TRUE) +
  facet_grid(~Location+Staph.culture, scales = "free_x")





