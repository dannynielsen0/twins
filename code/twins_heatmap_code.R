
#Heatmap plotting for Twins

rm(list=ls())
setwd("~/Desktop/twins/data")

#read in rarefied phyloseq object
physeq1 <- readRDS("physeq_rare_twins.rds")

#make the location_staph category
physeq1@sam_data$location_staph <- paste(physeq1@sam_data$Location, physeq1@sam_data$Staph.culture, sep="_")

#need to trim to 50 most abundant species, or OTUS?

top_50 <- prune_taxa(names(sort(taxa_sums(physeq1),TRUE)[1:50]), physeq1)
top_50_glom <- tax_glom(top_50, taxrank = 'Class')

#plot the heatmap
genus_heatmap <- plot_heatmap(top_50_glom, "PCoA", "bray", "Staph.culture", "Class", weighted=TRUE) +
  facet_grid(~Location+Staph.culture, scales = "free_x")


ggsave(plot=genus_heatmap, "../figures/class_level_heatmap.jpg", width = 14, height =6 , device='jpg', dpi=500)


