##2D PCoA plotting for Twins

library(phyloseq)
library(ggplot2)

rm(list=ls())
setwd("~/Desktop/twins/data")

#read in rarefied phyloseq object
physeq_rare <- readRDS("physeq_rare_twins.rds")


#create ordinations for plotting below
PCoA_bray <- ordinate(physeq_rare, method = "PCoA", distance = "bray", scale=TRUE, center=TRUE) #ordination using bray-curtis distances
PCoA_unifrac <- ordinate(physeq_rare, method = "PCoA", distance = "unifrac") #ordination using unifrac distances
PCoA_wunifrac <- ordinate(physeq_rare, method = "PCoA", distance = "wunifrac") #ordination using weighted unifrac distances


#build plots by replacing the desired ordination (From last step) into the below plotting code

PCoA_plot <- plot_ordination(physeq_rare, PCoA_wunifrac, shape= "Staph.culture", color = "Location", axes = c(2,3))

plot <- PCoA_plot + geom_point(size = 5) +
  scale_color_manual(values= c("forestgreen", "maroon", "blue")) +
  scale_shape_manual(values=c(1,16)) + theme_bw() + ggtitle("Weighted unifrac") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14))

print(plot)


ggsave(plot=plot, "../figures/Wunifrac_2_3.pdf", width = 12, height =10 , device='pdf', dpi=500)





