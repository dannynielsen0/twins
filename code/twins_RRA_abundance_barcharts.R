##RRA barchart plotting for Twins

rm(list=ls())
setwd("~/Desktop/twins/data")

library(phyloseq)
library(RColorBrewer)


#read in rarefied phyloseq object
physeq1 <- readRDS("physeq_rare_twins.rds")

#number of OTUs after rarefying
physeq1 #1347
#after rarefying, # of reads
sum(sample_sums(physeq1@otu_table))


physeq1@sam_data$location_staph <- paste(physeq1@sam_data$Location, physeq1@sam_data$Staph.culture, sep="_")


#organize data and plot RRA barcharts

########To plot the desired taxanomic level, change each level between this line and the next string of #s below

y1 <- tax_glom(physeq1, taxrank = 'Phylum') # agglomerate taxa at Genus level
y2 = merge_samples(y1, "location_staph") # merge samples on sample variable of interest
y3 <- transform_sample_counts(y2, function(x) x/sum(x)) #get abundance in %
y4 <- psmelt(y3) # create dataframe from phyloseq object
y4$Phylum <- as.character(y4$Phylum) #convert to character
y4$Phylum[y4$Abundance < 0.05] <- "Phyla < 5% abund." #rename genera with < .25% abundance

#split back out the location and staph into unique columns
y4 <- cbind(y4, read.table(text=y4$Sample, sep="_", header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
y4$Phylum <- factor(y4$Phylum, levels=rev(unique(y4$Phylum)))


#set color palette to accommodate the number of genera

mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(length(unique(y4$Phylum)))


#plot
p <- ggplot(data=y4, aes(x=col2, y=Abundance, fill=Phylum, Phylum =Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=mycolors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="bottom", legend.title=element_text(size=12), legend.text=element_text(size=12)) + guides(fill=guide_legend(nrow=5)) + 
  scale_x_discrete(limits=rev(levels(as.factor(y4$col2)))) + facet_wrap(~col1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  ylab("16S relative read abundance") + xlab("") +
  theme(strip.text = element_text(size=12, family = "sans")) +
  theme(text=element_text(size=12,family="sans"))

p


#family
ggsave(plot=p, "Phylum_RRA.jpg", width = 14, height =8 , device='jpg', dpi=500)
########


#This is just a quick plotting of data to dummy check my above resulting figures,
#change names and taxonomic levels as needed

clostridium <- subset_taxa(physeq1, Genus=="Pseudomonas")

plot_bar(clostridium) + facet_wrap(~location_staph)



#get top 50 otus per each location

#subset for each location, and glom to genus level
physeq_throat <- subset_samples(physeq1, Location =="Throat")
physeq_throat <- tax_glom(physeq_throat, taxrank = "Genus")

physeq_nose <- subset_samples(physeq1, Location =="Nose")
physeq_nose <- tax_glom(physeq_nose, taxrank = "Genus")

physeq_hand <- subset_samples(physeq1, Location =="Hand")
physeq_hand <- tax_glom(physeq_hand, taxrank = "Genus")


#list top 50 of each
throat_50 <- prune_taxa(names(sort(taxa_sums(physeq_throat), TRUE)) [1:50], physeq_throat)
nose_50 <- prune_taxa(names(sort(taxa_sums(physeq_nose), TRUE)) [1:50], physeq_nose)
hand_50 <- prune_taxa(names(sort(taxa_sums(physeq_hand), TRUE)) [1:50], physeq_hand)

write.table(physeq_throat@tax_table, "top50_throat.txt", sep="\t")
write.table(physeq_nose@tax_table, "top50_nose.txt", sep="\t")
write.table(physeq_hand@tax_table, "top50_hand.txt", sep="\t")



