###DESEQ2 analysis with the twin data


rm(list=ls())
setwd("~/Desktop/twins/data")


#load libraries

library(phyloseq)
library(lefser)
library(vegan)
library(DESeq2)
library(ggplot2)


#load physeq object
physeq <- readRDS("physeq_twins.rds")


#make a taxonomic rank that includes both family and genus - this will help when visualizing later
Fam_Genus <- paste(tax_table(physeq)[ ,"Family"],tax_table(physeq)[ ,"Genus"], sep = "_")
tax_table(physeq) <- cbind(tax_table(physeq), Fam_Genus)


#subset to only throat
physeq_throat <- subset_samples(physeq, physeq@sam_data$Location != "Hand")

#set no staph as reference level for diet treatment
physeq_throat@sam_data$Staph.culture = relevel(physeq_throat@sam_data$Location, "Throat")

#create the deseq object and run the function  
diff_abund = phyloseq_to_deseq2(physeq_throat, ~Location)
diff_abund = DESeq(diff_abund, test="Wald", fitType ="parametric", sfType = "poscounts") #the poscounts sfType flag can be used if many 0s in the data

#check the reference levels
resultsNames(diff_abund)



#set contrast to compare each so that comparing staph against no staph
res = results(diff_abund, contrast = list(c("Location_Throat_vs_Nose"))) # neg is greater in nose; positive is greater in throat
sigtab = res[which(res$padj < 0.05), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))



#ggplot settings
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab$Phylum = factor(as.character(sigtab$Phylum), levels=names(x))
# Genus order
x = tapply(sigtab$log2FoldChange, sigtab$Fam_Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab$Fam_Genus = factor(as.character(sigtab$Fam_Genus), levels=names(x))



### We can plot the results  
#positive values mean greater abundance in staph positive, negative values mean greater in staph negative

diff_abund_plot <- ggplot(sigtab, aes(x=Fam_Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + coord_flip()

print(diff_abund_plot)




clostridium <- subset_taxa(physeq, Genus=="Bergeyella")

plot_bar(clostridium) + facet_wrap(~Location)






