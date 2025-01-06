### Examples from a collaboration on a recent project studying the effects of staph infection on microbiome composition in a human twin study. For complete project code, see https://github.com/calacademy-research/Dalman_Twin_Microbiome_Scripts. Paper titled, Staphylococcus aureus carriage is associated with microbiome composition in the nares and oropharynx, not the hand, of monozygotic twins, currently in press at Frontiers in Microbiomes - doi: 10.3389/frmbi.2024.1457940.

## Below is a portion of code used to summarize microbiome composition at three sites (hand, throat, and nose) in indivuals with and without prior staph infection.
```
####Below shows the relative read abundance at the genus level. Change the taxonmic rank to any desirable level.

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
```
## Plot using ggplot

```
p <- ggplot(data=y4, aes(x=col2, y=Abundance, fill=Phylum, Phylum =Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values=mycolors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position="bottom", legend.title=element_text(size=12), legend.text=element_text(size=12)) + guides(fill=guide_legend(nrow=5)) + 
  scale_x_discrete(limits=rev(levels(as.factor(y4$col2)))) + facet_wrap(~col1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  ylab("16S relative read abundance") + xlab("") +
  theme(strip.text = element_text(size=12, family = "sans")) +
  theme(text=element_text(size=12,family="sans"))
```

![Relative read abundance of microbial phyla across sites and between staph postiive and negative samples](https://github.com/dannynielsen0/twins/blob/main/figures/Phylum_RRA.jpg)



