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
  fixed_effects = c("Staph.culture")),
  random_effects = c("Location"))
  #reference = "Location,Throat")

sig_res <- read.table("Maaslin_location/significant_results.tsv", sep="\t", header=TRUE)


#melt to long format the data with location and staph as the variables

crust_dat_long <- melt(data=crust_dat[,c(1,4:43)])

#change staph code 0,1 to meaningful name
levels(crust_dat_long$Staph.culture) <-c("Staph negative", "Staph positive")

Aggregate_crust_location_long_mean <- stats::aggregate(value ~ Location + variable, crust_dat_long, function(x) mean=mean(x))
Aggregate_crust_location_long_sd <- stats::aggregate(value ~ Location + variable, crust_dat_long, function(x) sd=sd(x))

crust_dat_long <- cbind(Aggregate_crust_location_long_mean, Aggregate_crust_location_long_sd$value)

colnames(crust_dat_long)[3] <- "mean"
colnames(crust_dat_long)[4] <- "sd"

crust_dat_long$mean <- sqrt(crust_dat_long$mean)
crust_dat_long$sd<- sqrt(crust_dat_long$sd)

# annotation table with adjusted pvals and y-position of the labels
anno_df = compare_means(value ~ Staph.culture, group.by = c("Location", "variable"), data = crust_dat_long, method="t.test") %>%
  mutate(y_pos = 2000)

ggplot(crust_dat_long, aes(x=variable, y=value, fill=Staph.culture)) + 
  geom_boxplot(position=position_dodge()) + 
  #geom_point(aes(color=supp), position=position_jitterdodge()) + 
  facet_wrap(~Location) + coord_flip() +
  ggsignif::geom_signif(
    data=anno_df, 
    aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
    manual=TRUE
  )


ggplot(df, aes(x=supp, y=len)) + 
  geom_boxplot(position=position_dodge()) + 
  geom_point(aes(color=supp), position=position_jitterdodge()) + 
  facet_wrap(~dose) + 
  ggsignif::geom_signif(
    data=anno_df, 
    aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
    manual=TRUE
  )


# Add 10% spaces between the p-value labels and the plot border
KEGG_2_location + stat_pvalue_manual(
  stat.test,  label = "p", tip.length = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

KEGG_2_location <- ggplot(data=crust_dat_long, aes(x=variable, y=sqrt(value), fill=Staph.culture)) +
  geom_boxplot(outlier.shape = NA, lwd=0.75) + coord_flip() + 
  theme_bw() + ylim(0,3000) +
  facet_wrap(~Location) + 
  scale_fill_manual(values= c("white", "grey")) +
  theme(text=element_text(size=12,  family="sans")) + xlab("KEGG Level 2 Pathway") +
  ylab("Abundance") +
  stat_summary(fun=mean, geom="point", shape=20, size=1, color="red",
               position = position_dodge2(width = 0.75,   
                                          preserve = "single")) +
  stat_compare_means(label="p.signif", method = "kruskal", 
                     hide.ns = TRUE, label.y = 1500, size=7, vjust = -.225) +
  guides(fill = guide_legend(reverse = TRUE)) 


KEGG_2_location <- ggplot(data=crust_dat_long, aes(x=variable, y=sqrt(value), fill=Staph.culture)) +
  geom_boxplot() + coord_flip() + 
  geom_signif(test="wilcox.test",
              comparisons = list("Staph positive", "Staph negative")) +
  geom_col(aes(fill=Location), position=position_dodge(0.85), stat = "identity") + theme_bw() +
  coord_flip() + #facet_wrap(~Staph.culture) + 
  geom_errorbar(aes(ymin = mean, ymax=mean+sd), position = position_dodge(.9)) +
  scale_fill_manual(values= c("forestgreen", "maroon", "blue")) +
  theme(text=element_text(size=12,  family="sans")) + xlab("KEGG Level 2 Pathway") +
  ylab("Abundance") + 
  guides(fill = guide_legend(reverse = TRUE)) 

ggsave(plot=KEGG_2_location, "../figures/KEGG_2_location.jpg", width = 12, height =8 , device='jpg', dpi=600)

  
KEGG_2_staph <- ggplot(data=crust_dat_long, aes(x=variable, y=sqrt(value), group=Staph.culture)) +
  geom_col(aes(fill=Staph.culture), position=position_dodge(0.85), stat = "identity") + theme_bw() +
  coord_flip() + facet_wrap(~variable) + 
  scale_fill_manual(values= c("grey","black")) +
  theme(text=element_text(size=12,  family="sans")) + xlab("KEGG Level 2 Pathway") +
  ylab("Abundance") + 
  guides(fill = guide_legend(reverse = TRUE))

ggsave(plot=KEGG_2_staph, "../figures/KEGG_2_staph.jpg", width = 12, height =8 , device='jpg', dpi=600)






