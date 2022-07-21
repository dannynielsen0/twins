###LEFSER analysis with the twin PICRUST2 data


rm(list=ls())
setwd("~/Desktop/twins/data")


#load libraries

library(phyloseq)
library(lefser)
library(vegan)

#load example data
zeller_data <- zeller14@assays@data



#load phyloseq object

physeq <- readRDS("physeq_twins.rds")

#convert staph culture to factor
physeq@sam_data$Staph.culture <- as.factor(physeq@sam_data$Staph.culture)
physeq@sam_data$Location <- as.factor(physeq@sam_data$Location)
physeq@sam_data$Description <- as.character(physeq@sam_data$Description)
physeq@sam_data$Twindex <- as.factor(physeq@sam_data$Twindex)



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

#create dataset for online galaxy LEFSE analysis
crust_dat$sample <- row.names(crust_dat)

#subset hand out
crust_dat <- subset(crust_dat, crust_dat$Location != "Hand")
#change staph to positive and negative
levels(crust_dat$Staph.culture) <- as.factor(c("negative", "positive"))
crust_dat$Location <- as.factor(crust_dat$Location)

#convert to percentage
crust_dat[,5:43] <- crust_dat[,5:43]/rowSums(crust_dat[,5:43])
crust_dat$sample <- as.factor(crust_dat$sample)

ggplot(data=crust_dat, aes(x=Staph.culture, y=crust_dat$Replication_and_Repair)) + geom_boxplot() + 
  stat_summary(fun = mean)

#make df
crust_dat_lefse <- data.frame(crust_dat[,c(1,44,5:43)])
row.names(crust_dat_lefse) <- NULL

write.table(crust_dat_lefse, "../data/crust_lefse.tsv", sep="\t", row.names = FALSE)


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

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# vegan::adonis2(formula = dist ~ Location + Staph.culture, data = crust_dat)
# Df SumOfSqs      R2       F Pr(>F)    
# Location        2   0.8481 0.06462 15.2420  0.001 ***
#   Staph.culture   1   0.0907 0.00691  3.2617  0.057 .  
# Residual      438  12.1856 0.92847                   
# Total         441  13.1244 1.00000                   
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1



#LEFSER analysis

#only 2 contrasts at a time,

#for some reason, I can only get it to run if I remove nose here

throat_dat <- crust_dat

#look at table for # of samples fitting each Staph condition in the throat, nose, hand
table(crust_dat$Location,crust_dat$Staph.culture)


#set up data for the summarized experiment
throat_counts <- data.frame(throat_dat[,5:43])
throat_colData <- throat_dat[,1:4]

#construct the summarized experiment with throat data

throat_exp <- SummarizedExperiment(assays=list(counts=t(throat_counts)),
                                   colData=throat_colData)

#remove hand

throat_exp <- throat_exp[,throat_exp$Location!="Nose"]
table(throat_exp$Location,throat_exp$Staph.culture)

#run lefser for throat data
throat_res <- lefser(throat_exp, groupCol = "Staph.culture")

#plot results
throat_plot <- lefserPlot(throat_res)
throat_plot


