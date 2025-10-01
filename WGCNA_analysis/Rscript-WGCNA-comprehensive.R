

library(pheatmap)
library(WGCNA)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
set.seed(65)









wd= "/home/pgsb/sepideh.jafarian/bmw_project/WGCNA_comprehensive_traits_isotope"
setwd(wd)
getwd()




vsd1 <- read.csv("vsd.csv", sep=",", header=TRUE)
vsd1[, 2:ncol(vsd1)] <- lapply(vsd1[, 2:ncol(vsd1)], as.numeric) # as.numeric = typecast
#str(vsd1)
dim(vsd1)





dim(vsd1)
head(vsd1)






vsd<- t(vsd1)
head(vsd)

colnames(vsd)<- vsd[1,]
vsd<- vsd[-1, ]
any(is.na(vsd))



vsd<- as.data.frame(vsd)

class(vsd)
vsd[1]
vsd$SampleName <- rownames(vsd)
head(vsd)


# Reading the first meta table for joining to vsd


meta <- read.csv("Wheat_treatments_meta.csv", sep = ",", header = TRUE)
head(meta)
meta_pheno<- read.csv("meta_all_phenotypic_combined.csv", sep=",", header=TRUE)
head(meta_pheno)
colnames(meta_pheno)[24]<- "SampleName"




meta_pheno<- meta_pheno %>%
	  mutate(SampleName = paste0("BMW_", SampleName))



meta$batches<- as.factor(meta$batches)
meta$Genotype<- as.factor(meta$Genotype)
meta$Treatment<- as.factor(meta$Treatment)
meta$CID <- as.factor(meta$CID)
str(meta)







# Merging isotop variables to the meta table



c13_data<- read.csv("merged_pheno_manual_sum_iso_factors.csv", sep=",", header=TRUE)
head(c13_data)
colnames(c13_data)[4]<- "SampleName"

c13_data<- c13_data %>%
	  mutate(SampleName = paste0("BMW_", SampleName))

  meta_pheno<-inner_join(meta_pheno, c13_data , by="SampleName")
  meta_pheno<-meta_pheno %>% select(-SampleName, SampleName)







# joining vsd and meta to have factores included as well for the median later on

final_meta<- dplyr::inner_join(meta_pheno,meta, by = "SampleName")


vsd_final<- inner_join(vsd,final_meta, by= "SampleName")

vsd_final<- vsd_final %>% select(-SampleName, SampleName)



vsd_final[, 1:73722] <- lapply(vsd_final[, 1:73722], as.numeric)


vsd_final$batches<- as.factor(vsd_final$batches)
vsd_final$Genotype<- as.factor(vsd_final$Genotype)
vsd_final$Treatment<- as.factor(vsd_final$Treatment)
vsd_final$CID<- as.factor(vsd_final$CID)



vsd_final$batches = NULL
vsd_final$max.Height<- NULL
colnames(vsd_final)[colnames(vsd_final) == "Shoots"] <- "Nr_shoots"




## median calculation for groups
library(magrittr)
m_vsd<- vsd_final %>%
	  group_by(Genotype,Treatment, CID) %>%
	   dplyr::summarise(across(.cols = 1:73721, .fns = median, na.rm = TRUE), .groups = 'drop')

   head(m_vsd)


 # Renaming group levels

 m_vsd<- m_vsd %>%
         mutate(CID = ifelse(CID == "H", "High", "Low"))
head(m_vsd)


any_missing <- any(is.na(m_vsd))
print(any_missing)


m_vsd$SampleName <- paste(m_vsd$Genotype, m_vsd$Treatment,  m_vsd$CID, sep = "_")
head(m_vsd)


m_vsd <- as.data.frame(m_vsd)



#Reordering the sample names
m_vsd <- m_vsd %>%
	  mutate(SampleName = str_replace(SampleName, "^(.*)_(control|stress)_(High|Low)$", "\\3_\\2_\\1")) %>%
	    arrange(factor(str_sub(SampleName, 1, 2)), # Ordering by control_H, control_L, stress_H, stress_L
		              SampleName)
  
  print(m_vsd$SampleName)
    head(m_vsd)



m_vsd <- column_to_rownames(m_vsd, var = "SampleName")
head(m_vsd)

rownames(m_vsd)

head(m_vsd)


##dandrogram including bars for the traits

head(vsd_final)

vsd_final<- column_to_rownames(vsd_final, var="SampleName")

vsd_bar<- vsd_final[ , 73697:73721]
rownames(vsd_bar)
#vsd_bar<- column_to_rownames(vsd_bar, var="SampleName")
print(vsd_bar)




png(file=paste(wd, "samples dandrogram with bars.png", sep='/') ,width=20,height=10,units="in",res=1200)
par(mar=c(5, 50, 3, 0))
sampleTree<- hclust(dist(vsd_final), method = "average")
traitColors<-numbers2colors(vsd_bar, signed = FALSE);
plotDendroAndColors(sampleTree, traitColors,
		                        groupLabels = names(vsd_bar),
					                    main="sample dendrogram and heatmap")
dev.off()



head(m_vsd)
m_vsd<- m_vsd[, c(-1,-2, -3)]
head




# check for sample quality and remove remaining bad samples or genes
check <- goodSamplesGenes(m_vsd, verbose = 3)
if (!check$allOK) {
  # Print the gene and sample names that were removed:
  if (sum(!check$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(m_vsd)[!check$goodGenes], collapse = ", ")));
  if (sum(!check$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(m_vsd)[!check$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  m_vsd <- m_vsd[check$goodSamples, check$goodGenes]
}
print(paste0("Dimensions of final normalized count table:", dim(m_vsd)))
head(m_vsd)

# Find the soft-thresholding power
#### Based on Pearson-correlation

powers <- c(c(1:14), seq(from = 15, to = 21, by = 2))
sft_w <- pickSoftThreshold(m_vsd, powerVector = powers, verbose = 5, networkType = "unsigned")
# open display window for plotting
options(repr.plot.width = 8, repr.plot.height = 5, repr.plot.res = 200)
par(mfrow = c(1, 2))
cex1 = 0.9
# Plot the results
# Scale-free topology fit index as a function of the soft-thresholding power
png(file=paste(wd, "st_R square.png", sep='/') ,width=10,height=10,units="in",res=1200)
plot(sft_w$fitIndices[,1], -sign(sft_w$fitIndices[,3])*sft_w$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,unsigned R^2",type="n",
     main = paste("Scale independence"))
text(sft_w$fitIndices[,1], -sign(sft_w$fitIndices[,3])*sft_w$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.87,col="red")
dev.off()
# Mean connectivity as a function of the soft-thresholding power

png(file=paste(wd, "st_mean connectivity.png", sep='/') ,width=10,height=10,units="in",res=1200)
plot(sft_w$fitIndices[,1], sft_w$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft_w$fitIndices[,1], sft_w$fitIndices[,5], labels=powers, cex=cex1,col="red")


dev.off()

softPower <- sft_w$powerEstimate
print(paste0("Chosen softPower:", softPower))

# Calculate adjacency matrix and cluster by topological overlap
### calculate the adjacency matrix
### this is a metric to determine which genes j and i are connected


adjacency <- adjacency(m_vsd, power = softPower)


# calculate the topological overlap matrix (TOM)
# this estimates how many connections gene j and i have in common



TOM <- TOMsimilarity(adjacency)
colnames(TOM) <- colnames(m_vsd)
rownames(TOM) <- colnames(m_vsd)
collectGarbage() # clean R

# Module detection and merging


# Clustering
geneTree = hclust(as.dist(1-TOM), method = "average");
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = geneTree,
                               distM = 1-TOM,
                               deepSplit = 2,
                               cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)
table(unmergedColors)
table(unmergedLabels)







png(file=paste(wd, "Dynamic tree cut.png", sep='/') ,width=15,height=10,units="in",res=1200)
plotDendroAndColors(geneTree, unmergedColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()



# Module Eigengene calculation
MElist <- moduleEigengenes(m_vsd, colors = unmergedColors)
MEs <- MElist$eigengenes
# calculate dissimilarity of eigengenes
MEdiss <- 1-cor(MEs)
MEtree <- hclust(as.dist(MEdiss), method = "average")
head(MEs, 2)
dim(MEs)



# clustering of Module Eigengenes


png(file=paste(wd, "wheat-BMW-Clustering of module eigengenes.png", sep='/') ,width=15,height=10,units="in",res=1200)
options(repr.plot.width = 10, repr.plot.height = 8, repr.plot.res = 100)
par(mfrow = c(1,1))
plot(MEtree, main = "wheat-BMW-Clustering of module eigengenes", xlab = "", sub = "")
abline(h=0.3, col = "red")
abline(h=0.2, col = "green")
dev.off()



# What are exactly moduls eigengens values?
## ** they are representative gene expression profile of a cluster ** our goal is to know is there any difference in the gene expression profiles between groups!!

# Here we merge modules for which the eigengenes have a distance smaller than 0.1
# this roughly corresponds to a correlation of 0.7
MEdissThres = 0.3
merge <- mergeCloseModules(m_vsd, unmergedColors, cutHeight = MEdissThres)

# for comparison
merge <- mergeCloseModules(m_vsd, unmergedColors, cutHeight = MEdissThres)
mergedColors <- merge$colors
table(mergedColors)


merge_l <- mergeCloseModules(m_vsd, unmergedLabels, cutHeight = 
MEdissThres) 
merge_l<- merge_l$dendro 
mergedLabels <- merge_l$labels
length(mergedLabels)





# plot gene dendogram with new & old module association
png(file=paste(wd, "wheat-BMW-Clustering of module eigengenes_smaller03.png", sep='/') ,width=15,height=10,units="in",res=1200)
options(repr.plot.width = 8, repr.plot.height = 5, repr.plot.res = 200)
plotDendroAndColors(geneTree, cbind(unmergedColors, mergedColors),
                    c("DynamicTree cut", "MergedTree cut"),
                    dendroLabels = FALSE, hang = 0.03, addGuide = FALSE, guideHang = 0.05, abHeight = 0.9)
dev.off()

#Create heatmap displaying the relationship of the identified modules
# using the topological overlap dissimilarity matrix(TOM)
### lighter shades of yellow indicate closer expression neighborhood
##notice:command uses to much RAM



png(file=paste(wd, "heatmap for identified modules.png", sep='/') ,width=25,height=25,units="in",res=1200)
par(mar = c(15, 13, 3, 2.2))
disTOM <- 1-TOM
TOMplot(disTOM^7, geneTree, as.character(mergedColors))
dev.off()















# making modified meta table for the rest of the analysis and getting heat map


final_meta<- final_meta %>%
	  mutate(CID = ifelse(CID == "H", "High","Low"))
  head(final_meta)


head(final_meta)
final_meta$SampleName <- paste(final_meta$Genotype, final_meta$Treatment,  final_meta$CID, sep = "_")



meta2 <- final_meta %>%
	  mutate(SampleName = str_replace(SampleName, "^(.*)_(control|stress)_(High|Low)$", "\\3_\\2_\\1")) %>%
	    arrange(factor(str_sub(SampleName, 1, 2)), 
		              SampleName)

head(meta2)


meta2$max.Height<- NULL
colnames(meta2)[colnames(meta2) == "Shoots"] <- "Nr_shoots"

library(magrittr)
meta2<- meta2 %>%
	  group_by(Genotype,Treatment, CID, SampleName) %>%
	   dplyr::summarise(across(.cols = 1:25, .fns = median, na.rm = TRUE), .groups = 'drop')
   head(meta2)


   meta2 <- meta2 %>%
	     mutate(SampleName = str_replace(SampleName, "^(.*)_(control|stress)_(High|Low)$", "\\3_\\2_\\1")) %>%
	       arrange(factor(str_sub(SampleName, 1, 2)), 
		                 SampleName)

   print(meta2)





class(meta2)
colnames(meta2)
rownames(meta2)<- NULL
rownames(meta2)



meta2$Name<- meta2$SampleName
meta2<- column_to_rownames(meta2, var = "SampleName")
head(meta2)






#. Relating modules to trait information

Treatment<- meta2$Name
sample <- rownames(meta2)
bin_t = binarizeCategoricalVariable(Treatment, # same as binarizeCategoricalColumns()
                                    includePairwise = FALSE,
                                    includeLevelVsAll = TRUE,
                                    minCount = 1,
                                    dropFirstLevelVsAll = FALSE);
bin_t <- data.frame(sample, bin_t)
rownames(bin_t) <- bin_t[,1]
bin_t[,1] <- NULL
sorted_treatment <- unique(sort(Treatment))
colnames(bin_t) <- sorted_treatment
print(head(bin_t, 3))
print(dim(bin_t))

moduleTraitCor = cor(merge$newMEs, bin_t, use = "p");
write.table(moduleTraitCor, "moduleTraitCor_all_samples.csv", sep=",", row.names = TRUE, col.names = TRUE)
moduleTraitPvalue = corPvalueFisher(moduleTraitCor, 544)
write.table(moduleTraitPvalue, "moduleTraitPvalue_all_samples.csv", sep=",", row.names = TRUE, col.names = TRUE)
MEColors = substring(names(merge$newMEs), 3);
MEColorNames = paste(MEColors, sep="")
length(MEColorNames)




# Plot the module-trait relationship table as in WGCNA tutorial



png(file=paste(wd, "heatmap_for_all_samples.png", sep='/') ,width=18,height=23,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(15, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bin_t),
               yLabels = MEColorNames,
               xLabelsAngle = 90,
               ySymbols = MEColorNames,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
	       verticalSeparator.x = 8,
	       verticalSeparator.col = 1,
	       verticalSeparator.lty = 1,
               verticalSeparator.lwd = 4,
               verticalSeparator.ext = 0,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4, #0.3,
               cex.lab.x = 1.5,
               cex.lab.y = 1.2,
               cex.legendLabel = 1,
               legendLabel = "Pearson's correlation",
               zlim = c(-1,1),
               main = " Module Eigengene -genotypes relationships")
dev.off()








png(file=paste(wd, "heatmap for all samples2.png", sep='/') ,width=20,height=15,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
	       xLabels = names(bin_t),
               yLabels = MEColorNames,
	       xLabelsAngle = 45,
	       ySymbols = MEColorNames,
	       colorLabels = TRUE,
               colors = blueWhiteRed(50),
	       #verticalSeparator.x = 10,
	       #verticalSeparator.col = 1,
	       #verticalSeparator.lty = 1,
	       #verticalSeparator.lwd = 4,
	       #verticalSeparator.ext = 0,
	       textMatrix = textMatrix,
	       setStdMargins = FALSE,
	       cex.text = 0.4, #0.3,
	       cex.lab.x = 1.5,
	       cex.lab.y = 1.2,
	       cex.legendLabel = 1,
	       legendLabel = "Pearson's correlation",
	       zlim = c(-1,1),
	       main = " Module Eigengene -relationships")
dev.off()






png(file=paste(wd, "heatmap_for_all_samples3.png", sep='/') ,width=18,height=23,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(15, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bin_t),
               yLabels = MEColorNames,
               xLabelsAngle = 90,
               ySymbols = MEColorNames,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
	       #verticalSeparator.x = 8,
	       #verticalSeparator.col = 1,
	       #verticalSeparator.lty = 1,
               #verticalSeparator.lwd = 4,
               #verticalSeparator.ext = 0,
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4, #0.3,
               cex.lab.x = 1.5,
               cex.lab.y = 1.2,
               cex.legendLabel = 1,
               legendLabel = "Pearson's correlation",
               zlim = c(-1,1),
               main = " Module Eigengene -genotypes relationships")
dev.off()





## Heatmap between drought and control condition


Treatment1<- meta2$Treatment
sample <- rownames(meta2)
bin_t = binarizeCategoricalVariable(Treatment1, # same as binarizeCategoricalColumns()
                                    includePairwise = FALSE,
                                    includeLevelVsAll = TRUE,
                                    minCount = 1,
                                    dropFirstLevelVsAll = FALSE);
bin_t <- data.frame(sample, bin_t)
rownames(bin_t) <- bin_t[,1]
bin_t[,1] <- NULL
sorted_treatment1 <- unique(sort(Treatment1))
colnames(bin_t) <- sorted_treatment1
print(head(bin_t, 3))
print(dim(bin_t))

moduleTraitCor = cor(merge$newMEs, bin_t, use = "p");
write.table(moduleTraitCor, "moduleTraitCor_Drought_vs_control.csv", sep=",", row.names = TRUE, col.names = TRUE)
moduleTraitPvalue = corPvalueFisher(moduleTraitCor, 544)
write.table(moduleTraitPvalue, "moduleTraitPvalue_Drought_vs_control.csv", sep=",", row.names = TRUE, col.names = TRUE)
MEColors = substring(names(merge$newMEs), 3);
MEColorNames = paste(MEColors, sep="")
length(MEColorNames)

png(file=paste(wd, "heatmap for  Drought vs contro.png", sep='/') ,width=13,height=20,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bin_t),
               yLabels = MEColorNames,
               xLabelsAngle = 0,
               ySymbols = MEColorNames,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4, #0.3,
               cex.lab.x = 1.5,
               cex.lab.y = 1.2,
               cex.legendLabel = 1,
               legendLabel = "Pearson's correlation",
               zlim = c(-1,1),
               main = " Module Eigengene - Treatments relationship")
dev.off()






png(file=paste(wd, "heatmap_for_Drought_vs_contro.png", sep='/') ,width=13,height=20,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(bin_t),
               yLabels = MEColorNames,
               xLabelsAngle = 0,
               ySymbols = MEColorNames,
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4, #0.3,
               cex.lab.x = 1.5,
               cex.lab.y = 1.2,
               cex.legendLabel = 1,
               legendLabel = "Pearson's correlation",
               zlim = c(-1,1),
               main = " Module Eigengene - Treatments relationship")
dev.off()








## Heatmap between low and high CID groups


Treatment2<- meta2$CID
sample <- rownames(meta2)
bin_t = binarizeCategoricalVariable(Treatment2, # same as binarizeCategoricalColumns()
                                    includePairwise = FALSE,
                                    includeLevelVsAll = TRUE,
                                    minCount = 1,
                                    dropFirstLevelVsAll = FALSE);
bin_t <- data.frame(sample, bin_t)
rownames(bin_t) <- bin_t[,1]
bin_t[,1] <- NULL
sorted_treatment2 <- unique(sort(Treatment2))
colnames(bin_t) <- sorted_treatment2
print(head(bin_t, 3))
print(dim(bin_t))

moduleTraitCor = cor(merge$newMEs, bin_t, use = "p");
write.table(moduleTraitCor, "moduleTraitCor_h_vs_l.csv", sep=",", row.names = TRUE, col.names = TRUE)
moduleTraitPvalue = corPvalueFisher(moduleTraitCor, 544)
write.table(moduleTraitPvalue, "moduleTraitPvalue_h_vs_l.csv", sep=",", row.names = TRUE, col.names = TRUE)
MEColors = substring(names(merge$newMEs), 3);
MEColorNames = paste(MEColors, sep="")
length(MEColorNames)




png(file=paste(wd, "heatmap for HCID vs LCID.png", sep='/') ,width=13,height=18,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
	       xLabels = names(bin_t),
               yLabels = MEColorNames,
	       ySymbols = MEColorNames,
	       xLabelsAngle = 0,
	       colorLabels = TRUE,
	       colors = blueWhiteRed(50),
	       #textMatrix = textMatrix,
	       setStdMargins = FALSE,
               cex.text = 0.5, #0.3,
               cex.lab.x = 1,   
	       cex.lab.y = 1.2, 
	       cex.legendLabel = 1.5, 
	       legendLabel = "Pearson's correlation",
	       zlim = c(-1,1),
	       main = " Module Eigengene - CID groups relationship")
dev.off()





png(file=paste(wd, "heatmap_for_HCID_vs_LCID.png", sep='/') ,width=13,height=18,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
	       xLabels = names(bin_t),
               yLabels = MEColorNames,
	       ySymbols = MEColorNames,
	       xLabelsAngle = 0,
	       colorLabels = TRUE,
	       colors = blueWhiteRed(50),
	       textMatrix = textMatrix,
	       setStdMargins = FALSE,
               cex.text = 0.4, #0.3,
               cex.lab.x = 1,   
	       cex.lab.y = 1.2, 
	       cex.legendLabel = 1.5, 
	       legendLabel = "Pearson's correlation",
	       zlim = c(-1,1),
	       main = " Module Eigengene - CID groups relationship")
dev.off()





# modul trait relation for phenotyping traits

meta3<- meta2[,4:28]



moduleTraitCor = cor(merge$newMEs , meta3, use = "p");
write.table(moduleTraitCor, "moduleTraitCor_phenotypes.csv", sep=",", row.names = TRUE, col.names = TRUE)
moduleTraitPvalue = corPvalueFisher(moduleTraitCor, 544)
write.table(moduleTraitPvalue, "moduleTraitPvalue_phenotypes.csv", sep="," ,row.names = TRUE, col.names = TRUE)
MEColors = substring(names(merge$newMEs), 3);
MEColorNames = paste(MEColors, sep="")



png(file=paste(wd, "heatmap for phenotyping traits.png", sep='/') ,width=20,height=18,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
	       xLabels = names(meta3),
	       yLabels = MEColorNames,
	       xLabelsAngle = 45,
	       ySymbols = MEColorNames,
	       colorLabels = TRUE,
	       colors = blueWhiteRed(50),
	       textMatrix = textMatrix,
	       setStdMargins = FALSE,
	       cex.text = 0.4, #0.3,
	       cex.lab.x = 1.2,
	       cex.lab.y = 1.5,
	       cex.legendLabel = 1,
	       legendLabel = "Pearson's correlation",
	       zlim = c(-1,1),
	       main = " Module Eigengene -relationships")
dev.off()











png(file=paste(wd, "heatmap for phenotyping traits2.png", sep='/') ,width=20,height=18,units="in",res=1200)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 13, 3, 2.2)); #bottom, left, top, and right
labeledHeatmap(Matrix = moduleTraitCor,
	       xLabels = names(meta3),
	       yLabels = MEColorNames,
	       xLabelsAngle = 45,
	       ySymbols = MEColorNames,
	       colorLabels = TRUE,
	       colors = blueWhiteRed(50),
	       #textMatrix = textMatrix,
	       setStdMargins = FALSE,
	       cex.text = 0.4, #0.3,
	       cex.lab.x = 1.2,
	       cex.lab.y = 1.5,
	       cex.legendLabel = 1,
	       legendLabel = "Pearson's correlation",
	       zlim = c(-1,1),
	        main = " Module Eigengene -relationships")
dev.off()













# Hub genes detection

hubs= chooseTopHubInEachModule(m_vsd, mergedColors)
hubs
  hubs <- as.data.frame(hubs)
  hubs$modules<- rownames(hubs)
    colnames(hubs)
    print(hubs)
      rownames(hubs) <- NULL
    colnames(hubs) <- c("gene_id", "modules")
write.csv(hubs, "Hub_genes.csv", sep=",", row.names=FALSE)






##Extract Gene-Module information
genelist <- row.names(TOM)
gene_mod <- as.data.frame(cbind(genelist, mergedLabels, mergedColors))
colnames(gene_mod) <- c("ID","module_label","module_color")
str(gene_mod)
write.csv(gene_mod, "genes_modules_info.csv", sep=",") 




#calculate module sizes
module_sizes <- table(mergedColors)
print(module_sizes)
# Create a data frame with module colors and their sizes
module_size_df <- as.data.frame(module_sizes)
colnames(module_size_df) <- c("module_color", "size")
print(module_size_df)

# Sort the data frame based on sizes in descending order
module_size_df <- module_size_df[order(-module_size_df$size), ]
print(module_size_df)


write.table(module_size_df, "ordered_modules_sizes.csv", sep=",", row.names = FALSE)







