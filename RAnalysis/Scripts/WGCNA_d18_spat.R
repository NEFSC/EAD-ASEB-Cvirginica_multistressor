---
# title: "Day18_WGCNA_all"
# author: "Samuel Gurr"
# date: "January 8, 2021"
---
  
# LOAD PACKAGES
library(WGCNA) # note: this was previously installed with the command `BiocManager::install("WGCNA")`
library(dplyr)
library(zoo)
library(DESeq2)

# for heatmap 
# library(devtools)
# install_github("jokergoo/ComplexHeatmap") first run these - commented out to avoid running a second time...
library(ComplexHeatmap)
library(circlize)
library(reshape)
library(ggplot2)
library(hrbrthemes)


# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")
# load in the count matrix with rownames as sample ID and colnames as gene ID

d18.filtered_count_tble  <- read.csv(file="Data/TagSeq/Filtered_counts/filtered_counts_5cpm_50perc/day18.filtered_5cpm50perc.csv", sep=',', header=TRUE) 

# d18.filtered_count_data          <- read.csv(file ="Presentations/URI/2022_URI_Puritz_Genomics_class/Day18_countmatrix_WGCNA.csv", sep =',', header=TRUE)
# Master.Treatment_Phenotype.data <- data.frame(read.csv(file="Presentations/URI/2022_URI_Puritz_Genomics_class/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE))

# The following setting is important, do not omit. (as recommended by WGCNA authors - view tutorial)
options(stringsAsFactors = FALSE)

# ===================================================================================
#
#
# Day 18 WGCNA - PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/14181-2105-9-559.pdf
# ===================================================================================

# count data  ========================================================== #

d18.data_matrix <- data.frame(d18.filtered_count_tble[,-1], row.names=d18.filtered_count_tble[,1]) 
dim(d18.data_matrix) # 4936   11 - 11 samples and 4820 total genes 


# trait data ========================================================== #

# Treatment data

d18.Treatment.data <- read.csv(file="Data/TagSeq/day18.exp.data.csv", sep=',', header=TRUE) %>%  
  dplyr::mutate_if(is.character, as.factor) %>% 
  dplyr::rename('Sample.Name' = 'SapleName_readmatrix') %>% 
  dplyr::rename('pCO2' = 'OA') %>% 
  dplyr::select(c('Sample.Name','Temperature','pCO2','Salinity')) %>% 
  dplyr::mutate(All_treatment = paste( (substr(Temperature,1,1)), 
                                       (substr(pCO2,1,1)), 
                                       (substr(Salinity,1,1)), sep = '')) %>% 
  dplyr::mutate(pCO2_Salinity = substr(All_treatment, 2,3))

View(d18.Treatment.data) # notice that there is ONE temperature for day 18


# sanity check =========== # 

dim(d18.Treatment.data)[1] ==  dim(d18.data_matrix)[2]# TRUE - each contains all 36 samples sequenced for Day 18 of the experiment 


# ===================================================================================
#
# Day 18 DESeqDataSet or 'dds' object (using DESeq2) 
#       vst transform the 'dds' onject for WGCNA 
# 
# ===================================================================================

# create dds  ========================================================== #

# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an object to transform, we do not need a design here...
dds.d18 <- DESeqDataSetFromMatrix(countData = d18.data_matrix,
                                 colData = d18.Treatment.data, design = ~ 1) # DESeq Data Set (dds)
dds.d18 # view the DESeqDataSet - notice the colData containg our critical treatment and sample ID data, rownames, etc. 


# transform the data  ========================================================== #
# run in order (kept name as dds.d18_vst)
dds.d18_vst <- vst(dds.d18) # transform it vst (variance stabilized transformation)
dds.d18_vst <- assay(dds.d18_vst) # call only the transformed coutns in the dds object
#fix(dds.d18_vst)
dds.d18_vst <- t(dds.d18_vst) # transpose columns to rows and vice versa

# ===================================================================================
#
# Day 18 WGCNA Sample tree - omit outlier sample(s)
#
# ===================================================================================

# checks before we start....
dim(dds.d18_vst) #  4936 genes; 11  samples
gsg = goodSamplesGenes(dds.d18_vst, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 

# determine outlier(s) with a sample tree    ========================================================== #

# call the cluster and set window dimenstions to view..
sampleTree = hclust(dist(dds.d18_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))
# output the tree as .png file..
png("Output/WGCNA/day18_spat/day18_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG59; can remove by hand or an automatic appraoch 
abline(h = 35, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

# view your tree   ========================================================== #
# NOTE: I inserted a horizontal line at 110 demonstrating where we likely want to cut our tree - omitting a single outlier 

# cut the tree and omit  ========================================================== #
clust = cutreeStatic(sampleTree, cutHeight = 35, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 1 and save 35; exactly what we want!  
keepSamples = (clust==1) # 'keepsamples' boolean to call the main dataset - notice there are TWO occurrences of FALSE - these are C6.larva and B12.larva

# integrate keepsamples  ========================================================== #
dds.d18_vst = dds.d18_vst[keepSamples, ] # integrate the boolean 'keepsamples' to omit outliers determined in the sample tree above 
nGenes = ncol(dds.d18_vst) # number of genes == 4936 
nSamples = nrow(dds.d18_vst) # number of samples == 10  - the cut tree removed 1 sample 

# plot the tree with the 'keep samples'  =========================================== #
sampleTree2 = hclust(dist(dds.d18_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
png("Output/WGCNA/day18_spat/day18_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()


# based on outlier removal... call Trait data ===================================================== #
dim(dds.d18_vst) # 10 4936- transformed count data now  has 35 samples - 1 cut in the tree step above 
dim(d18.Treatment.data) # 11 6 - trait data has 11  samples - not yet cut! 

# Form a data frame analogous to expression data that will hold the clinical traits.
d18.Samples    = rownames(dds.d18_vst); # start new variable 'd18.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows      = match(d18.Samples, d18.Treatment.data$Sample.Name); # match the names - calls the list number of 'd18.Samples' matching 'd18.Treatment.data$Sample.Name'
d18.Traits     = d18.Treatment.data[TreatRows, -1]; # insert TreatRows as the row numbers in 'd18.Treatment.data'
rownames(d18.Traits)     = d18.Treatment.data[TreatRows, 1]; # inserts the new TreatRows - matches sample IDs
all(rownames(d18.Traits) == rownames(dds.d18_vst)) # should be TRUE
dim(d18.Traits) # 10 Samples 5 columns; now we have 22 samples! - colnames are all treatment, primary and second treatment


# ===================================================================================
#
# Prepare Trait data (Phys and Treatment groups)
# ===================================================================================

# Salinity groups  ===================================================== #
d18.Traits.Salinity      <-  d18.Traits %>% dplyr::select('Salinity')  %>% # primary treatment as Ambient (A) vs. Moderate (M)
  dplyr::mutate(High = as.factor(as.numeric(Salinity == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
  dplyr::mutate(Low = as.factor(as.numeric(Salinity == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
  dplyr::select(-Salinity)
d18.Traits.Salinity  # final dataset of 0,1 for treatment groups - Primary only!


# pCO2 groups  ===================================================== #
d18.Traits.pCO2          <-  d18.Traits %>% dplyr::select('pCO2')  %>% # primary treatment as Ambient (A) vs. Moderate (M)
  dplyr::mutate(High = as.factor(as.numeric(pCO2 == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
  dplyr::mutate(Low = as.factor(as.numeric(pCO2 == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
  dplyr::select(-pCO2)
d18.Traits.pCO2  # final dataset of 0,1 for treatment groups - Primary only!

# oCO2_Salinity (as _ _ pCO2 and salinity)  ================================================================ #
d18.Traits.pCO2Salinity   <- d18.Traits %>% 
  dplyr::select('pCO2_Salinity') %>% # primary treatment as Ambient (A) vs. Moderate (M)
  dplyr::mutate(HH = as.factor(as.numeric(pCO2_Salinity == "HH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::mutate(HL = as.factor(as.numeric(pCO2_Salinity == "HL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::mutate(LL = as.factor(as.numeric(pCO2_Salinity == "LL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::mutate(LH = as.factor(as.numeric(pCO2_Salinity == "LH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::select(-pCO2_Salinity)
d18.Traits.pCO2Salinity


# ===================================================================================
#
# cluster samples by Trait
# ===================================================================================

# pCO2 ONLY
png("Output/WGCNA/day18_spat/day18_ClusterTree_pCO2.png", 1000, 1000, pointsize=20)
traitColors_pCO2 = labels2colors(d18.Traits.pCO2); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_pCO2, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d18.Traits.pCO2), 
                    main = "Sample dendrogram and trait heatmap (pCO2)")
dev.off()

# Salinity ONLY
png("Output/WGCNA/day18_spat/day18_ClusterTree_Salinity.png", 1000, 1000, pointsize=20)
traitColors_Salinity = labels2colors(d18.Traits.Salinity); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Salinity, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d18.Traits.Salinity), 
                    main = "Sample dendrogram and trait heatmap (Salinity)")
dev.off()


# pCO2 and Salinity only 
png("Output/WGCNA/day18_spat/day18_ClusterTree_pCO2Salinity.png", 1000, 1000, pointsize=20)
traitColors_pCO2Sal = labels2colors(d18.Traits.pCO2Salinity); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_pCO2Sal, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d18.Traits.pCO2Salinity), 
                    main = "Sample dendrogram and trait heatmap (pCO2_Salinity)")
dev.off()


# write some of te data tus far.... ================================== #
# save(dds.d18_vst, d18.Traits, file = "Presentations/URI/2022_URI_Puritz_Genomics_class/d.18-dataInput.RData")
# write.csv(dds.d18_vst, "Presentations/URI/2022_URI_Puritz_Genomics_class/Day18_vstTransformed_WGCNAdata.csv") # # write the vst transformed data 



# ===================================================================================
#
# soft threshold
# ===================================================================================

dim(dds.d18_vst) #  22 4820 - again double check you have the correct data...
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.d18_vst, powerVector = powers, verbose = 5) #...wait for this to finish
#  pickSoftThreshold 
#  performs the analysis of network topology and aids the
#  user in choosing a proper soft-thresholding power.
#  The user chooses a set of candidate powers (the function provides suitable default values)
#  function returns a set of network indices that should be inspected

# png to output 
sizeGrWindow(9, 5) # set window size 
png("Output/WGCNA/day18_spat/day18_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], # Scale-free topology fit index as a function of the soft-thresholding power
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") # look at at cut off at power of 3 - this line corresponds to using an R^2 cut-off of h

plot(sft$fitIndices[,1], sft$fitIndices[,5], # Mean connectivity as a function of the soft-thresholding power
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off() # output 


# The left panel... shows the scale-free fit index (y-axis) 
# as a function of the soft-thresholding power (x-axis). 
# We choose the power 9, which is the lowest power for which the scale-free 
# topology index curve attens out upon reaching a high value (in this case, roughly 0.90).
# The right panel.... displays the mean connectivity
# (degree, y-axis) as a function of the soft-thresholding power (x-axis).


#=====================================================================================
#
#  Satrt the step-wise module construction:  
# Step 1 = create adjacency matrix 
# https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
# https://www.rdocumentation.org/packages/WGCNA/10cpm/versions/1.69/topics/adjacency
# https://ramellose.github.io/networktutorials/wgcna.html
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/10cpm/TechnicalReports/signedTOM.pdf
#=====================================================================================
softPower = 12 # set your soft threshold based on the plots above 

# signed - must call te type, defaults to unsigned
adjacency_sign = adjacency(dds.d18_vst, power = softPower, type="signed") # this takes a long time.. just wait...

#=====================================================================================
#
#  Step 2: Turn adjacency into topological overlap matrix
# Calculation of the topological overlap matrix, (TOM)
# and the corresponding dissimilarity, from a given adjacency matrix.
#=====================================================================================

# signed matrix
TOM_sign       = TOMsimilarity(adjacency_sign, TOMType="signed")  # this takes a long time.. just wait...
dissTOM_sign   = 1-TOM_sign

#=====================================================================================
#
#  Step 3:Call the hierarchical clustering function - plot the tree
#
#=====================================================================================
# Call the hierarchical clustering function
geneTree_sign   = hclust(as.dist(dissTOM_sign), method = "average");

# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree_sign, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity - SIGNED",
     labels = FALSE, hang = 0.04);

#=====================================================================================
#
#  Step 4: Set module size and 'cutreeDynamic' to create clusters 
#
#=====================================================================================

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100; # set this for the subseqent call - WGCNA authors recomend diligence when calling module size to avoid too many/too few modules...
# Module identification using dynamic tree cut:
dynamicMods_sign = cutreeDynamic(dendro = geneTree_sign, distM = dissTOM_sign,
                                 deepSplit = 1, pamRespectsDendro = FALSE,
                                 minClusterSize = minModuleSize);
table(dynamicMods_sign) # view the number of genes per module 
# dynamicMods_sign
# 1    2    3    4    5    6    7    8    9   10   11   12   13   14 
# 1040  808  568  439  291  277  237  225  198  194  178  178  163  140 



#=====================================================================================
#
#  Step 5: convert numeric network to colors and plot the dendrogram
#
#=====================================================================================

# Convert numeric lables into colors
dynamicColors_sign = labels2colors(dynamicMods_sign) # add colors to module labels (previously numbers)
table(dynamicColors_sign) # lets look at this table...
# dynamicColors_sign
# black        blue       brown        cyan       green greenyellow     magenta        pink 
# 237         808         568         140         291         178         198         225 
# purple         red      salmon         tan   turquoise      yellow 
# 194         277         163         178        1040         439 
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
png("Output/WGCNA/day18_spat/day18_Dendrogram.png", 1000, 1000, pointsize=20)
plotDendroAndColors(geneTree_sign, dynamicColors_sign, "Dynamic Tree Cut - SIGNED",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors 'SIGNED'")
dev.off()

#=====================================================================================
#
#  Step 6: Calculate Eigengenes - view thier connectivity based on 'MEDiss = 1-cor(MEs)'
#
#=====================================================================================

# Calculate eigengenes ========================================================== # 

# MEList = moduleEigengenes(dds.d18_vst, colors = dynamicColors)
MEList = moduleEigengenes(dds.d18_vst, colors = dynamicColors_sign)
MEs    = MEList$eigengenes # you can view MEs, condenses gene counts down to a single number for each sample representive of that expressoin pattern 

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average") # Cluster module eigengenes

# Plot the result ================================================================ #
sizeGrWindow(18, 6)
png("Output/WGCNA/day18_spat/day18_ClusterEigengenes.png", 1000, 1000, pointsize=20)
plot(METree, main = "Clustering of module eigengenes - SIGNED (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
dev.off() # we will keep all of these , no apparent need to condense related nodes (merge modules)


#=====================================================================================
#
#  Step 7: Specify the cut line for the dendrogram (module) - Calc MODULE EIGENGENES (mergeMEs)
#
#=====================================================================================


# change the cut height for the eignengene tree! 
MEDissThres = 0.4 # call the cut height

# save with the abline added
plot(METree, main = "Clustering of module eigengenes (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "") +
  abline(h=MEDissThres, col = "red") # add the line to the tree

# Call an automatic merging function
merge = mergeCloseModules(dds.d18_vst, dynamicColors_sign, cutHeight = MEDissThres, verbose = 3) # signed TOM 

# The merged module colors
mergedColors = merge$colors;

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

# Calculate dissimilarity of module eigengenes
MEDiss_merged = 1-cor(mergedMEs);
METree_merged = hclust(as.dist(MEDiss_merged), method = "average") # Cluster module eigengenes

# count of merged eigengene modules
table(merge$colors)
# black        blue       brown        cyan       green greenyellow      purple         red         tan 
# 435         808         568         140         516         178        1836         277         178 

nrow(geneInfo_GROUPS %>% dplyr::filter(moduleColor %in% 'blue' & p.MM.blue < 0.05))
View(geneInfo_GROUPS)
sizeGrWindow(18, 6)
png("Output/WGCNA/day18_spat/day18_ClusterEigengenes_merged.png", 1000, 1000, pointsize=20)
plot(METree_merged, main = "Clustering of module eigengenes - SIGNED (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
dev.off() # we will keep all of these , no apparent need to condense related nodes (merge modules)




#=====================================================================================
#
#  Step 8: Plot dendrogram with the cut line 'MEDissThres' 
#
#=====================================================================================

# use 'mergedColors' we needed to merge related eigengene modules together 

sizeGrWindow(12, 9)
png("Output/WGCNA/day18_spat/day18_ClusterDendrogram_merged.png", 1000, 1000, pointsize=20)
plotDendroAndColors(geneTree_sign, cbind(dynamicColors_sign, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


#=====================================================================================
#
#  Step 9: Commit to mergedcolors as 'MEs' and 'moduleColors'
#
#=====================================================================================
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;

# Save module colors and labels for use in subsequent parts
save(mergedMEs, dds.d18_vst, moduleLabels, moduleColors, file = "Output/WGCNA/day18_spat/day18-networkConstruction-stepByStep.RData")
# write csv - save the module eigengenes
write.csv(mergedMEs, file = "Output/WGCNA/day18_spat/d18.WGCNA_ModulEigengenes.csv")


#=====================================================================================
#
#  Prepare for  module trait associations - Eigengene calc - trait data as numeric
#
#=====================================================================================
# identify modules that are signiFcantly associated with the measured  traits (here as treatment history)

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the  significant associations:

load("Output/WGCNA/day18_spat/day18-networkConstruction-stepByStep.RData")

# Define numbers of genes and samples
nGenes   = ncol(dds.d18_vst); # 4936
nSamples = nrow(dds.d18_vst); # 10
# Recalculate MEs with color labels
MEs0           <-  read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulEigengenes.csv") # read merged eigengene dataset
rownames(MEs0) <- MEs0[,1] # make first column into row names
MEs0           <- MEs0[,-1] # omit the first column (now inserted as rownames...)
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)


#=====================================================================================
#
# Module trait correlation
#
#=====================================================================================
# ALL TRAIT DATA

dim(d18.Traits)  # 10 5
dim(MEs)  # 10 9
# moduleTraitCor = cor(MEs, d18.Traits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# pCO2
d18.Traits.pCO2.asnum  <- data.frame(lapply(d18.Traits.pCO2, function(x) as.numeric(as.character(x))),
                                    check.names=F, row.names = row.names(d18.Traits.pCO2))
moduleTraitCor_pCO2    = cor(MEs, d18.Traits.pCO2.asnum, use = "p");
moduleTraitPvalue_pCO2 = corPvalueStudent(moduleTraitCor_pCO2, nSamples);
moduleTraitPvalue_pCO2


# Salinity
d18.Traits.Salinity.asnum  <- data.frame(lapply(d18.Traits.Salinity, function(x) as.numeric(as.character(x))),
                                        check.names=F, row.names = row.names(d18.Traits.Salinity))
moduleTraitCor_Salinity    = cor(MEs, d18.Traits.Salinity.asnum, use = "p");
moduleTraitPvalue_Salinity = corPvalueStudent(moduleTraitCor_Salinity, nSamples);
moduleTraitPvalue_Salinity

# pCO2 and Salinity 
d18.Traits.pCO2Salinity.asnum  <- data.frame(lapply(d18.Traits.pCO2Salinity, function(x) as.numeric(as.character(x))),
                                            check.names=F, row.names = row.names(d18.Traits.pCO2Salinity))
moduleTraitCor_pCO2Salinity    = cor(MEs, d18.Traits.pCO2Salinity.asnum, use = "p");
moduleTraitPvalue_pCO2Salinity = corPvalueStudent(moduleTraitCor_pCO2Salinity, nSamples);
moduleTraitPvalue_pCO2Salinity


#=====================================================================================
#
# Heatmaps
#
#=====================================================================================



# pCO2 ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d18.pCO2.matrix <-  paste(signif(moduleTraitCor_pCO2, 2), "\n(",
                         signif(moduleTraitPvalue_pCO2, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day18_spat/heatmaps/day18_pCO2_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_pCO2,
               xLabels = names(d18.Traits.pCO2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d18.pCO2.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - pCO2"))
dev.off()

# this heatmap looks better
d18.pCO2.text <-  as.matrix(signif(moduleTraitPvalue_pCO2, 3))
pa                  = cluster::pam(d18.pCO2.text, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day18_spat/heatmaps/day18_pCO2_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_pCO2, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 2 WGCNA - pCO2", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 1,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d18.pCO2.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()











# Salinity ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d18.Salinity.matrix <-  paste(signif(moduleTraitCor_Salinity, 2), "\n(",
                             signif(moduleTraitPvalue_Salinity, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day18_spat/heatmaps/day18_Salinity_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Salinity,
               xLabels = names(d18.Traits.Salinity),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d18.Salinity.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Salinity"))
dev.off()

# this heatmap looks better
d18.Salinity.text <-  as.matrix(signif(moduleTraitPvalue_Salinity, 3))
pa                  = cluster::pam(d18.Salinity.text, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day18_spat/heatmaps/day18_Salinity_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Salinity, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 2 WGCNA - Salinity", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 1,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d18.Salinity.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()












# pCO2 and Salinity  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d18.pCO2Salinity.matrix <-  paste(signif(moduleTraitCor_pCO2Salinity, 2), "\n(",
                                 signif(moduleTraitPvalue_pCO2Salinity, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day18_spat/heatmaps/day18_pCO2Salinity_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_pCO2Salinity,
               xLabels = names(d18.Traits.pCO2Salinity),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d18.pCO2Salinity.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - pCO2Salinity"))
dev.off()

# this heatmap looks better
d18.pCO2Salinity.text <-  as.matrix(signif(moduleTraitPvalue_pCO2Salinity, 3))
pa                  = cluster::pam(d18.pCO2Salinity.text, k = 4)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day18_spat/heatmaps/day18_pCO2Salinity_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_pCO2Salinity, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 2 WGCNA - pCO2Salinity", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 2,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d18.pCO2Salinity.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()


#=====================================================================================
#
# Module eigengene -  MEs boxplots by treatment group
#
#=====================================================================================
MEs_table             <- mergedMEs # new table for plotting 
MEs_table$Sample.Name <- rownames(mergedMEs) # call rows as coolumn to merge with treatment data
MEsPlotting           <- merge(d18.Treatment.data, MEs_table, by = 'Sample.Name') # merge
MEsPlotting           <- MEsPlotting[,-c(2,5)] # ommit the all treatments column
MEsPlotting_melt      <- melt(MEsPlotting, id=c('Sample.Name', 'pCO2', 'Salinity', 'pCO2_Salinity'))

#plot it
png("Output/WGCNA/day18_spat/Day18_ME_Boxplot.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Salinity, y=value, fill = factor(pCO2), shape=pCO2)) +
  geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.185),
               geom = "point", shape = 19, size = 3,
               show.legend = FALSE) +
  ylab("ModuleEigengene") +
  ylim(-0.5,0.5) +
  scale_fill_manual(values=c("#D55E00","#56B4E9")) +
  geom_hline(yintercept=0, linetype='dotted', col = 'black', size = 1)+
  theme_bw() +
  #theme(legend.position = "none") +
  facet_wrap(~variable)
dev.off()





#=====================================================================================
#
# geneModuleMembership - geneTraitSignificance - GSPvalue
#
#=====================================================================================
# Gene relationship to trait and important modules: 
# Gene Significance and Module Membership

# We quantify associations of individual genes with our trait of interest (TAOC)
load("Output/WGCNA/day18_spat/day18-networkConstruction-stepByStep.RData")
# names (colors) of the modules
modNames = substring(names(mergedMEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(dds.d18_vst, mergedMEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


# HH treatment group
HH = as.data.frame(as.numeric(d18.Traits.pCO2Salinity$HH)); # Define variable containing the desired column 
names(HH) = "HH"
HH_geneTraitSignificance = as.data.frame(cor(dds.d18_vst, HH, use = "p"));
HH_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(HH_geneTraitSignificance), nSamples));
names(HH_geneTraitSignificance) = paste("GS.", names(HH), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(HH_GSPvalue) = paste("p.GS.", names(HH), sep=""); # corPvalueStudent

#  PLOT mean.µmol.CRE.g.protein in the MAGENTA module
# unique(moduleColors)
# module = "brown"
# column = match(module, modNames);
# moduleGenes = moduleColors==module;
# 
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(MA_geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for 'MA' treatment group",
#                    main = paste("day14 'MA' treatment group  'BROWN': Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#=====================================================================================
#
#   COUNT GENES OF INTEREST IN  MODULES (i.e. violet- refer to heatmap)
#
#=====================================================================================

length(colnames(dds.d18_vst)[moduleColors=="blue"]) # 808 total genes in the blue module

#=====================================================================================
#
#  Call annotation data to get module gene data (prep for downstream GO)
#
#=====================================================================================
annot = read.csv(file = "Data/TagSeq/Seq_details/Seq_Reference_Master.csv",header = T) 
dim(annot) # 66437     8
names(annot) # view the column names to call 
probes = names(as.data.frame(t(d18.data_matrix[, -(1)])))
probes2annot = match(probes, annot$Cvirginica_TranscriptID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot)) # 380
# Should return 0.
#=====================================================================================
#
#  BUILD GENE INFO DATAFRAMES
#
#=====================================================================================
# Create the starting data frame
names(annot)

#   dataframe  --------------------------------------------------------------------------- # 

geneInfo_GROUPS = data.frame(geneSymbol       = annot$Cvirginica_GeneID[probes2annot],
                             TranscriptID     = annot$Cvirginica_TranscriptID[probes2annot],
                             moduleColor      = moduleColors,
                             KEGG_ID          = annot$Cgigas_KEGGID[probes2annot],
                             Protein_name     = annot$Cvirginica_Protein_name[probes2annot],
                             gene_length      = annot$Cvirginica_length[probes2annot],
                             GO.terms         = annot$Annotation_GO_ID[probes2annot],
                             HH_geneTraitSignificance, HH_GSPvalue)
# call this specific to the module and trait of interest
View(geneInfo_GROUPS)
modOrder = order(-abs(cor(mergedMEs, HH, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) { # Add module membership information in the chosen order
  
  oldNames = names(geneInfo_GROUPS)
  geneInfo_GROUPS = data.frame(geneInfo_GROUPS, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]]);
  names(geneInfo_GROUPS) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_GROUPS$moduleColor, -abs(geneInfo_GROUPS$GS.HH));
geneInfo_GROUPS = geneInfo_GROUPS[geneOrder, ]
View(geneInfo_GROUPS)


#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 

write.csv(geneInfo_GROUPS, file = "Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv")

