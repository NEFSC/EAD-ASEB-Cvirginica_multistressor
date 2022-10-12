---
# title: "Day7_WGCNA_all"
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

d2.filtered_count_tble  <- read.csv(file="Data/TagSeq/Filtered_counts/filtered_counts_5cpm_50perc/day2.filtered_5cpm50perc.csv", sep=',', header=TRUE) 

# d7.filtered_count_data          <- read.csv(file ="Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_countmatrix_WGCNA.csv", sep =',', header=TRUE)
# Master.Treatment_Phenotype.data <- data.frame(read.csv(file="Presentations/URI/2022_URI_Puritz_Genomics_class/Master_Phyenotype.and.Exp.Treatment_Metadata.csv", sep=',', header=TRUE))

# The following setting is important, do not omit. (as recommended by WGCNA authors - view tutorial)
options(stringsAsFactors = FALSE)

# ===================================================================================
#
#
# Day 2 WGCNA - PREPROCESSING THE DATA INPUT 
#
#  Read here for the pre processing steps using WGCNA!
#  https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-9-559.pdf
# ===================================================================================

# count data  ========================================================== #

d2.data_matrix <- data.frame(d2.filtered_count_tble[,-1], row.names=d2.filtered_count_tble[,1]) 
dim(d2.data_matrix) # 4820   24 - 24 samples and 4820 total genes 


# trait data ========================================================== #

# Treatment data

d2.Treatment.data <- read.csv(file="Data/TagSeq/day2.exp.data.csv", sep=',', header=TRUE) %>%  
                        dplyr::mutate_if(is.character, as.factor) %>% 
                        dplyr::rename('Sample.Name' = 'SapleName_readmatrix') %>% 
                        dplyr::rename('pCO2' = 'OA') %>% 
                        dplyr::select(c('Sample.Name','Temperature','pCO2','Salinity', 'Aragonite_saturation')) %>% 
                        dplyr::mutate(All_treatment = paste( (substr(Temperature,1,1)), 
                                                             (substr(pCO2,1,1)), 
                                                             (substr(Salinity,1,1)), sep = '')) %>% 
                        dplyr::mutate(pCO2_Salinity = substr(All_treatment, 2,3)) %>% 
                        dplyr::mutate(Aragonite_saturation = case_when(Aragonite_saturation < 0.5 ~ 'Low', 
                                                                       (Aragonite_saturation > 0.5 & Aragonite_saturation < 1.0 & Salinity == 'High') ~ 'Mid_CO2', # assign 'mid' category as  pCO2 driven (Salinity is HIGH, pCO2 HIGH is driving the reduced Arag)
                                                                       (Aragonite_saturation > 0.5 & Aragonite_saturation < 1.0 & Salinity == 'Low') ~ 'Mid_Sal', # assign 'mid' category as  Salinity driven (Salinity is LOW and pCO2 LOW, thus Salinty is driving low  Arag)
                                                                       Aragonite_saturation > 1.0 ~ 'High'))

View(d2.Treatment.data)


# sanity check =========== # 

dim(d2.Treatment.data)[1] ==  dim(d2.data_matrix)[2]# TRUE - each contains all 36 samples sequenced for Day 7 of the experiment 


# ===================================================================================
#
# Day 2 DESeqDataSet or 'dds' object (using DESeq2) 
#       vst transform the 'dds' onject for WGCNA 
# 
# ===================================================================================

# create dds  ========================================================== #

# NOTE: ~1 stands for no design; user will need to add a design for differential testing
# however for our purpose of just creating an object to transform, we do not need a design here...
dds.d2 <- DESeqDataSetFromMatrix(countData = d2.data_matrix,
                                 colData = d2.Treatment.data, design = ~ 1) # DESeq Data Set (dds)
dds.d2 # view the DESeqDataSet - notice the colData containg our critical treatment and sample ID data, rownames, etc. 

# transform the data  ========================================================== #
# run in order (kept name as dds.d2_vst)
dds.d2_vst <- vst(dds.d2) # transform it vst (variance stabilized transformation)
dds.d2_vst <- assay(dds.d2_vst) # call only the transformed coutns in the dds object
#fix(dds.d2_vst)
dds.d2_vst <- t(dds.d2_vst) # transpose columns to rows and vice versa

# ===================================================================================
#
# Day 2 WGCNA Sample tree - omit outlier sample(s)
#
# ===================================================================================

# checks before we start....
dim(dds.d2_vst) #  8548 genes; 36  samples
gsg = goodSamplesGenes(dds.d2_vst, verbose = 3);  # We first check for genes and samples with too many missing values:
gsg$allOK # If the statement returns TRUE, all genes have passed the cuts. 

# determine outlier(s) with a sample tree    ========================================================== #

# call the cluster and set window dimenstions to view..
sampleTree = hclust(dist(dds.d2_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
sizeGrWindow(12,9) # The user should change the dimensions if the window is too large or too small.
par(cex = 0.6);
par(mar = c(0,4,2,0))
# output the tree as .png file..
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2) # appears there are two outliers SG59; can remove by hand or an automatic appraoch 
abline(h = 11, col = "red") # add line to plot to show the cut-off od outlier samples (40000) SG105 and SG55
dev.off()

# view your tree   ========================================================== #
# NOTE: I inserted a horizontal line at 110 demonstrating where we likely want to cut our tree - omitting a single outlier 

# cut the tree and omit  ========================================================== #
clust = cutreeStatic(sampleTree, cutHeight = 11, minSize = 10) # Determine cluster under the line
table(clust) # 0 = cut; 1 = kept; says it will cut 1 and save 35; exactly what we want!  
keepSamples = (clust==1) # 'keepsamples' boolean to call the main dataset - notice there are TWO occurrences of FALSE - these are C6.larva and B12.larva

# integrate keepsamples  ========================================================== #
dds.d2_vst = dds.d2_vst[keepSamples, ] # integrate the boolean 'keepsamples' to omit outliers determined in the sample tree above 
nGenes = ncol(dds.d2_vst) # number of genes == 4820 
nSamples = nrow(dds.d2_vst) # number of samples == 22  - the cut tree removed 2 samples 

# plot the tree with the 'keep samples'  =========================================== #
sampleTree2 = hclust(dist(dds.d2_vst), method = "average") # Next we cluster the samples (in contrast to clustering genes that will come later)  to see if there are any obvious outliers.
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
dev.off()


# based on outlier removal... call Trait data ===================================================== #
dim(dds.d2_vst) # 22 4820- transformed count data now  has 35 samples - 1 cut in the tree step above 
dim(d2.Treatment.data) # 24 4 - trait data has  36 samples - not yet cut! 

# Form a data frame analogous to expression data that will hold the clinical traits.
d2.Samples    = rownames(dds.d2_vst); # start new variable 'd7.Samples' calling the row names of the gene data (sample as 'Mouse' in trait data)
TreatRows     = match(d2.Samples, d2.Treatment.data$Sample.Name); # match the names - calls the list number of 'd7.Samples' matching 'd7.Treatment.data$Sample.Name'
d2.Traits     = d2.Treatment.data[TreatRows, -1]; # insert TreatRows as the row numbers in 'd7.Treatment.data'
rownames(d2.Traits)     = d2.Treatment.data[TreatRows, 1]; # inserts the new TreatRows - matches sample IDs
all(rownames(d2.Traits) == rownames(dds.d2_vst)) # should be TRUE
dim(d2.Traits) # 22 Samples 6 columns; now we have 22 samples! - colnames are all treatment, primary and second treatment


# ===================================================================================
#
# Prepare Trait data (Phys and Treatment groups)
# ===================================================================================


# temperature groups ======================================================== #
d2.Traits.Temperature  <-  d2.Traits %>% dplyr::select('Temperature') %>% # primary treatment as Ambient (A) vs. Moderate (M)
                               dplyr::mutate(High = as.factor(as.numeric(Temperature == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
                               dplyr::mutate(Low = as.factor(as.numeric(Temperature == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
                               dplyr::select(-Temperature)
d2.Traits.Temperature  # final dataset of 0,1 for treatment groups - temperature only!


# Salinity groups  ===================================================== #
d2.Traits.Salinity      <-  d2.Traits %>% dplyr::select('Salinity')  %>% # primary treatment as Ambient (A) vs. Moderate (M)
  dplyr::mutate(High = as.factor(as.numeric(Salinity == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
  dplyr::mutate(Low = as.factor(as.numeric(Salinity == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
  dplyr::select(-Salinity)
d2.Traits.Salinity  # final dataset of 0,1 for treatment groups - Primary only!


# pCO2 groups  ===================================================== #
d2.Traits.pCO2          <-  d2.Traits %>% dplyr::select('pCO2')  %>% # primary treatment as Ambient (A) vs. Moderate (M)
                                dplyr::mutate(High = as.factor(as.numeric(pCO2 == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
                                dplyr::mutate(Low = as.factor(as.numeric(pCO2 == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
                                dplyr::select(-pCO2)
d2.Traits.pCO2  # final dataset of 0,1 for treatment groups - Primary only!


# pCO2 groups  ===================================================== #
d2.Traits.AragoniteSat <-  d2.Traits %>% dplyr::select('Aragonite_saturation')  %>% # primary treatment as Ambient (A) vs. Moderate (M)
  dplyr::mutate(High    = as.factor(as.numeric(Aragonite_saturation == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
  dplyr::mutate(Mid_CO2 = as.factor(as.numeric(Aragonite_saturation == "Mid_CO2")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
  dplyr::mutate(Mid_Sal = as.factor(as.numeric(Aragonite_saturation == "Mid_Sal")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
  dplyr::mutate(Low     = as.factor(as.numeric(Aragonite_saturation == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
  dplyr::select(-Aragonite_saturation)
d2.Traits.AragoniteSat  # final dataset of 0,1 for treatment groups - Primary only!



# oCO2_Salinity (as _ _ pCO2 and salinity)  ================================================================ #
d2.Traits.pCO2Salinity   <- d2.Traits %>% 
  dplyr::select('pCO2_Salinity') %>% # primary treatment as Ambient (A) vs. Moderate (M)
  dplyr::mutate(HH = as.factor(as.numeric(pCO2_Salinity == "HH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::mutate(HL = as.factor(as.numeric(pCO2_Salinity == "HL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::mutate(LL = as.factor(as.numeric(pCO2_Salinity == "LL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
  dplyr::mutate(LH = as.factor(as.numeric(pCO2_Salinity == "LH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
dplyr::select(-pCO2_Salinity)
d2.Traits.pCO2Salinity


# all treatment grousp (as _ _ _ for temp, pCO2 and salinity)  ================================================================ #
d2.Traits.Group         <- d2.Traits %>% 
                              dplyr::select('All_treatment') %>% # primary treatment as Ambient (A) vs. Moderate (M)
                              dplyr::mutate(HHH = as.factor(as.numeric(All_treatment == "HHH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(HHL = as.factor(as.numeric(All_treatment == "HHL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(HLL = as.factor(as.numeric(All_treatment == "HLL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(HLH = as.factor(as.numeric(All_treatment == "HLH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(LLL = as.factor(as.numeric(All_treatment == "LLL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(LHH = as.factor(as.numeric(All_treatment == "LHH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(LHL = as.factor(as.numeric(All_treatment == "LHL")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::mutate(LLH = as.factor(as.numeric(All_treatment == "LLH")))  %>%  # call occurrence of 'AA' as 0s and 1s (factor)
                              dplyr::select(-All_treatment)
d2.Traits.Group


# ===================================================================================
#
# cluster samples by Trait
# ===================================================================================

# Temperature ONLY
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_Temperature.png", 1000, 1000, pointsize=20)
traitColors_Temperature = labels2colors(d2.Traits.Temperature); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Temperature, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d2.Traits.Temperature), 
                    main = "Sample dendrogram and trait heatmap (temperature)")
dev.off()

# pCO2 ONLY
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_pCO2.png", 1000, 1000, pointsize=20)
traitColors_pCO2 = labels2colors(d2.Traits.pCO2); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_pCO2, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d2.Traits.pCO2), 
                    main = "Sample dendrogram and trait heatmap (pCO2)")
dev.off()

# Salinity ONLY
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_Salinity.png", 1000, 1000, pointsize=20)
traitColors_Salinity = labels2colors(d2.Traits.Salinity); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Salinity, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d2.Traits.Salinity), 
                    main = "Sample dendrogram and trait heatmap (Salinity)")
dev.off()

# Aragonite saturation ONLY
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_AragoniteSat.png", 1000, 1000, pointsize=20)
traitColors_AragoniteSat = labels2colors(d2.Traits.AragoniteSat); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_AragoniteSat, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d2.Traits.AragoniteSat), 
                    main = "Sample dendrogram and trait heatmap (AragoniteSat)")
dev.off()

# pCO2 and Salinity only 
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_pCO2Salinity.png", 1000, 1000, pointsize=20)
traitColors_pCO2Sal = labels2colors(d2.Traits.pCO2Salinity); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_pCO2Sal, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d2.Traits.pCO2Salinity), 
                    main = "Sample dendrogram and trait heatmap (pCO2_Salinity)")
dev.off()


# Group ONLY
png("Output/WGCNA/day2_larvae/Day2_ClusterTree_Group.png", 1000, 1000, pointsize=20)
traitColors_Group = labels2colors(d2.Traits.Group); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Group, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d2.Traits.Group), 
                    main = "Sample dendrogram and trait heatmap (Group)")
dev.off()


# write some of te data tus far.... ================================== #
# save(dds.d7_vst, d7.Traits, file = "Presentations/URI/2022_URI_Puritz_Genomics_class/d.7-dataInput.RData")
# write.csv(dds.d7_vst, "Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_vstTransformed_WGCNAdata.csv") # # write the vst transformed data 



# ===================================================================================
#
# soft threshold
# ===================================================================================

dim(dds.d2_vst) #  22 4820 - again double check you have the correct data...
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.d2_vst, powerVector = powers, verbose = 5) #...wait for this to finish
#  pickSoftThreshold 
#  performs the analysis of network topology and aids the
#  user in choosing a proper soft-thresholding power.
#  The user chooses a set of candidate powers (the function provides suitable default values)
#  function returns a set of network indices that should be inspected

# png to output 
sizeGrWindow(9, 5) # set window size 
png("Output/WGCNA/day2_larvae/Day2_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)
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
softPower = 10 # set your soft threshold based on the plots above 

# signed - must call te type, defaults to unsigned
adjacency_sign = adjacency(dds.d2_vst, power = softPower, type="signed") # this takes a long time.. just wait...


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
# 1   2   3   4   5   6   7   8 
# 881 770 690 645 518 506 431 379

#=====================================================================================
#
#  Step 5: convert numeric network to colors and plot the dendrogram
#
#=====================================================================================

# Convert numeric lables into colors
dynamicColors_sign = labels2colors(dynamicMods_sign) # add colors to module labels (previously numbers)
table(dynamicColors_sign) # lets look at this table...
# dynamicColors_sign
# black      blue     brown     green      pink       red turquoise    yellow 
# 431       770       690       518       379       506       881       645 
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
png("Output/WGCNA/day2_larvae/day2_Dendrogram.png", 1000, 1000, pointsize=20)
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

# MEList = moduleEigengenes(dds.d7_vst, colors = dynamicColors)
MEList = moduleEigengenes(dds.d2_vst, colors = dynamicColors_sign)
MEs    = MEList$eigengenes # you can view MEs, condenses gene counts down to a single number for each sample representive of that expressoin pattern 

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average") # Cluster module eigengenes

# Plot the result ================================================================ #
sizeGrWindow(7, 6)
png("Output/WGCNA/day2_larvae/Day2_ClusterEigengenes.png", 1000, 1000, pointsize=20)
plot(METree, main = "Clustering of module eigengenes - SIGNED (dissimilarity calc = MEDiss = 1-cor(MEs))",
     xlab = "", sub = "")
dev.off()


#=====================================================================================
#
#  Step 7: Specify the cut line for the dendrogram (module) - Calc MODULE EIGENGENES (mergeMEs)
#
#=====================================================================================


# not necessary for this dataset - METree above looks to be parsed well into the 5 modules 

#=====================================================================================
#
#  Step 8: Plot dendrogram with the cut line 'MEDissThres' 
#
#=====================================================================================

# use 'mergedColors' we needed to merge related eigengene modules together 


# sizeGrWindow(12, 9)
# png("Analysis/Output/WGCNA/10cpm/Day7/Day7_ClusterDendrogram_signed.png", 1000, 1000, pointsize=20)
# plotDendroAndColors(geneTree_sign, cbind(dynamicColors_sign, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()


#=====================================================================================
#
#  Step 9: Commit to mergedcolors as 'MEs' and 'moduleColors'
#
#=====================================================================================
# Rename to moduleColors
moduleColors = dynamicColors_sign
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;

# Save module colors and labels for use in subsequent parts
save(MEs, dds.d2_vst, moduleLabels, moduleColors, file = "Output/WGCNA/day2_larvae/Day2-networkConstruction-stepByStep.RData")
# write csv - save the module eigengenes
write.csv(MEs, file = "Output/WGCNA/day2_larvae/d2.WGCNA_ModulEigengenes.csv")


#=====================================================================================
#
#  Prepare for  module trait associations - Eigengene calc - trait data as numeric
#
#=====================================================================================
# identify modules that are signiFcantly associated with the measured  traits (here as treatment history)

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the  significant associations:

# Define numbers of genes and samples
nGenes = ncol(dds.d2_vst); # 8548
nSamples = nrow(dds.d2_vst); # 35
# Recalculate MEs with color labels
# Recalculate MEs with color labels
MEs0           <-  read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulEigengenes.csv") # read merged eigengene dataset
rownames(MEs0) <- MEs0[,1] # make first column into row names
MEs0           <- MEs0[,-1] # omit the first column (now inserted as rownames...)
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)

#=====================================================================================
#
# Module trait correlation
#
#=====================================================================================
# ALL TRAIT DATA

dim(d2.Traits)  # 22  6
dim(MEs)  # 22  8
# moduleTraitCor = cor(MEs, d7.Traits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# temperature 
d2.Traits.Temperature.asnum    <- data.frame(lapply(d2.Traits.Temperature, function(x) as.numeric(as.character(x))),
                                check.names=F, row.names = row.names(d2.Traits.Temperature))
moduleTraitCor_Temperature      = cor(MEs, d2.Traits.Temperature.asnum, use = "p");
moduleTraitPvalue_Temperature   = corPvalueStudent(moduleTraitCor_Temperature, nSamples);
moduleTraitPvalue_Temperature


# pCO2
d2.Traits.pCO2.asnum  <- data.frame(lapply(d2.Traits.pCO2, function(x) as.numeric(as.character(x))),
                                       check.names=F, row.names = row.names(d2.Traits.pCO2))
moduleTraitCor_pCO2    = cor(MEs, d2.Traits.pCO2.asnum, use = "p");
moduleTraitPvalue_pCO2 = corPvalueStudent(moduleTraitCor_pCO2, nSamples);
moduleTraitPvalue_pCO2


# Salinity
d2.Traits.Salinity.asnum  <- data.frame(lapply(d2.Traits.Salinity, function(x) as.numeric(as.character(x))),
                                    check.names=F, row.names = row.names(d2.Traits.Salinity))
moduleTraitCor_Salinity    = cor(MEs, d2.Traits.Salinity.asnum, use = "p");
moduleTraitPvalue_Salinity = corPvalueStudent(moduleTraitCor_Salinity, nSamples);
moduleTraitPvalue_Salinity


# Aragonite.Sat
d2.Traits.AragoniteSat.asnum  <- data.frame(lapply(d2.Traits.AragoniteSat, function(x) as.numeric(as.character(x))),
                                        check.names=F, row.names = row.names(d2.Traits.AragoniteSat))
moduleTraitCor_AragoniteSat    = cor(MEs, d2.Traits.AragoniteSat.asnum, use = "p");
moduleTraitPvalue_AragoniteSat = corPvalueStudent(moduleTraitCor_AragoniteSat, nSamples);
moduleTraitPvalue_AragoniteSat


# pCO2 and Salinity 
d2.Traits.pCO2Salinity.asnum  <- data.frame(lapply(d2.Traits.pCO2Salinity, function(x) as.numeric(as.character(x))),
                                        check.names=F, row.names = row.names(d2.Traits.pCO2Salinity))
moduleTraitCor_pCO2Salinity    = cor(MEs, d2.Traits.pCO2Salinity.asnum, use = "p");
moduleTraitPvalue_pCO2Salinity = corPvalueStudent(moduleTraitCor_pCO2Salinity, nSamples);
moduleTraitPvalue_pCO2Salinity


# Group
d2.Traits.Group.asnum  <- data.frame(lapply(d2.Traits.Group, function(x) as.numeric(as.character(x))),
                                        check.names=F, row.names = row.names(d2.Traits.Group))
moduleTraitCor_Group    = cor(MEs, d2.Traits.Group.asnum, use = "p");
moduleTraitPvalue_Group = corPvalueStudent(moduleTraitCor_Group, nSamples);
moduleTraitPvalue_Group

#=====================================================================================
#
# Heatmaps
#
#=====================================================================================



# TEMPERATURE ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d2.TEMPERATURE.matrix <-  paste(signif(moduleTraitCor_Temperature, 2), "\n(",
                                      signif(moduleTraitPvalue_Temperature, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day2_larvae/heatmaps/Day2_Temperature_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Temperature,
               xLabels = names(d2.Traits.Temperature),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d2.TEMPERATURE.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - temperature"))
dev.off()

# this heatmap looks better
d2.TEMPERATURE.text <-  as.matrix(signif(moduleTraitPvalue_Temperature, 3))
d2.TEMPERATURE.COR  <-  as.matrix(signif(moduleTraitCor_Temperature, 3))
pa                  = cluster::pam(d2.TEMPERATURE.COR, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day2_larvae/heatmaps/Day2_Temperature_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Temperature, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 2 WGCNA - Temperature", 
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
          grid.text(sprintf("%.1f", d2.TEMPERATURE.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()








# pCO2 ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d2.pCO2.matrix <-  paste(signif(moduleTraitCor_pCO2, 2), "\n(",
                                signif(moduleTraitPvalue_pCO2, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day2_larvae/heatmaps/Day2_pCO2_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_pCO2,
               xLabels = names(d2.Traits.pCO2),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d2.pCO2.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - pCO2"))
dev.off()

# this heatmap looks better
d2.pCO2.text  <-  as.matrix(signif(moduleTraitPvalue_pCO2, 3))
d2.pCO2.COR   <-  as.matrix(signif(moduleTraitCor_pCO2, 3))
pa                  = cluster::pam(d2.pCO2.COR, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day2_larvae/heatmaps/Day2_pCO2_heatmap.pdf", width=5, height=6)
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
          grid.text(sprintf("%.1f", d2.pCO2.COR[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()











# Salinity ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d2.Salinity.matrix <-  paste(signif(moduleTraitCor_Salinity, 2), "\n(",
                         signif(moduleTraitPvalue_Salinity, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day2_larvae/heatmaps/Day2_Salinity_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Salinity,
               xLabels = names(d2.Traits.Salinity),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d2.Salinity.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Salinity"))
dev.off()

# this heatmap looks better
d2.Salinity.text     <-  as.matrix(signif(moduleTraitPvalue_Salinity, 3))
d2.Salinity.COR      <-  as.matrix(signif(moduleTraitCor_Salinity, 3))
pa                  = cluster::pam(d2.Salinity.COR, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day2_larvae/heatmaps/Day2_Salinity_heatmap.pdf", width=5, height=6)
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
          grid.text(sprintf("%.1f", d2.Salinity.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()











# AragoniteSat ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d2.AragoniteSat.matrix <-  paste(signif(moduleTraitCor_AragoniteSat, 2), "\n(",
                         signif(moduleTraitPvalue_AragoniteSat, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day2_larvae/heatmaps/Day2_AragoniteSat_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_AragoniteSat,
               xLabels = names(d2.Traits.AragoniteSat),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d2.AragoniteSat.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - AragoniteSat"))
dev.off()

# this heatmap looks better
d2.AragoniteSat.text <-  as.matrix(signif(moduleTraitPvalue_AragoniteSat, 3))
d2.AragoniteSat.COR  <-  as.matrix(signif(moduleTraitCor_AragoniteSat, 3))
pa                  = cluster::pam(d2.AragoniteSat.COR, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day2_larvae/heatmaps/Day2_AragoniteSat_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_AragoniteSat, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 2 WGCNA - AragoniteSat", 
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
          grid.text(sprintf("%.1f", d2.AragoniteSat.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()












# pCO2 and Salinity  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d2.pCO2Salinity.matrix <-  paste(signif(moduleTraitCor_pCO2Salinity, 2), "\n(",
                             signif(moduleTraitPvalue_pCO2Salinity, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day2_larvae/heatmaps/Day2_pCO2Salinity_heatmap.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_pCO2Salinity,
               xLabels = names(d2.Traits.pCO2Salinity),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d2.pCO2Salinity.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - pCO2Salinity"))
dev.off()

# this heatmap looks better
d2.pCO2Salinity.text <-  as.matrix(signif(moduleTraitPvalue_pCO2Salinity, 3))
d2.pCO2Salinity.COR  <-  as.matrix(signif(moduleTraitCor_pCO2Salinity, 3))
pa                  = cluster::pam(d2.pCO2Salinity.COR, k = 4)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day2_larvae/heatmaps/Day2_pCO2Salinity_heatmap.pdf", width=5, height=6)
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
          grid.text(sprintf("%.1f", d2.pCO2Salinity.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()








# Group ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d2.Group.matrix <-  paste(signif(moduleTraitCor_Group, 2), "\n(",
                             signif(moduleTraitPvalue_Group, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Output/WGCNA/day2_larvae/heatmaps/Day2_Group_heatmap.png", 1000, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Group,
               xLabels = names(d2.Traits.Group),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d2.Group.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Group"))
dev.off()

# this heatmap looks better
d2.Group.text <-  as.matrix(signif(moduleTraitPvalue_Group, 3))
d2.Group.COR  <-  as.matrix(signif(moduleTraitCor_Group, 3))
pa                  = cluster::pam(d2.Group.COR, k = 3)
col_fun             = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Output/WGCNA/day2_larvae/heatmaps/Day2_Group_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Group, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 2 WGCNA - Group", 
        column_title_gp = gpar(fontsize = 12, fontface = "bold"),
        # row_title = "WGCNA modules",
        #row_km = 4, 
        column_km = 3,
        row_split = paste0("clstr", pa$clustering),
        row_gap = unit(5, "mm"),
        column_gap = unit(5, "mm"),
        # grid.text(matrix(textMatrix)),
        # border = TRUE,
        border = TRUE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", d2.Group.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()




#=====================================================================================
#
# Module eigengene -  MEs boxplots by treatment group
#
#=====================================================================================
MEs_table             <- MEs # new table for plotting 
MEs_table$Sample.Name <- rownames(MEs) # call rows as coolumn to merge with treatment data
MEsPlotting           <- merge(d2.Treatment.data, MEs_table, by = 'Sample.Name') # merge
MEsPlotting_melt      <- melt(MEsPlotting, id=c('Sample.Name', 'Temperature', 'pCO2', 'Salinity', 'pCO2_Salinity', 'All_treatment'))

#plot it - salinity and temperature
png("Output/WGCNA/day2_larvae/Day2_ME_Boxplot_Salinity_Temperature.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Salinity, y=value, fill = factor(Temperature), shape=Temperature)) +
  geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
          stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
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

#plot it - salinity and pCO2
png("Output/WGCNA/day2_larvae/Day2_ME_Boxplot_Salinity_pCO2.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Salinity, y=value, fill = factor(pCO2), shape=pCO2)) +
  geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
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

#plot it - salinity and pCO2
png("Output/WGCNA/day2_larvae/Day2_ME_Boxplot_pCO2_Temperature.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=pCO2, y=value, fill = factor(Temperature), shape=Temperature)) +
  geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
  stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
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

#plot it - all treatment
library(forcats)
png("Output/WGCNA/day2_larvae/Day2_ME_Boxplot_All_Treatment.png", 600, 1000, pointsize=20)
MEsPlotting_melt %>% 
  dplyr::mutate(All_treatment = fct_relevel(All_treatment, 'HHH','LHH','HLH','LLH',
                                                           'HHL','LHL','HLL','LLL')) %>% 
    ggplot(aes(x=All_treatment, y=value, fill = factor(Salinity), shape=Salinity)) +
      geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
      stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
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

# names (colors) of the modules
modNames = substring(names(MEs), 3) # name all the modules, from 3rd character on (first two are ME)

geneModuleMembership = as.data.frame(cor(dds.d2_vst, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


# HHH treatment group
HHH = as.data.frame(as.numeric(d2.Traits.Group$HHH)); # Define variable containing the desired column 
names(HHH) = "HHH"
HHH_geneTraitSignificance = as.data.frame(cor(dds.d2_vst, HHH, use = "p"));
HHH_GSPvalue = as.data.frame(corPvalueStudent(as.matrix(HHH_geneTraitSignificance), nSamples));
names(HHH_geneTraitSignificance) = paste("GS.", names(HHH), sep=""); # MA_geneTraitSignificance - pearsons correlation between reads and the MA grop
names(HHH_GSPvalue) = paste("p.GS.", names(HHH), sep=""); # corPvalueStudent

#  PLOT mean.?mol.CRE.g.protein in the MAGENTA module
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

length(colnames(dds.d2_vst)[moduleColors=="blue"]) # 48 total genes in the violet module

#=====================================================================================
#
#  Call annotation data to get module gene data (prep for downstream GO)
#
#=====================================================================================
annot = read.csv(file = "Data/TagSeq/Seq_details/Seq_Reference_Master.csv",header = T) 
dim(annot) # 66437     8
names(annot) # view the column names to call 
probes = names(as.data.frame(t(d2.data_matrix[, -(1)])))
probes2annot = match(probes, annot$Cvirginica_TranscriptID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot)) # 251
# Should return 0.
#=====================================================================================
#
#  BUILD GENE INFO DATAFRAMES
#
#=====================================================================================
# Create the starting data frame
names(annot)

#   TAOC dataframe - grey60 --------------------------------------------------------------------------- # 
names(HHH_geneTraitSignificance)
names(HHH_GSPvalue)
geneInfo_GROUPS = data.frame(geneSymbol       = annot$Cvirginica_GeneID[probes2annot],
                             TranscriptID     = annot$Cvirginica_TranscriptID[probes2annot],
                             moduleColor      = moduleColors,
                             KEGG_ID          = annot$Cgigas_KEGGID[probes2annot],
                             Protein_name     = annot$Cvirginica_Protein_name[probes2annot],
                             gene_length      = annot$Cvirginica_length[probes2annot],
                             GO.terms         = annot$Annotation_GO_ID[probes2annot],
                             HHH_geneTraitSignificance, HHH_GSPvalue)
# call this specific to the module and trait of interest
View(geneInfo_GROUPS)
modOrder = order(-abs(cor(MEs, HHH, use = "p"))); # order by the strength of the correlation between module and trait values for each sample

for (mod in 1:ncol(geneModuleMembership)) { # Add module membership information in the chosen order

  oldNames = names(geneInfo_GROUPS)
  geneInfo_GROUPS = data.frame(geneInfo_GROUPS, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]]);
  names(geneInfo_GROUPS) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo_GROUPS$moduleColor, -abs(geneInfo_GROUPS$GS.HHH));
geneInfo_GROUPS = geneInfo_GROUPS[geneOrder, ]
View(geneInfo_GROUPS)


#=====================================================================================
#
#  Write csv for the modules and corresponding raw read counts
#
#=====================================================================================
# call the module of interest for follow-up GO analysis 

write.csv(geneInfo_GROUPS, file = "Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv")



#===================================================================================== 
# 
# EXPLORE THE Expression of each module (for loop plots!) BY TREATMENT
#
# evidence from the heatmaps and the regression shwo strong module membership and eigenenge cor strength with lighcyan 
# in response to treatment (primary specifically) and Total Antioxidant capacity 
#
#=====================================================================================
# load necessary data 
load("Output/WGCNA/day2_larvae/day2-networkConstruction-stepByStep.RData") # load dds.d2_vst (in addition to other datasets we do not need...) 

day2_ModuleMembership  <- read.csv(file="Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv", sep=',', header=TRUE) 

d2.Treatment.data <- read.csv(file="Data/TagSeq/day2.exp.data.csv", sep=',', header=TRUE) %>%  
  dplyr::mutate_if(is.character, as.factor) %>% 
  dplyr::rename('Sample.Name' = 'SapleName_readmatrix') %>% 
  dplyr::rename('pCO2' = 'OA') %>% 
  dplyr::select(c('Sample.Name','Temperature','pCO2','Salinity')) %>% 
  dplyr::mutate(All_treatment = paste( (substr(Temperature,1,1)), 
                                       (substr(pCO2,1,1)), 
                                       (substr(Salinity,1,1)), sep = '')) %>% 
  dplyr::mutate(pCO2_Salinity = substr(All_treatment, 2,3)) # experiment treatment data

# cut the samples we do not need from dds.d2 for rlog transofrmation
dim(d2.Treatment.data) # REMEBER! sample tree cut some of these samples out (2 here!) so we need to do the same beofre plotting 
length(rownames(dds.d2_vst)) # 22 total samples


dds.d2_rlogtrans <- as.data.frame(rlogTransformation(assay(dds.d2))) # rlog transoform the expression data matrix (dds object)
dim(dds.d2_rlogtrans) # 4820   24 - note there are 24 samples here, not ommitted from when run WGCNA
dds.d2_rlogtrans <- dds.d2_rlogtrans[,rownames(dds.d2_vst)]
dim(dds.d2_rlogtrans) # 4820   22 - we now have the correct number of samples!
dds.d2_rlogtrans <- tibble::rownames_to_column(dds.d2_rlogtrans,"TranscriptID") # rownames as first column

# write out this rlog tranformed master data (used for plotting below1) 
write.csv(dds.d2_rlogtrans, file = "Output/WGCNA/day2_larvae/d2_rlog_transformed.csv")

# call the module colors 
modcolor <- as.data.frame(unique(day2_ModuleMembership$moduleColor))
names(modcolor)[1] <- "color"


for(i in 1:nrow(modcolor)) {
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs     <- day2_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) %>%  dplyr::select("TranscriptID") %>%  na.omit()
  d2_rlog_Mod      <- dds.d2_rlogtrans %>% dplyr::filter(TranscriptID %in% Mod_geneIDs[,1])
  d2_rlog_Mod_MELT <- melt(d2_rlog_Mod, id=("TranscriptID")) # melt using reshape2
  names(d2_rlog_Mod_MELT)[(2:3)] <-  c('Sample.Name', 'rlog_Expression') # change column names
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d2_rlog_Mod_MELT, d2.Treatment.data, by ='Sample.Name')
  
  # mean Exp response table 
  meanEXp_Mod <- merged_Expdata_Mod %>% 
    select(c('Sample.Name','rlog_Expression','Temperature', 'Salinity', 'pCO2')) %>% 
    group_by(Sample.Name, Temperature, Salinity, pCO2) %>%
    dplyr::summarize(mean.rlogExp = mean(rlog_Expression), 
                     sd.rlogtExp = sd(rlog_Expression),
                     na.rm=TRUE)
  
   
  # summarize datasets further by treatment period  =========================================================================================== #
  # remember:this is a mean of a mean!! First we complete mean vst exp by sample id (compiling all red module genes) - next all sample IDs by the treatment period (below
  # I will use these for mean SE plots 
  
  # Temperature ========================== #
  
  meanEXp_Summary.Temperature <- meanEXp_Mod %>% 
    group_by(Temperature) %>%
    dplyr::summarize(mean = mean(mean.rlogExp), 
                     sd = sd(sd.rlogtExp),
                     n = n(), 
                     se = sd/sqrt(n))
  
  # Salinity treatment ========================== #
  
  meanEXp_Summary.Salinity <- meanEXp_Mod %>% 
    group_by(Salinity) %>%
    dplyr::summarize(mean = mean(mean.rlogExp), 
                     sd = sd(sd.rlogtExp),
                     n = n(), 
                     se = sd/sqrt(n))
  
  # pCO2 treatment ========================== #
  
  meanEXp_Summary.pCO2 <- meanEXp_Mod %>% 
    group_by(pCO2) %>%
    dplyr::summarize(mean = mean(mean.rlogExp), 
                     sd = sd(sd.rlogtExp),
                     n = n(), 
                     se = sd/sqrt(n))
  
  # Salinity treatment ========================== #
  
  meanEXp_Summary.All.Treatment <- meanEXp_Mod %>% 
    group_by(Salinity, Temperature, pCO2) %>%
    dplyr::summarize(mean = mean(mean.rlogExp), 
                     sd = sd(sd.rlogtExp),
                     n = n(), 
                     se = sd/sqrt(n))
  
  # PLOT =========================================================================================== #
  # The errorbars overlapped, so use position_dodge to move them horizontally
  pd <- position_dodge(0.3) # move them .05 to the left and right
   
  # Temperature mean sd plot ========================== #
  
  min_p1 <- min(meanEXp_Summary.Temperature$mean) - max(meanEXp_Summary.Temperature$se)
  max_p1 <- max(meanEXp_Summary.Temperature$mean) + max(meanEXp_Summary.Temperature$se)
  
  Temperature.rlog.Mod <- meanEXp_Summary.Temperature %>% 
    dplyr::mutate(Temperature    = forcats::fct_relevel(Temperature, 'Low', 'High')) %>%
      ggplot(aes(x=Temperature, y=mean, fill=Temperature)) +  # , colour=supp, group=supp))
      theme_classic() +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
      geom_line(position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Temperature") +
      ylab("rlog gene expression") +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("grey", "grey")) +
      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
      # ggtitle(paste("Day 7 WGCNA", modcolor[i,], "Module VST GeneExp", sep =' ')) +
      # expand_limits(y=0) +                                                    # Expand y range
      scale_y_continuous(limits=c((min_p1), (max_p1))) +
      theme(text = element_text(size=10), legend.position="none")
  
  
  # Salinity mean sd plot ========================== #
  
  min_p2 <- min(meanEXp_Summary.Salinity$mean) - max(meanEXp_Summary.Salinity$se)
  max_p2 <- max(meanEXp_Summary.Salinity$mean) + max(meanEXp_Summary.Salinity$se)
  
  Salinity.rlog.Mod <- meanEXp_Summary.Salinity %>% 
    dplyr::mutate(Salinity    = forcats::fct_relevel(Salinity, 'Low', 'High')) %>%
      ggplot(aes(x=Salinity, y=mean, fill=Salinity)) +
      theme_classic() +
      geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
      geom_line(position=pd) +
      geom_point(position=pd, size = 4, shape=21) +            
      xlab("Salinity") +
      ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
      # ylab(paste(modcolor[i,]," Module rlog Gene Expression (Mean +/- SE)", sep = ' ')) +                 # note the mean was first by sample ID THEN by treatment
      scale_fill_manual(values=c("grey", "grey")) +
      # scale_color_manual(values=c("#56B4E9","#E69F00")) +
      # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
      # expand_limits(y=0) +                                                    # Expand y range
      scale_y_continuous(limits=c((min_p2), (max_p2))) +
      theme(text = element_text(size=10), legend.position="none")
  
  
  # pCO2 mean sd plot ========================== #
  
  min_p3 <- min(meanEXp_Summary.pCO2$mean) - max(meanEXp_Summary.pCO2$se)
  max_p3 <- max(meanEXp_Summary.pCO2$mean) + max(meanEXp_Summary.pCO2$se)
  
  pCO2.rlog.Mod <- meanEXp_Summary.pCO2 %>% 
    dplyr::mutate(pCO2    = forcats::fct_relevel(pCO2, 'Low', 'High')) %>%
        ggplot(aes(x=pCO2, y=mean, fill=pCO2)) +
        theme_classic() +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
        geom_line(position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("pCO2") +
        ylab(NULL) +                 # note the mean was first by sample ID THEN by treatment
        # ylab(paste(modcolor[i,]," Module rlog Gene Expression (Mean +/- SE)", sep = ' ')) +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("grey", "grey")) +
        # scale_color_manual(values=c("#56B4E9","#E69F00")) +
        # ggtitle("Day 21 WGCNA red' Module VST GeneExp") +
        # expand_limits(y=0) +                                                    # Expand y range
        scale_y_continuous(limits=c((min_p3), (max_p3))) +
        theme(text = element_text(size=10), legend.position="none")
  
  
  
  
  
  # Assemble these together =========================================================================================== #
  #library(ggpubr)
  single.factor.plot <-  ggarrange(Temperature.rlog.Mod, Salinity.rlog.Mod,  pCO2.rlog.Mod,      
                                    plotlist = NULL,
                                    ncol = 3,
                                    nrow = 1,
                                    labels = NULL)
                          
  
  # Summary plot of all treatments ==================================================================================== #
  # All.Treatment mean sd plot
  min_p4 <- min(meanEXp_Summary.All.Treatment$mean) - max(meanEXp_Summary.All.Treatment$se)
  max_p4 <- max(meanEXp_Summary.All.Treatment$mean) + max(meanEXp_Summary.All.Treatment$se)
  
  AllTreatment.rlog.Mod <- meanEXp_Summary.All.Treatment %>% 
    dplyr::mutate(Salinity    = forcats::fct_relevel(Salinity, 'Low', 'High')) %>%
    dplyr::mutate(pCO2        = forcats::fct_relevel(pCO2, 'Low', 'High')) %>%
    dplyr::mutate(Temperature = forcats::fct_relevel(Temperature, 'Low', 'High')) %>%
        ggplot(aes(x=pCO2, y=mean, group=Temperature)) + # group aesthetic connect line (Slaintiy) and color - the x axis in this case is pCO2 
        theme_classic() +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
        geom_line(aes(linetype=Temperature), position=pd) +
        geom_point(position=pd, size = 4, shape=21) +            
        xlab("pCO2") +
        ylab("rlog gene expression") +                 # note the mean was first by sample ID THEN by treatment
         # ylab(paste(modcolor[i,]," Module rlog Gene Expression (Mean +/- SE)", sep = ' ')) +                 # note the mean was first by sample ID THEN by treatment
        scale_fill_manual(values=c("#56B4E9", "#E69F00")) +
        scale_y_continuous(limits=c((min_p4), (max_p4))) +
        theme(text = element_text(size=15)) + 
        facet_wrap(~Salinity) # facetted by temperature
  
  # output   ======================================================================================================== #
  pdf(paste("Output/WGCNA/day2_larvae/ModuleExpression_Treatment/day2_Exp_Module_",modcolor[i,],".pdf", sep = ''), width=9, height=8)
  print(ggarrange(single.factor.plot, AllTreatment.rlog.Mod,         
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 2,
                  labels = NULL))
  dev.off()
  
}
 
