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

d18.filtered_count_tble  <- read.csv(file="Data/TagSeq/Filtered_counts/filtered_counts_5cpm_50perc/day18.counts.filtered_5cpm50perc.csv", sep=',', header=TRUE) 
d18.filtered_count_mtx   <- data.frame(d18.filtered_count_tble[,-1], row.names=d18.filtered_count_tble[,1]) 


d18.exp_data   <- read.csv(file="Data/TagSeq/day18.exp.data.csv", sep=',', header=TRUE) %>%   mutate_if(is.character, as.factor)


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
                        dplyr::select(c('Sample.Name','Temperature','pCO2','Salinity')) %>% 
                        dplyr::mutate(All_treatment = paste( (substr(Temperature,1,1)), 
                                                             (substr(pCO2,1,1)), 
                                                             (substr(Salinity,1,1)), sep = '')) %>% 
                        dplyr::mutate(pCO2_Salinity = substr(All_treatment, 2,3))

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
png("Output/WGCNA/Day2_larvae/Day2_ClusterTree_Precut.png", 1000, 1000, pointsize=20)
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
png("Output/WGCNA/Day2_larvae/Day2_ClusterTree_Postcut.png", 1000, 1000, pointsize=20)
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
dim(d2.Traits) # 22 Samples 3 columns; now we have 22 samples! - colnames are all treatment, primary and second treatment


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

# pCO2 groups  ===================================================== #
d2.Traits.pCO2          <-  d2.Traits %>% dplyr::select('pCO2')  %>% # primary treatment as Ambient (A) vs. Moderate (M)
                                dplyr::mutate(High = as.factor(as.numeric(pCO2 == "High")))  %>%  # call occurrence of 'A' as 0s and 1s (factor)
                                dplyr::mutate(Low = as.factor(as.numeric(pCO2 == "Low")))    %>%  # call occurrence of 'M'  as 0s and 1s (factor)
                                dplyr::select(-pCO2)
d2.Traits.pCO2  # final dataset of 0,1 for treatment groups - Primary only!

# Group  ================================================================ #
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

# Primary treatment Only
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_ClusterTree_PrimaryTreatment.png", 1000, 1000, pointsize=20)
traitColors_Primary = labels2colors(d7.Traits.Primary); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Primary, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.Primary), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Second treatment 
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_ClusterTree_SecondTreatment.png", 1000, 1000, pointsize=20)
traitColors_Second = labels2colors(d7.Traits.Secondary); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColors_Second, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.Secondary), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Group treatment (primary_second)
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_ClusterTree_GroupTreatment.png", 1000, 1000, pointsize=20)
traitColorsGroup = labels2colors(d7.Traits.Group); # Convert traits to a color representation: white means low, red means high, grey means missing entry
plotDendroAndColors(sampleTree2, traitColorsGroup, # Plot the sample dendrogram and the colors underneath.
                    groupLabels = names(d7.Traits.Group), 
                    main = "Sample dendrogram and trait heatmap")
dev.off()


# write some of te data tus far.... ================================== #
# save(dds.d7_vst, d7.Traits, file = "Presentations/URI/2022_URI_Puritz_Genomics_class/d.7-dataInput.RData")
# write.csv(dds.d7_vst, "Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_vstTransformed_WGCNAdata.csv") # # write the vst transformed data 



# ===================================================================================
#
# soft threshold
# ===================================================================================

dim(dds.d7_vst) #  35 8548 - again double check you have the correct data...
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(dds.d7_vst, powerVector = powers, verbose = 5) #...wait for this to finish
#  pickSoftThreshold 
#  performs the analysis of network topology and aids the
#  user in choosing a proper soft-thresholding power.
#  The user chooses a set of candidate powers (the function provides suitable default values)
#  function returns a set of network indices that should be inspected

# png to output 
sizeGrWindow(9, 5) # set window size 
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_ScaleTopology_SoftThresh.png", 1000, 1000, pointsize=20)
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
softPower = 3 # set your soft threshold based on the plots above 

# unsigned -do not want to run this, remember pros and cons from the presentation; default when unspecified runs as unsigned
adjacency_unsign = adjacency(dds.d7_vst, power = softPower, type="unsigned") # this takes a long time.. just wait...

# signed - must call te type, defaults to unsigned
adjacency_sign = adjacency(dds.d7_vst, power = softPower, type="signed") # this takes a long time.. just wait...


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
# 1    2    3    4    5 
# 3951 2821  862  610  304 

#=====================================================================================
#
#  Step 5: convert numeric network to colors and plot the dendrogram
#
#=====================================================================================

# Convert numeric lables into colors
dynamicColors_sign = labels2colors(dynamicMods_sign) # add colors to module labels (previously numbers)
table(dynamicColors_sign) # lets look at this table...
# dynamicColors_sign
# blue     brown     green turquoise    yellow 
# 2821       862       304      3951       610 
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree_sign, dynamicColors_sign, "Dynamic Tree Cut - SIGNED",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors 'SIGNED'")

#=====================================================================================
#
#  Step 6: Calculate Eigengenes - view thier connectivity based on 'MEDiss = 1-cor(MEs)'
#
#=====================================================================================

# Calculate eigengenes ========================================================== # 

# MEList = moduleEigengenes(dds.d7_vst, colors = dynamicColors)
MEList = moduleEigengenes(dds.d7_vst, colors = dynamicColors_sign)
MEs    = MEList$eigengenes # you can view MEs, condenses gene counts down to a single number for each sample representive of that expressoin pattern 

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average") # Cluster module eigengenes

# Plot the result ================================================================ #
sizeGrWindow(7, 6)
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_ClusterEigengenes.png", 1000, 1000, pointsize=20)
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
save(MEs, moduleLabels, moduleColors, file = "Presentations/URI/2022_URI_Puritz_Genomics_class/Day7-networkConstruction-stepByStep.RData")
# write csv - save the module eigengenes
write.csv(MEs, file = "Presentations/URI/2022_URI_Puritz_Genomics_class/d7.WGCNA_ModulEigengenes.csv")


#=====================================================================================
#
#  Prepare for  module trait associations - Eigengene calc - trait data as numeric
#
#=====================================================================================
# identify modules that are signiFcantly associated with the measured  traits (here as treatment history)

# Since we already have a summary profile (eigengene) for each module,
# we simply correlate eigengenes with external traits and look for the  significant associations:

# Define numbers of genes and samples
nGenes = ncol(dds.d7_vst); # 8548
nSamples = nrow(dds.d7_vst); # 35
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dds.d7_vst, moduleColors)$eigengenes
MEs = orderMEs(MEs0) # reorders the columns (colors/modules)


#=====================================================================================
#
# Module trait correlation
#
#=====================================================================================
# ALL TRAIT DATA

dim(d7.Traits)  # 35  8
dim(MEs)  # 35  6
# moduleTraitCor = cor(MEs, d7.Traits, use = "p");
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# primary 
d7.Traits.Primary.asnum    <- data.frame(lapply(d7.Traits.Primary, function(x) as.numeric(as.character(x))),
                                check.names=F, row.names = row.names(d7.Traits.Primary))
moduleTraitCor_Primary      = cor(MEs, d7.Traits.Primary.asnum, use = "p");
moduleTraitPvalue_Primary   = corPvalueStudent(moduleTraitCor_Primary, nSamples);
moduleTraitPvalue_Primary


# Second
d7.Traits.Secondary.asnum  <- data.frame(lapply(d7.Traits.Secondary, function(x) as.numeric(as.character(x))),
                                       check.names=F, row.names = row.names(d7.Traits.Secondary))
moduleTraitCor_Secondary    = cor(MEs, d7.Traits.Secondary.asnum, use = "p");
moduleTraitPvalue_Secondary = corPvalueStudent(moduleTraitCor_Secondary, nSamples);
moduleTraitPvalue_Secondary


# Group
d7.Traits.Group.asnum      <- data.frame(lapply(d7.Traits.Group, function(x) as.numeric(as.character(x))),
                                        check.names=F, row.names = row.names(d7.Traits.Group))
moduleTraitCor_Group        = cor(MEs, d7.Traits.Group.asnum, use = "p");
moduleTraitPvalue_Group     = corPvalueStudent(moduleTraitCor_Group, nSamples);
moduleTraitPvalue_Group

#=====================================================================================
#
# Heatmaps
#
#=====================================================================================



# PRRIMARY TRETMENT ONLY  ------------------------------------------------------------------ # 

sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.PRIMARYTreatments.matrix <-  paste(signif(moduleTraitCor_Primary, 2), "\n(",
                                      signif(moduleTraitPvalue_Primary, 1), ")", sep = "")
#dim(textMatrix) == dim(moduleTraitCor_treatonly)
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_Treatments_Primary_heatmap2.png", 500, 1000, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Primary,
               xLabels = names(d7.Traits.Primary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.PRIMARYTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Primary Treatment"))
dev.off()

# this heatmap looks better
d7.PRIMARYtreatment.text <-  as.matrix(signif(moduleTraitPvalue_Primary, 3))
pa = cluster::pam(d7.PRIMARYtreatment.text, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
pdf("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_Treatments_Primary_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Primary, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Primary Treatment", 
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
          grid.text(sprintf("%.1f", d7.PRIMARYtreatment.text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()






# SECOND TREATMENT ONLY ------------------------------------------------------------------------------------- # 






sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.SECONDTreatments.matrix <-  paste(signif(moduleTraitCor_Secondary, 2), "\n(",
                                     signif(moduleTraitPvalue_Secondary, 3), ")", sep = "")
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_Treatments_Second_heatmap2.png", 1000, 500, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Secondary,
               xLabels = names(d7.Traits.Secondary),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.SECONDTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Treatment Groups"))
dev.off()

# this heatmap looks better
d7.SECONDtreatment_text<-  as.matrix(signif(moduleTraitPvalue_Secondary, 4))
d7.SECONDtreatment_cor <-  as.matrix(signif(moduleTraitCor_Secondary, 4))
pa = cluster::pam(d7.SECONDtreatment_cor, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

pdf("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_Treatments_Second_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Secondary, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Treatment Groups", 
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
          grid.text(sprintf("%.1f", d7.SECONDtreatment_text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()








# GROUPS ONLY ------------------------------------------------------------------------------------- # 






sizeGrWindow(10,10)
# Will display correlations and their p-values
d7.GROUPTreatments.matrix <-  paste(signif(moduleTraitCor_Group, 2), "\n(",
                                    signif(moduleTraitPvalue_Group, 3), ")", sep = "")
par(mar = c(8, 9.5, 5, 3));
# Display the correlation values within a heatmap plot
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_Group_heatmap2.png", 1000, 500, pointsize=20)
labeledHeatmap(Matrix = moduleTraitCor_Group,
               xLabels = names(d7.Traits.Group),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = TRUE,
               colors = blueWhiteRed(50),
               textMatrix = d7.GROUPTreatments.matrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-0.6,0.6),
               main = paste("Module-trait relationships - Treatment Groups"))
dev.off()

# this heatmap looks better
d7.GROUPtreatment_text<-  as.matrix(signif(moduleTraitPvalue_Group, 4))
d7.GROUPtreatment_cor <-  as.matrix(signif(moduleTraitCor_Group, 4))
pa = cluster::pam(d7.GROUPtreatment_cor, k = 3)
col_fun = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

pdf("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_Group_heatmap.pdf", width=5, height=6)
Heatmap(moduleTraitCor_Group, 
        name = "gene_cor", 
        rect_gp = gpar(col = "grey", lwd = 1),
        column_title = "Day 7 WGCNA - Treatment Groups", 
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
          grid.text(sprintf("%.1f", d7.GROUPtreatment_text[i, j]), x, y, gp = gpar(fontsize = 10))
        })
dev.off()


#=====================================================================================
#
# Module eigengene -  MEs boxplots by treatment group
#
#=====================================================================================
MEs_table             <- MEs # new table for plotting 
MEs_table$Sample.Name <- rownames(MEs) # call rows as coolumn to merge with treatment data
MEsPlotting           <- merge(d7.Treatment.data, MEs_table, by = 'Sample.Name') # merge
MEsPlotting           <- MEsPlotting[,-2] # ommit the all treatments column
MEsPlotting_melt      <- melt(MEsPlotting, id=c('Sample.Name', 'Primary_Treatment', 'Second_Treament'))

#plot it
png("Presentations/URI/2022_URI_Puritz_Genomics_class/Day7_ME_Boxplot.png", 600, 1000, pointsize=20)
ggplot(MEsPlotting_melt, aes(x=Second_Treament, y=value, fill = factor(Primary_Treatment), shape=Primary_Treatment)) +
          geom_boxplot(aes(middle = mean(value)), position=position_dodge(0.8), outlier.size = 0, alpha = 0.5) + 
          stat_summary(fun.y = mean, color = "black", position = position_dodge(0.75),
                       geom = "point", shape = 19, size = 3,
                       show.legend = FALSE) +
          ylab("ModuleEigengene") +
          ylim(-0.5,0.5) +
          scale_fill_manual(values=c("#56B4E9","#D55E00")) +
          geom_hline(yintercept=0, linetype='dotted', col = 'black', size = 1)+
          theme_bw() +
          theme(legend.position = "none") +
          facet_wrap(~variable)
dev.off()

