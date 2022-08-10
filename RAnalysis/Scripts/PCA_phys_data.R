
# LOAD PACKAGES :::::::::::::::::::::::::::::::::::::::::::::::::::::::
library(rlang)
library(dplyr)
library(ggplot2)
library(ggfortify)
library(DESeq2)
library(devtools)
library(ggbiplot)
library(VennDiagram)# venn diagrams
library(eulerr) #venn diagrams -  check out the R shiny app (http://eulerr.co/) 

# SET WORKING DIRECTORY :::::::::::::::::::::::::::::::::::::::::::::::

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/") # personal computer
#setwd("C:/Users/samuel.gurr/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis") # Work computer

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
SampleKey  <- read.csv(file="Data/TagSeq/Seq_details/Sample_Key.csv", sep=',', header=TRUE) %>%  
  dplyr::rename('Sample.Name' = 'SapleName_readmatrix') %>% 
  dplyr::rename('Chamber_tank' = 'ID') %>%   
  dplyr::select(-c('SampleName', 'Treatment'))  %>% 
  dplyr::filter(Age.days %in% 2)





# respiration data ------------------------ #
resp_master         <- read.csv(file="Output/Respiration/RespirationMaster.csv", header=T) # treatment assignments to 'Chamber_Tank'
resp_ref            <- resp_master %>% 
                          select(c('Chamber_tank','Temp','pCO2','Salinity')) %>% 
                          unique()
resp_master_RepMean <- resp_master %>% 
                          dplyr::group_by(Chamber_tank, Date) %>% 
                          dplyr::summarise(mean_resp = mean(resp_ng_L_indiv_hr)) %>% 
                          dplyr::mutate(Day = case_when( # convert Date to Day to merge with theLength and Survival data (below) 
                            Date=="4/30/2021" ~ 1, 
                            Date=="5/7/2021" ~ 8)) %>% 
                          dplyr::select(-Date) #remove date
nrow(resp_master_RepMean) == nrow(resp_ref) # must be true - sanity check before merging
resp_master_RepMean <- merge(resp_master_RepMean, resp_ref) # merged the values 






# length and survival data  ------------ #
length_survival_master  <- read.csv(file="Data/Length_Survival/LengthSurvival_master.csv", header=T) 
length_survival_D1.8 <- length_survival_master %>% 
                              dplyr::rename(Chamber_tank = Id.) %>% 
                              dplyr::select(c('Day', 'Chamber_tank', 
                                              'Survival', 'Average.Length')) %>% 
                              dplyr::filter(Day %in% c(1,8)) # choose only dates that respiration was also measured






# rlog gene expression data ----------------- # 


rlog_WGCNA_D1     <- read.csv(file="Output/WGCNA/day2_larvae/d2_rlog_transformed.csv", header=T) %>% dplyr::select(-X)
d2.Seq_SampleKey  <- read.csv(file="Data/TagSeq/Seq_details/Sample_Key.csv", sep=',', header=TRUE) %>%  
                        dplyr::rename('SampleID' = 'SapleName_readmatrix') %>% 
                        dplyr::rename('Chamber_tank' = 'ID') %>%   
                        dplyr::select(-c('SampleName', 'Treatment', 'Replicate')) %>% 
                        dplyr::filter(Age.days %in% 2)
d2.Aragsat        <- read.csv(file="Data/TagSeq/day2.exp.data.csv", sep=',', header=TRUE) %>% 
                          dplyr::select(c(2,7)) %>% 
                          dplyr::rename('SampleID' = 'SapleName_readmatrix') %>%
                          dplyr::mutate(Aragonite_saturation = case_when(Aragonite_saturation < 0.5 ~ 'Low', 
                                                                         (Aragonite_saturation > 0.5 & Aragonite_saturation < 1.0) ~ 'Mid', 
                                                                         Aragonite_saturation > 1.0 ~ 'High'))
d2.Seq_SampleKey  <- merge(d2.Seq_SampleKey,d2.Aragsat, by = 'SampleID')
ModMem_D1         <- read.csv(file="Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv", header=T) %>%  na.omit()


# data frames and loop sets for the for statement below;
D1_modCols        <- as.data.frame(unique(ModMem_D1$moduleColor)) %>% dplyr::filter(.[[1]] %in% c('black', 'blue', 'brown', 'pink', 'red', 'turquoise')) # yellow and green were NOT significant
meanExp_total     <- data.frame()
meanExp_statsloop <- data.frame(matrix(nrow = 1, ncol = 5)) # create dataframe to save cumunalitively during for loop
colnames(meanExp_statsloop) <- c('Day', 'modColor', 'Gene.count', 'Gene.count.MM<0.5', 'Percent_MM<0.05') # names for comuns in the for loop
meanExp_stats     <- data.frame()

for (i in 1:nrow(D1_modCols)) {
   loopModCol     <- D1_modCols[i,]
   loopModCol_cor <- paste("MM.", loopModCol, sep = '')
   loopModCol_p   <- paste("p.MM.", loopModCol, sep = '')
   
   # all modules per mod color (with significant eigengene-treatment interaction) - no Module Membership threshold
   ModMem         <- ModMem_D1 %>% 
                       dplyr::select(c('TranscriptID',moduleColor, loopModCol_p, loopModCol_cor)) %>% 
                       dplyr::filter(moduleColor %in% loopModCol)
   # all modules per mod color (with significant eigengene-treatment interaction) - Module Membership p < 0.05 based on DEG overalap (view R script)
   ModMem_0.05    <- ModMem %>% 
                       dplyr::filter(.[[3]] < 0.05 & .[[4]] > 0.6) 

   
   MM_0.5_meanExp <- as.data.frame(colMeans(merge(ModMem_0.05, rlog_WGCNA_D1, by = 'TranscriptID')[,-c(1:3)])) %>%  # mean expression by sampleID for this reduced gene pool (Module membership p < 0.05)
                              dplyr::mutate(modcolor = loopModCol) %>% 
                              tibble::rownames_to_column("SampleID") %>% 
                              dplyr::rename(meanExp = 2)
  # print this loop Rbdin for each module 
   meanExp_total <- rbind(meanExp_total,MM_0.5_meanExp) #bind to a cumulative list dataframe
   print(meanExp_total) # print to monitor progress
   
   
   # print stats for each module - these will assist stats for the 0.05 threshold 
   meanExp_statsloop$Day                 <- 1
   meanExp_statsloop$modColor            <- loopModCol
   meanExp_statsloop$Gene.count          <- nrow(ModMem)
   meanExp_statsloop$`Gene.count.MM<0.5` <- nrow(ModMem_0.05)
   meanExp_statsloop$`Percent_MM<0.05`   <- (nrow(ModMem_0.05) / nrow(ModMem)) * 100

   df            <- data.frame(meanExp_statsloop) # name dataframe for this single row
   meanExp_stats <- rbind(meanExp_stats,df) #bind to a cumulative list dataframe
   print(meanExp_stats) # print to monitor progress
   
}
meanExp_total_wide <- reshape2::dcast(meanExp_total, SampleID ~ modcolor, value.var="meanExp") # convert to wide format to merge for OCA analysis with physiological variables 
meanExp_Master     <- merge(d2.Seq_SampleKey, meanExp_total_wide) %>% select(-c(Temperature, OA, Salinity, Age.days))


meanExp_stats # percent contribution of Module member threshold cutoff to the total module membership 
colMeans(meanExp_stats[c(4:5)])
meanExp_stats %>% summarise(sd_Gene_count = sd(meanExp_stats$Gene.count.MM.0.5),
                            sd_Perc = sd(meanExp_stats$Percent_MM.0.05))
# > 0.8 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
# 70.50000          12.01422

# > 0.6 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
# 240.00000          42.05936 + - 5.293175

# > 0.4 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
# 417.00000          73.28806






# master phys (merge resp and shell length + survival data )
Master_Days1.8_phys  <- merge(resp_master_RepMean,length_survival_D1.8) %>%  
                            dplyr::mutate(AllTreat = paste(Temp, pCO2, Salinity, sep = ''))
Master_Day1          <- Master_Days1.8_phys %>%  dplyr::filter(Day %in% 1)
Master_Day8          <- Master_Days1.8_phys %>%  dplyr::filter(Day %in% 8)

Master_Days1.8_phys_cors <- Master_Days1.8_phys %>% # Master_Days1.8_phys_cors = corrected for the treatment day - to run PCA regardless of time point!
                      dplyr::mutate(Length_Cor = ifelse( (Day == 1) | (Day == 8), 
                                                         Average.Length/mean(Master_Day1$Average.Length), 
                                                         Average.Length/mean(Master_Day8$Average.Length))) %>% 
                      dplyr::mutate(Surv_Cor = ifelse( (Day == 1) | (Day == 8), 
                                                       Survival/mean(Master_Day1$Survival), 
                                                       Survival/mean(Master_Day8$Survival))) %>%  
                      dplyr::mutate(Resp_Cor = ifelse( (Day == 1) | (Day == 8), 
                                                       mean_resp/mean(Master_Day1$mean_resp), 
                                                       mean_resp/mean(Master_Day8$mean_resp))) 


# Run a PCA for Day 1 
Master_Day1$Day <- as.factor(Master_Day1$Day) # convert Day into a factor
Day1PCA         <- (merge(Master_Day1, meanExp_Master))[-9,]

# PCA Day 1 
phys_pca1   <- prcomp(Day1PCA[,c(3,7,8,12:17)], # all numeric (phys + all modules) - PCA 1 = 0.4133 , PCA 2 0.1786  (cumulative 0.5919)
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(3,7,8)],   # phys only  - PCA 1 = 0.5805  PCA 2 0.2946   (cumulative 0.8751 )
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(11:16)],   # modules only   (cumulative 0.6990 )
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(3,7,8,12,13,17)],   # main effect modules only (brown blur, turquoise)   (cumulative 96.28)
                      center = TRUE,
                      scale. = TRUE)

print(phys_pca1)

summary(phys_pca1)

# > 0.8 Pearson's cor and > 0.05 P value 
# Cumulative Proportion PC1+PC2 == 0.6097 (0.6392 without row 9 outlier) < 1% to 0.6 cor coeff

# > 0.6 Pearson's cor and > 0.05 P value 
# Cumulative Proportion PC1+PC2 == 0.6032 (0.6373 without row 9 outlier) > 1% to 0.4 cor coeff; < 1% to 0.8 cor coeff

# > 0.4 Pearson's cor and > 0.05 P value 
# Cumulative Proportion PC1+PC2 == 0.5919 (0.6268 without row 9 outlier) > 1% to 0.6 cor coeff

# summary: the histogram of cor coeffs for DEGs overlapped with WGCNA module indicates 0.6 as a good cut off
# we confirm here that 0.6 cor coeff and < 0.05 p value can increase explanatory variance relative to 0.4 cor coeff
# however more conservative threshold (0.8) does not make a difference, albeit these (cor coeff cutoffs) change the number of genes contr'ing to the mean by >2 fold!
p1pCO2 <- ggbiplot(phys_pca1,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = as.factor(Day1PCA$pCO2),
                  ellipse = TRUE,
                  circle = TRUE,
                  ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +   ggtitle("day1, pCO2") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

p1Sal <- ggbiplot(phys_pca1,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = Day1PCA$Salinity,
                  ellipse = TRUE,
                  circle = TRUE,
                  ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, salinity") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

p1Temp <- ggbiplot(phys_pca1,
                   obs.scale = 1,
                   var.scale = 1,
                   groups = Day1PCA$Temp,
                   ellipse = TRUE,
                   circle = TRUE,
                   ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, temp") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

p1Arag <- ggbiplot(phys_pca1,
                   obs.scale = 1,
                   var.scale = 1,
                   groups = Day1PCA$Aragonite_saturation.y,
                   ellipse = TRUE,
                   circle = TRUE,
                   ellipse.prob = 0.5) +
  scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, aragonite saturation") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

library(ggpubr)
ggarrange(p1Sal, p1pCO2, p1Temp, ncol = 3,nrow = 1)


pdf("Output/PCA_day1_phys.pdf", 
    width = 10, height = 10)
print(ggarrange(p1Sal, p1pCO2, p1Temp ))
dev.off()



















# Run a PCA for day 8

Master_Day8 <- Master_Days1.8_phys %>% dplyr::filter(Day %in% 8)
Master_Day8$Day <- as.factor(Master_Day8$Day)
phys_pca8 <- prcomp(Master_Day8[,c(3,7,8)], #all numeric data for phys values
                             center = TRUE,
                             scale. = TRUE)
print(phys_pca8)
#                   PC1        PC2        PC3
# mean_resp      -0.4912969 0.80057018 -0.3430959
# Survival        0.6565117 0.08151057 -0.7498989
# Average.Length  0.5723808 0.59366945  0.5656296
# PC1 decreases mean resp which increasing the others 
# PC2 increase mean resp whole not changing survival
summary(phys_pca8)
# Importance of components:
#                         PC1    PC2    PC3
# Standard deviation     1.3332 0.9035 0.6374
# Proportion of Variance 0.5925 0.2721 0.1354
# Cumulative Proportion  0.5925 0.8646 1.0000
p8Sal <- ggbiplot(phys_pca8,
              obs.scale = 1,
              var.scale = 1,
              groups = Master_Day8$Salinity,
              ellipse = TRUE,
              circle = TRUE,
              ellipse.prob = 0.68) +
      scale_color_discrete(name = '') + 
      theme_classic() + 
  ggtitle("day8, salinity") +
  theme(legend.direction = 'horizontal',
                   legend.position = 'top')

p8pCO2 <- ggbiplot(phys_pca8,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = Master_Day8$pCO2,
                  ellipse = TRUE,
                  circle = TRUE,
                  ellipse.prob = 0.68) +
  scale_color_discrete(name = '') + 
  theme_classic() + 
  ggtitle("day8, pCO2") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

p8Temp <- ggbiplot(phys_pca8,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = Master_Day8$Temp,
                  ellipse = TRUE,
                  circle = TRUE,
                  ellipse.prob = 0.68) +
  scale_color_discrete(name = '') + 
  theme_classic() + 
  ggtitle("day8, temperature") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

library(ggpubr)
ggarrange(p8Sal, p8pCO2, p8Temp )
# arrows together indicate a high correlation
# exaple - reapiration and the average survival and length are poorly correlated 

# PC1 is weakly correlated with respiration
# PC2 is highly positively correlated with mean resp
# PC1 explains most of the variation
pdf("Output/PCA_day8_phys.pdf", 
    width = 10, height = 10)
print(ggarrange(p8Sal, p8pCO2, p8Temp ))
dev.off()












# Run a PCA for all length and survival data 



length_survival_all <- length_survival_master %>% 
  dplyr::rename(Chamber_tank = Id.) %>% 
  dplyr::select(c('Day', 'Chamber_tank', 
                  'Survival', 'Average.Length'))
length_survival_all  <- merge(length_survival_all, SampleKey) %>% dplyr::select(-Age.days)
length_survival_D1   <- length_survival_all %>% dplyr::filter(Day == 1)
length_survival_D4   <- length_survival_all %>% dplyr::filter(Day == 4)
length_survival_D8   <- length_survival_all %>% dplyr::filter(Day == 8)
length_survival_D11  <- length_survival_all %>% dplyr::filter(Day == 11)
length_survival_D15  <- length_survival_all %>% dplyr::filter(Day == 15)



LengthSurvival_Master <- length_survival_all %>% 
  dplyr::mutate(Length_Cor = case_when(Day == '1' ~ Average.Length/mean(length_survival_D1$Average.Length, na.rm = TRUE),
                                       Day == '4'  ~ Average.Length/mean(length_survival_D4$Average.Length, na.rm = TRUE),
                                       Day == '8'  ~ Average.Length/mean(length_survival_D8$Average.Length, na.rm = TRUE),
                                       Day == '11'  ~ Average.Length/mean(length_survival_D11$Average.Length, na.rm = TRUE),
                                       Day == '15'  ~ Average.Length/mean(length_survival_D15$Average.Length, na.rm = TRUE))) %>% 
  dplyr::mutate(Survival_Cor = case_when(Day == '1' ~ Survival/mean(length_survival_D1$Survival, na.rm = TRUE),
                                     Day == '4'  ~ Survival/mean(length_survival_D4$Survival, na.rm = TRUE),
                                     Day == '8'  ~ Survival/mean(length_survival_D8$Survival, na.rm = TRUE),
                                     Day == '11'  ~ Survival/mean(length_survival_D11$Survival, na.rm = TRUE),
                                     Day == '15'  ~ Survival/mean(length_survival_D15$Survival, na.rm = TRUE)))             
                
LengthSurvival_Master_OM <- LengthSurvival_Master %>% na.omit()

BiocManager::install('GEOquery')          
library(PCAtools)
library(Biobase)
library(GEOquery)



metadata_1     <- LengthSurvival_Master_OM[,c(1:2,6:8,10:11)] %>% 
                                         dplyr::mutate(ID = paste(Chamber_tank, Day, sep = '_'))
metadata_1$Chamber_tank <- as.character(metadata_1$Chamber_tank)
metadata_2     <- rbind((metadata_1 %>% 
                       # dplyr::select(-Survival_Cor) %>% 
                        dplyr::mutate(X = paste('L', Chamber_tank, sep = '_'))), 
                     (metadata_1 %>% 
                       # dplyr::select(-Length_Cor) %>% 
                       dplyr::mutate(X = paste('S', Chamber_tank, sep = '_'))) )
metadata_2MTX  <- as.matrix(metadata_1)

rownames(metadata_2MTX) <- metadata_1$Chamber_tank
metadata_3         <- as.data.frame(metadata_2MTX)
mat                <- as.data.frame(t(metadata_3 %>% dplyr::select(Length_Cor,Survival_Cor)))
metadata_final     <- metadata_3 %>% dplyr::select(-c(ID, Chamber_tank ,Length_Cor, Survival_Cor))
# filter the expression data to match the samples in our pdata
mat <- as.matrix(mat[,which(colnames(mat) %in% rownames(metadata_final))])
# check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(metadata_final))
str(mat)
metadata_final <- metadata_final %>% mutate_if(is.character,as.factor)

p <- pca(mat, metadata_final, removeVar = 0.1)

biplot(p,
       colby = 'Salinity', colkey = c('High' = 'forestgreen', 'Low' = 'purple'),
       # ellipse config
       ellipse = TRUE,
       ellipseType = 't',
       ellipseLevel = 0.95,
       ellipseFill = TRUE,
       ellipseAlpha = 1/4,
       ellipseLineSize = 1.0,
       xlim = c(-125,125), ylim = c(-50, 80),
       hline = 0, vline = c(-25, 0, 25),
       legendPosition = 'top', legendLabSize = 16, legendIconSize = 8.0)


p$rotated
phys_all_PCA   <- prcomp(LengthSurvival_Master_OM[,c(10:11)], #all numeric data for phys values
                      center = TRUE,
                      scale. = TRUE)
print(phys_all_PCA)

# plot the pca
pCA_all_time <- ggbiplot(phys_all_PCA,
                   obs.scale = 1,
                   var.scale = 1,
                   groups = as.factor(LengthSurvival_Master_OM$Salinity),
                   ellipse = TRUE,
                   circle = TRUE,
                   ellipse.prob = 0.61) +
  scale_color_discrete(name = '') + 
  theme_classic() + 
  ggtitle("all, time") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')
pCA_all_time






























mat <- exprs(gset[[1]])

# remove Affymetrix control probes
mat <- mat[-grep('^AFFX', rownames(mat)),]
str(mat)
# extract information of interest from the phenotype data (pdata)
idx <- which(colnames(pData(gset[[1]])) %in%
               c('relation', 'age:ch1', 'distant rfs:ch1', 'er:ch1',
                 'ggi:ch1', 'grade:ch1', 'size:ch1',
                 'time rfs:ch1'))
metadata <- data.frame(pData(gset[[1]])[,idx],
                       row.names = rownames(pData(gset[[1]])))

# tidy column names
colnames(metadata) <- c('Study', 'Age', 'Distant.RFS', 'ER', 'GGI', 'Grade',
                        'Size', 'Time.RFS')

# prepare certain phenotypes of interest
metadata$Study <- gsub('Reanalyzed by: ', '', as.character(metadata$Study))
metadata$Age <- as.numeric(gsub('^KJ', NA, as.character(metadata$Age)))
metadata$Distant.RFS <- factor(metadata$Distant.RFS,
                               levels = c(0,1))
metadata$ER <- factor(gsub('\\?', NA, as.character(metadata$ER)),
                      levels = c(0,1))
metadata$ER <- factor(ifelse(metadata$ER == 1, 'ER+', 'ER-'),
                      levels = c('ER-', 'ER+'))
metadata$GGI <- as.numeric(as.character(metadata$GGI))
metadata$Grade <- factor(gsub('\\?', NA, as.character(metadata$Grade)),
                         levels = c(1,2,3))
metadata$Grade <- gsub(1, 'Grade 1', gsub(2, 'Grade 2', gsub(3, 'Grade 3', metadata$Grade)))
metadata$Grade <- factor(metadata$Grade, levels = c('Grade 1', 'Grade 2', 'Grade 3'))
metadata$Size <- as.numeric(as.character(metadata$Size))
metadata$Time.RFS <- as.numeric(gsub('^KJX|^KJ', NA, metadata$Time.RFS))

# remove samples from the pdata that have any NA value
discard <- apply(metadata, 1, function(x) any(is.na(x)))
metadata <- metadata[!discard,]

# filter the expression data to match the samples in our pdata
mat <- mat[,which(colnames(mat) %in% rownames(metadata))]

# check that sample names match exactly between pdata and expression data 
all(colnames(mat) == rownames(metadata))
## [1] TRUE
str(metadata)
typeof(metadata$Age)
  
p <- pca(mat, metadata = metadata, removeVar = 0.1)
