
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
resp_master         <- read.csv(file="Output/Respiration/RespirationMaster.csv", header=T) %>% # treatment assignments to 'Chamber_Tank'
                          dplyr::rename(Sample.ID = Chamber_tank)


resp_ref            <- resp_master %>% 
                          select(c('Sample.ID',
                                   'Temp',
                                   'pCO2',
                                   'Salinity')) %>% 
                          unique()

resp_master_RepMean <- resp_master %>% 
                          dplyr::group_by(Sample.ID, Date) %>% 
                          dplyr::summarise(mean_resp = mean(resp_ng_L_indiv_hr)) %>% 
                          dplyr::mutate(Age = case_when( # convert Date to Day to merge with theLength and Survival data (below) 
                            Date=="4/30/2021" ~ 1, 
                            Date=="5/7/2021" ~ 8)) %>% 
                          dplyr::select(-Date) #remove date


nrow(resp_master_RepMean) == nrow(resp_ref) # must be true - sanity check before merging
resp_master_RepMean <- merge(resp_master_RepMean, resp_ref) # merged the values 






# length and survival data  ------------ #

# survival - days 1 and 18 (no survival data for day 2@)
survival_master_d1_d18  <- read.csv(file="Data/Survival/Survival_master.csv", header=T) %>% 
                              dplyr::rename(Chamber_tank = Id.) %>% 
                              dplyr::select(c('Day', 
                                              'Chamber_tank', 
                                              'Survival')) %>% 
                              dplyr::rename(Sample.ID = Chamber_tank) %>% 
                              dplyr::rename(Age = Day) %>% 
                              dplyr::filter(Age %in% c(1,18)) %>%  # choose only dates that respiration was also measured
                              dplyr::mutate(Age = case_when(Age %in% 18 ~ 22, # call the day 18 survival data day 22 to allow merge (note! this is NOT day 22 survival data!)
                                                            TRUE ~ as.numeric(Age)))


survival_master_d1_d18



# length data for day 1 (larvae) and day 22 (spat)
# upload
length_master_d1_d22    <- read.csv(file="Data/Length/cumulative_raw/Length_cumulative_raw.csv", header=T) %>% 
                              dplyr::select(c('Age', 
                                              'Sample.ID', 
                                              'length_um',
                                              'stage')) %>% 
                              dplyr::filter(Age %in% c(1,22))
# mean for day 1 larvae
length_master_d1MEAN    <- length_master_d1_d22 %>% 
                              dplyr::select(-stage) %>%  # all are larvae
                              dplyr::filter(Age %in% 1) %>% 
                              dplyr::group_by(Age, Sample.ID) %>% 
                              dplyr::summarise(meanLength_um  = mean(length_um))
# mean for day 22 spat
length_master_d22MEAN    <- length_master_d1_d22 %>% 
                              dplyr::filter(stage %in% 'spat') %>% 
                              dplyr::select(-stage) %>%  # do not need anymore..
                              dplyr::filter(Age %in% 22) %>% 
                              dplyr::group_by(Age, Sample.ID) %>% 
                              dplyr::summarise(meanLength_um  = mean(length_um))
# merge these
length_master_d1_d22     <- rbind(as.data.frame(length_master_d1MEAN), as.data.frame(length_master_d22MEAN))



# FINAL STEP merge the survival and length data 
Master_length_survival_d1_d22 <- merge(length_master_d1_d22, survival_master_d1_d18, by = c('Age', 'Sample.ID'))

Master_length_survival_d1     <- Master_length_survival_d1_d22 %>% dplyr::filter(Age %in% 1)
Master_length_survival_d22     <- Master_length_survival_d1_d22 %>% dplyr::filter(Age %in% 22)
                    




#  DAY 1 PCA --------------------------------------------------------------------------- # 



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
ModMem_D1         <- read.csv(file="Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv", header=T) #%>%  na.omit()


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

   
   MM_0.5_meanExp <- as.data.frame(colMeans(merge(ModMem_0.05, rlog_WGCNA_D1, by = 'TranscriptID')[,-c(1:4)])) %>%  # mean expression by sampleID for this reduced gene pool (Module membership p < 0.05)
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
   # Day  modColor Gene.count Gene.count.MM.0.5 Percent_MM.0.05
   # 1     black        431               206        47.79582
   # 1      blue        770               330        42.85714
   # 1     brown        690               342        49.56522
   # 1      pink        379               145        38.25858
   # 1       red        506               199        39.32806
   # 1 turquoise        881               312        35.41430
   
}
meanExp_total_wide <- reshape2::dcast(meanExp_total, SampleID ~ modcolor, value.var="meanExp") # convert to wide format to merge for OCA analysis with physiological variables 
meanExp_Master     <- merge(d2.Seq_SampleKey, meanExp_total_wide) %>% 
  select(-c(Temperature, OA, Salinity, Age.days)) %>% 
  dplyr::rename(Sample.ID = Chamber_tank)


meanExp_stats # percent contribution of Module member threshold cutoff to the total module membership 
colMeans(meanExp_stats[c(4:5)]) # 42.05936 +- 5.293175
meanExp_stats %>% summarise(sd_Gene_count = sd(meanExp_stats$Gene.count.MM.0.5),
                            sd_Perc = sd(meanExp_stats$Percent_MM.0.05))
# > 0.8 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
# 70.50000          12.01422

# > 0.6 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
# 240.00000  +- 79.84485        42.05936 + - 5.293175  *** this is the one we are using

# > 0.4 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
# 417.00000          73.28806






# master phys (merge resp and shell length + survival data )
Master_Days1.8_phys  <- merge(resp_master_RepMean,  Master_length_survival_d1) %>%  # contains resp for day 8 and redundant length andsurvival copied from day 1 to 8 (do not)
                            dplyr::mutate(AllTreat = paste(Temp, pCO2, Salinity, sep = ''))
Master_Day1          <- Master_Days1.8_phys %>%  dplyr::filter(Age %in% 1)
# Master_Day8          <- Master_Days1.8_phys %>%  dplyr::filter(Day %in% 8) # run this if you want a PCA for just the physiology data on day 8 - NO RNA this date!

Master_Days1_phys_cors <- Master_Day1 %>% # Master_Days1.8_phys_cors = corrected for the treatment day - to run PCA regardless of time point!
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
Master_Day1$Age <- as.factor(Master_Day1$Age) # convert Day into a factor
Day1PCA         <- (merge(Master_Day1, meanExp_Master)) %>% #[-9,] %>%  # outlier omit, chose not to
                      dplyr::mutate(Aragonite_saturation = 
                      case_when(Aragonite_saturation == 'Low' ~ 'Low', 
                               (Aragonite_saturation == 'Mid' & Salinity == 'L') ~ 'Mid_Sal', 
                               (Aragonite_saturation == 'Mid' & pCO2     == 'H') ~ 'Mid_pCO2', 
                                Aragonite_saturation == 'High' ~ 'High'))
                    

# PCA Day 1 
phys_pca1   <- prcomp(Day1PCA[,c(3,7,8,12:17)], # all numeric (phys + all modules) - PCA 1 = 0.4298  , PCA 2 0.1810   (cumulative 0.6108 )
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(3,7,8)],   # phys only  - PCA 1 = 0.5849   PCA 2 0.2895   (cumulative 0.8744 )
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(12:17)],   # modules only   (cumulative 0.7143  )
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(3,7,8,12,13,17)],   # main effect modules only (brown blur, turquoise)   (cumulative 96.28)
                      center = TRUE,
                      scale. = TRUE)
phys_pca1   <- prcomp(Day1PCA[,c(7,8,12:17)],   # without resp 0.452 0.2012  0.6532 
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
                   groups = Day1PCA$Aragonite_saturation,
                   ellipse = TRUE,
                   circle = TRUE,
                   ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, aragonite saturation") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

library(ggpubr)
ggarrange(p1Sal, p1pCO2, p1Temp, p1Arag, ncol = 2,nrow = 2)


pdf("Output/PCA_day1_phys.pdf", 
    width = 10, height = 10)
print(ggarrange(p1Sal, p1pCO2, p1Temp,p1Arag ))
dev.off()



































#  DAY 22 PCA --------------------------------------------------------------------------- # 



# rlog gene expression data ----------------- # 
rlog_WGCNA_D22     <- read.csv(file="Output/WGCNA/day18_spat/d18_rlog_transformed.csv", header=T) %>% dplyr::select(-X)
d22.Seq_SampleKey    <- read.csv(file="Data/TagSeq/Seq_details/Sample_Key.csv", sep=',', header=TRUE) %>%  
  dplyr::rename('SampleID' = 'SapleName_readmatrix') %>% 
  dplyr::rename('Chamber_tank' = 'ID') %>%   
  dplyr::select(-c('SampleName', 'Treatment', 'Replicate')) %>% 
  dplyr::filter(Age.days %in% 18)
d22.Aragsat        <- read.csv(file="Data/TagSeq/day18.exp.data.csv", sep=',', header=TRUE) %>% 
  dplyr::select(c(2,7)) %>% 
  dplyr::rename('SampleID' = 'SapleName_readmatrix') %>%
  dplyr::mutate(Aragonite_saturation = case_when(Aragonite_saturation < 0.5 ~ 'Low', 
                                                 (Aragonite_saturation > 0.5 & Aragonite_saturation < 1.0) ~ 'Mid', 
                                                 Aragonite_saturation > 1.0 ~ 'High'))
d22.Seq_SampleKey  <- merge(d22.Seq_SampleKey,d22.Aragsat, by = 'SampleID')
ModMem_D22         <- read.csv(file="Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv", header=T)
nrow(ModMem_D22 %>% filter(moduleColor %in% 'turquoise'))

# data frames and loop sets for the for statement below;
D22_modCols            <- as.data.frame(unique(ModMem_D22$moduleColor)) %>% dplyr::filter(.[[1]] %in% c('blue', 'red', 'salmon', 'tan', 'green', 'turquoise')) # yellow and green were NOT significant
D22_meanExp_total      <- data.frame()
D22_meanExp_statsloop  <- data.frame(matrix(nrow = 1, ncol = 5)) # create dataframe to save cumunalitively during for loop
colnames(D22_meanExp_statsloop) <- c('Day', 'modColor', 'Gene.count', 'Gene.count.MM<0.5', 'Percent_MM<0.05') # names for comuns in the for loop
D22_meanExp_stats      <- data.frame()

for (i in 1:nrow(D22_modCols)) {
  loopModCol     <- D22_modCols[i,]
  loopModCol_cor <- paste("MM.", loopModCol, sep = '')
  loopModCol_p   <- paste("p.MM.", loopModCol, sep = '')
  
  # all modules per mod color (with significant eigengene-treatment interaction) - no Module Membership threshold
  ModMem         <- ModMem_D22 %>% 
    dplyr::select(c('TranscriptID',moduleColor, loopModCol_p, loopModCol_cor)) %>% 
    dplyr::filter(moduleColor %in% loopModCol)
  # all modules per mod color (with significant eigengene-treatment interaction) - Module Membership p < 0.05 based on DEG overalap (view R script)
  ModMem_0.05    <- ModMem %>% 
    dplyr::filter(.[[3]] < 0.05 & .[[4]] > 0.6) 
  
  ModMem_outliers    <- ModMem %>% # exists non significant with high correlation
    dplyr::filter(.[[3]] > 0.05 & .[[4]] > 0.6) 
  
  ModMem_outliers    <- ModMem %>% 
    dplyr::filter(.[[3]] < 0.05 & .[[4]] < 0.6) 
  
  MM_0.5_meanExp <- as.data.frame(colMeans(merge(ModMem_0.05, rlog_WGCNA_D22, by = 'TranscriptID')[,-c(1:4)])) %>%  # mean expression by sampleID for this reduced gene pool (Module membership p < 0.05)
    dplyr::mutate(modcolor = loopModCol) %>% 
    tibble::rownames_to_column("SampleID") %>% 
    dplyr::rename(meanExp = 2)
  # print this loop Rbdin for each module 
  D22_meanExp_total <- rbind(D22_meanExp_total,MM_0.5_meanExp) #bind to a cumulative list dataframe
  print(D22_meanExp_total) # print to monitor progress
  
  
  # print stats for each module - these will assist stats for the 0.05 threshold 
  D22_meanExp_statsloop$Day                 <- 1
  D22_meanExp_statsloop$modColor            <- loopModCol
  D22_meanExp_statsloop$Gene.count          <- nrow(ModMem)
  D22_meanExp_statsloop$`Gene.count.MM<0.5` <- nrow(ModMem_0.05)
  D22_meanExp_statsloop$`Percent_MM<0.05`   <- (nrow(ModMem_0.05) / nrow(ModMem)) * 100
  
  df            <- data.frame(D22_meanExp_statsloop) # name dataframe for this single row
  D22_meanExp_stats <- rbind(D22_meanExp_stats,df) #bind to a cumulative list dataframe
  print(D22_meanExp_stats) # print to monitor progress
  
}
D22_meanExp_total_wide <- reshape2::dcast(D22_meanExp_total, SampleID ~ modcolor, value.var="meanExp") # convert to wide format to merge for OCA analysis with physiological variables 
D22_meanExp_Master     <- merge(d22.Seq_SampleKey, D22_meanExp_total_wide) %>% 
  dplyr::mutate(pCO2 = substr(OA, 1,1)) %>% 
  dplyr::mutate(Salinity = substr(Salinity,1,1)) %>% 
  dplyr::mutate(Temp = substr(Temperature,1,1)) %>% 
  select(-c(Temperature, OA, Age.days)) %>% 
  dplyr::rename(Sample.ID = Chamber_tank)


D22_meanExp_stats # percent contribution of Module member threshold cutoff to the total module membership 
colMeans(D22_meanExp_stats[c(4:5)]) # 364.50000          70.95545 + - 4.791734
D22_meanExp_stats %>% summarise(sd_Gene_count = sd(D22_meanExp_stats$Gene.count.MM.0.5),
                                sd_Perc = sd(D22_meanExp_stats$Percent_MM.0.05))
# get full mean of d1 and d22
colMeans(rbind((meanExp_stats[c(4:5)]),(D22_meanExp_stats[c(4:5)]))) # mean 317.83333   perc of full 56.17219           70.95545 + - 4.791734
rbind((meanExp_stats[c(4:5)]),(D22_meanExp_stats[c(4:5)])) %>% # sd dev. 178.6642 percent of full 15.50647
    dplyr::summarise(sd_Gene_count = sd(Gene.count.MM.0.5),
                                sd_Perc = sd(Percent_MM.0.05))
# > 0.6 Pearson's cor and > 0.05 P value 
# Gene.count.MM.0.5   Percent_MM.0.05 
#  364.50000 +-   206.5224    70.95545 +-  +- 4.496341 *** this is the one we are using








# Run a PCA for Day 1 
Master_length_survival_d22$Age <- as.factor(Master_length_survival_d22$Age) # convert Day into a factor
Day22PCA         <- (merge(Master_length_survival_d22, D22_meanExp_Master)) %>% #[-9,] %>%  # outlier omit, chose not to
  dplyr::mutate(Aragonite_saturation = 
                  case_when(Aragonite_saturation == 'Low' ~ 'Low', 
                            (Aragonite_saturation == 'Mid' & Salinity == 'L') ~ 'Mid_Sal', 
                            (Aragonite_saturation == 'Mid' & pCO2     == 'H') ~ 'Mid_pCO2', 
                            Aragonite_saturation == 'High' ~ 'High'))


# PCA Day 1 
phys_pca22   <- prcomp(Day22PCA[,c(3,4,8:13)], # all numeric (phys + all modules) - PCA 1 = 0.5477    , PCA 2 0.2290     (cumulative 0.7767   )
                      center = TRUE,
                      scale. = TRUE)
print(phys_pca22)

summary(phys_pca22)


# > 0.6 Pearson's cor and > 0.05 P value 
# Cumulative Proportion PC1+PC2 == 0.7767  

# summary: the histogram of cor coeffs for DEGs overlapped with WGCNA module indicates 0.6 as a good cut off
# we confirm here that 0.6 cor coeff and < 0.05 p value can increase explanatory variance relative to 0.4 cor coeff
# however more conservative threshold (0.8) does not make a difference, albeit these (cor coeff cutoffs) change the number of genes contr'ing to the mean by >2 fold!
p22pCO2 <- ggbiplot(phys_pca22,
                   obs.scale = 1,
                   var.scale = 1,
                   groups = as.factor(Day22PCA$pCO2),
                   ellipse = TRUE,
                   circle = TRUE,
                   ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +   ggtitle("day1, pCO2") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

p22Sal <- ggbiplot(phys_pca22,
                  obs.scale = 1,
                  var.scale = 1,
                  groups = Day22PCA$Salinity,
                  ellipse = TRUE,
                  circle = TRUE,
                  ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, salinity") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

# p22Temp <- ggbiplot(phys_pca22,  # all high temperature!!!!!
#                    obs.scale = 1,
#                    var.scale = 1,
#                    groups = Day22PCA$Temp,
#                    ellipse = TRUE,
#                    circle = TRUE,
#                    ellipse.prob = 0.67) +
#   scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, temp") +
#   theme(legend.direction = 'horizontal',
#         legend.position = 'top')

p22Arag <- ggbiplot(phys_pca22,
                   obs.scale = 1,
                   var.scale = 1,
                   groups = Day22PCA$Aragonite_saturation,
                   ellipse = TRUE,
                   circle = TRUE,
                   ellipse.prob = 0.67) +
  scale_color_discrete(name = '') +  theme_classic() +  ggtitle("day1, aragonite saturation") +
  theme(legend.direction = 'horizontal',
        legend.position = 'top')

library(ggpubr)
ggarrange(p22Sal, p22pCO2, p22Arag, ncol = 2,nrow = 2)


pdf("Output/PCA_day22_phys.pdf", 
    width = 10, height = 10)
print(ggarrange(p22Sal, p22pCO2, p22Arag))
dev.off()




