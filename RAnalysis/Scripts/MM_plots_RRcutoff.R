---
# title: "Module Membershop plots RRcutoff"
# author: "Samuel Gurr"
# date: "September 12, 2022"
---
  
# LOAD PACKAGES
library(dplyr)
library(zoo)
library(DESeq2)
library(reshape2)
library(ggplot2)
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")


##################################################################################### #
###################3  Day 2 RRcutoff Mod Mem plots ################################## #
##################################################################################### #

# Load data 
day2_ModuleMembership  <- read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership_RRcutoff.csv") 
d2_rlog                <- read.csv(file="Output/WGCNA/day2_larvae/d2_rlog_transformed.csv", sep=',', header=TRUE) %>% dplyr::select(-X)


d2.Treatment.data <- read.csv(file="Data/TagSeq/day2.exp.data.csv", sep=',', header=TRUE) %>%  
  dplyr::mutate_if(is.character, as.factor) %>% 
  dplyr::rename('Sample.Name' = 'SapleName_readmatrix') %>% 
  dplyr::rename('pCO2' = 'OA') %>% 
  dplyr::select(c('Sample.Name','Temperature','pCO2','Salinity', 'Aragonite_saturation')) %>% 
  dplyr::mutate(Aragonite_saturation = case_when(Aragonite_saturation < 0.5 ~ 'Low', 
                                                 (Aragonite_saturation > 0.5 & Aragonite_saturation < 1.0) ~ 'Mid', 
                                                 Aragonite_saturation > 1.0 ~ 'High')) %>% 
  dplyr::mutate(All_treatment = paste( (substr(Temperature,1,1)), 
                                       (substr(pCO2,1,1)), 
                                       (substr(Salinity,1,1)), sep = '')) %>% 
  dplyr::filter(Sample.Name %in% (colnames(d2_rlog_transformed)[-1])) %>%  # call only samples that are in the rlog transformed (cut during WGCNA processing)
  dplyr::mutate(pCO2_Salinity = substr(All_treatment, 2,3)) # experiment treatment data
dim(d2.Treatment.data)

# call the module colors 
modcolor <- as.data.frame(unique(day2_MM_RRcutoff$moduleColor))
names(modcolor)[1] <- "color"


for(i in 1:nrow(modcolor)) {
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs      <- day2_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) %>%  dplyr::select("TranscriptID") %>%  na.omit()
  d2_rlog_Mod      <- d2_rlog %>% dplyr::filter(TranscriptID %in% Mod_geneIDs[,1])
  d2_rlog_Mod_MELT <- melt(d2_rlog_Mod, id=("TranscriptID")) # melt using reshape2
  names(d2_rlog_Mod_MELT)[(2:3)] <-  c('Sample.Name', 'rlog_Expression') # change column names
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d2_rlog_Mod_MELT, d2.Treatment.data, by ='Sample.Name')
  
  # mean Exp response table 
  meanEXp_Mod <- merged_Expdata_Mod %>% 
    select(c('Sample.Name','rlog_Expression','Temperature', 'Salinity', 'pCO2', 'Aragonite_saturation')) %>% 
    group_by(Sample.Name, Temperature, Salinity, pCO2, Aragonite_saturation) %>%
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
  
  # Aragonite_saturation ========================== #
  
  meanEXp_Summary.Aragonite_saturation <- meanEXp_Mod %>% 
    group_by(Aragonite_saturation) %>%
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
  
  # Aragonite sat mean sd plot ========================== #
  
  min_p4 <- min(meanEXp_Summary.Aragonite_saturation$mean) - max(meanEXp_Summary.Aragonite_saturation$se)
  max_p4 <- max(meanEXp_Summary.Aragonite_saturation$mean) + max(meanEXp_Summary.Aragonite_saturation$se)
  
  Aragonite_saturation.rlog.Mod <- meanEXp_Summary.Aragonite_saturation %>% 
    dplyr::mutate(Aragonite_saturation    = forcats::fct_relevel(Aragonite_saturation, 'Low', 'Mid', 'High')) %>%
    ggplot(aes(x=Aragonite_saturation, y=mean, fill=Aragonite_saturation)) +  # , colour=supp, group=supp))
    theme_classic() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Aragonite_saturation") +
    ylab("rlog gene expression") +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("grey", "grey", "grey")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    # ggtitle(paste("Day 7 WGCNA", modcolor[i,], "Module VST GeneExp", sep =' ')) +
    # expand_limits(y=0) +                                                    # Expand y range
    scale_y_continuous(limits=c((min_p4), (max_p4))) +
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
  # library(ggpubr)
  single.factor.plot <-  ggarrange(Temperature.rlog.Mod, Salinity.rlog.Mod,  pCO2.rlog.Mod,  Aragonite_saturation.rlog.Mod,    
                                   plotlist = NULL,
                                   ncol = 4,
                                   nrow = 1,
                                   labels = NULL)
  
  
  # Summary plot of all treatments ==================================================================================== #
  # All.Treatment mean sd plot
  min_p5 <- min(meanEXp_Summary.All.Treatment$mean) - max(meanEXp_Summary.All.Treatment$se)
  max_p5 <- max(meanEXp_Summary.All.Treatment$mean) + max(meanEXp_Summary.All.Treatment$se)
  
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
    scale_y_continuous(limits=c((min_p5), (max_p5))) +
    theme(text = element_text(size=15)) + 
    facet_wrap(~Salinity) # facetted by temperature
  
  # output   ======================================================================================================== #
  pdf(paste("Output/WGCNA/day2_larvae/ModuleExpression_Treatment_RRcutoff/day2_Exp_Module_",modcolor[i,],".pdf", sep = ''), width=9, height=8)
  print(ggarrange(single.factor.plot, AllTreatment.rlog.Mod,         
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 2,
                  labels = NULL))
  dev.off()
  
}















##################################################################################### #
####################  Day 2 RRcutoff Mod Mem plots ################################## #
##################################################################################### #

# Load data 
day18_ModuleMembership  <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership_RRcutoff.csv") %>% dplyr::select(-X)
d18_rlog                <- read.csv(file="Output/WGCNA/day18_spat/d18_rlog_transformed.csv", sep=',', header=TRUE) %>% dplyr::select(-X)


d18.Treatment.data <- read.csv(file="Data/TagSeq/day18.exp.data.csv", sep=',', header=TRUE) %>%  
  dplyr::mutate_if(is.character, as.factor) %>% 
  dplyr::rename('Sample.Name' = 'SapleName_readmatrix') %>% 
  dplyr::rename('pCO2' = 'OA') %>% 
  dplyr::select(c('Sample.Name','Temperature','pCO2','Salinity', 'Aragonite_saturation')) %>% 
  dplyr::mutate(Aragonite_saturation = case_when(Aragonite_saturation < 0.5 ~ 'Low', 
                                                 (Aragonite_saturation > 0.5 & Aragonite_saturation < 1.0) ~ 'Mid', 
                                                 Aragonite_saturation > 1.0 ~ 'High')) %>% 
  dplyr::mutate(All_treatment = paste( (substr(Temperature,1,1)), 
                                       (substr(pCO2,1,1)), 
                                       (substr(Salinity,1,1)), sep = '')) %>% 
  dplyr::filter(Sample.Name %in% (colnames(d18_rlog)[-1])) %>%  # call only samples that are in the rlog transformed (cut during WGCNA processing)
  dplyr::mutate(pCO2_Salinity = substr(All_treatment, 2,3)) # experiment treatment data
dim(d18.Treatment.data)

# call the module colors 
modcolor <- as.data.frame(unique(day18_ModuleMembership$moduleColor))
names(modcolor)[1] <- "color"


for(i in 1:nrow(modcolor)) {
  
  # vst read count date - narrow the columns - reshape and rename
  Mod_geneIDs       <- day18_ModuleMembership %>% dplyr::filter(moduleColor %in% modcolor[i,]) %>%  dplyr::select("TranscriptID") %>%  na.omit()
  d18_rlog_Mod      <- d18_rlog %>% dplyr::filter(TranscriptID %in% Mod_geneIDs[,1])
  d18_rlog_Mod_MELT <- melt(d18_rlog_Mod, id=("TranscriptID")) # melt using reshape2
  names(d18_rlog_Mod_MELT)[(2:3)] <-  c('Sample.Name', 'rlog_Expression') # change column names
  
  # merge by common row values 'Sample.Name'
  merged_Expdata_Mod <- merge(d18_rlog_Mod_MELT, d18.Treatment.data, by ='Sample.Name')
  
  # mean Exp response table 
  meanEXp_Mod <- merged_Expdata_Mod %>% 
    select(c('Sample.Name','rlog_Expression','Temperature', 'Salinity', 'pCO2', 'Aragonite_saturation')) %>% 
    group_by(Sample.Name, Temperature, Salinity, pCO2, Aragonite_saturation) %>%
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
  
  # Aragonite_saturation ========================== #
  
  meanEXp_Summary.Aragonite_saturation <- meanEXp_Mod %>% 
    group_by(Aragonite_saturation) %>%
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
    group_by(Salinity, Temperature, pCO2, Aragonite_saturation) %>%
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
  
  
  # Aragonite sat mean sd plot ========================== #
  
  Aragonite_saturation.rlog.Mod <- meanEXp_Summary.Aragonite_saturation %>% 
    dplyr::mutate(Aragonite_saturation    = forcats::fct_relevel(Aragonite_saturation, 'Low', 'Mid', 'High')) %>%
    ggplot(aes(x=Aragonite_saturation, y=mean, fill=Aragonite_saturation)) +  # , colour=supp, group=supp))
    theme_classic() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("Aragonite_saturation") +
    ylab("rlog gene expression") +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("grey", "grey", "grey")) +
    # scale_color_manual(values=c("#56B4E9","#E69F00")) +
    # ggtitle(paste("Day 7 WGCNA", modcolor[i,], "Module VST GeneExp", sep =' ')) +
    # expand_limits(y=0) +                                                    # Expand y range
    scale_y_continuous(limits=c((min_p), (max_p))) +
    theme(text = element_text(size=10), legend.position="none")
  
  
  # Assemble these together =========================================================================================== #
  # library(ggpubr)
  single.factor.plot <-  ggarrange(Temperature.rlog.Mod, Salinity.rlog.Mod,  pCO2.rlog.Mod,  Aragonite_saturation.rlog.Mod,    
                                   plotlist = NULL,
                                   ncol = 4,
                                   nrow = 1,
                                   labels = NULL)
  
  
  # Summary plot of all treatments ==================================================================================== #
  # All.Treatment mean sd plot
  min_p <- min(meanEXp_Summary.All.Treatment$mean) - max(meanEXp_Summary.All.Treatment$se)
  max_p <- max(meanEXp_Summary.All.Treatment$mean) + max(meanEXp_Summary.All.Treatment$se)

  AllTreatment.rlog.Mod <- meanEXp_Summary.All.Treatment %>% 
    dplyr::mutate(Salinity    = forcats::fct_relevel(Salinity, 'Low', 'High')) %>%
    dplyr::mutate(pCO2        = forcats::fct_relevel(pCO2, 'Low', 'High')) %>%
    dplyr::mutate(Aragonite_saturation = forcats::fct_relevel(Aragonite_saturation, 'High', 'Mid', 'Low')) %>%
    ggplot(aes(x=pCO2, y=mean, fill=Aragonite_saturation, group=Aragonite_saturation)) + # group aesthetic connect line (Slaintiy) and color - the x axis in this case is pCO2 
    theme_classic() +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), colour="black", width=.1, position=pd) +
    geom_line(position=pd) +
    geom_point(position=pd, size = 4, shape=21) +            
    xlab("pCO2") +
    ylab("rlog gene expression") +                 # note the mean was first by sample ID THEN by treatment
    # ylab(paste(modcolor[i,]," Module rlog Gene Expression (Mean +/- SE)", sep = ' ')) +                 # note the mean was first by sample ID THEN by treatment
    scale_fill_manual(values=c("#D55E00", "#E69F00", "#56B4E9")) +
    scale_y_continuous(limits=c((min_p), (max_p))) +
    theme(text = element_text(size=15)) + 
    facet_wrap(~Salinity) # facetted by temperature
  
  # output   ======================================================================================================== #
  pdf(paste("Output/WGCNA/day18_spat/ModuleExpression_Treatment_RRcutoff/day18_Exp_Module_",modcolor[i,],".pdf", sep = ''), width=9, height=8)
  print(ggarrange(single.factor.plot, AllTreatment.rlog.Mod,         
                  plotlist = NULL,
                  ncol = 1,
                  nrow = 2,
                  labels = NULL))
  dev.off()
  
}
