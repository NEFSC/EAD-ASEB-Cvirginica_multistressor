library(dplyr)
library(forcats)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(tidyr)
library(tidyverse)

# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")
path_out = 'C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/WGCNA/' # personnal computer




# master reference view te master ref R script for details
Master_ref  <- read.csv(file= "Data/TagSeq/Seq_details/Seq_Reference_Master.csv", sep=',', header=TRUE) %>% 
                  dplyr::select(c('Cvirginica_TranscriptID','Annotation_GO_ID', 'Cvirginica_length'))


Cvirginica_annot_reference  <- read.csv(file="Data/TagSeq/Seq_details/seq_id_master.csv", sep=',', header=TRUE) %>% 
                                dplyr::select(c('TranscriptID','Function','GeneID')) %>% 
                                dplyr::mutate(TranscriptID = gsub(" ", "", TranscriptID)) %>% # remove the space at the end of each transcript ID
                                dplyr::mutate(Protein_name = gsub("\\s\\(LOC.*|\\sLOC111.*", "", perl=TRUE, Function)) %>% 
                                dplyr::select(!Function)


# WGCNA results -fromat and merge with the GO terms

# Day 2 - all including low nad high temperature :::::::::::::::::::::::::::::::::::::::::::::::; #

d2_Annot_ModuleMembership      <- read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv") %>% 
                                    dplyr::select(c('TranscriptID','geneSymbol','Protein_name','moduleColor')) %>%  
                                    na.omit() %>% 
                                    dplyr::mutate(Day = "Day2")

d2ModCols                      <- data.frame(moduleColor = unique(d2_Annot_ModuleMembership$moduleColor)) %>% # unique module colors 
                                    dplyr::filter(moduleColor %in% c('pink', 
                                                                     'blue', 
                                                                     'turquoise', 
                                                                     'brown', 
                                                                     'black', 
                                                                     'red')) %>% # sig modules
                                    dplyr::mutate(Day = "Day2")



# Day 2 - high temperature ONLY ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::; #

d2_Annot_ModuleMembership_hightemp  <- read.csv("Output/WGCNA/day2_larvae_hightemp/d2.WGCNA_ModulMembership.csv") %>% 
                                          dplyr::select(c('TranscriptID','geneSymbol','Protein_name','moduleColor')) %>%  
                                          na.omit() %>% 
                                          dplyr::mutate(Day = "Day2_hightemp")

d2ModCols_hightemp                  <- data.frame(moduleColor = unique(d2_Annot_ModuleMembership_hightemp$moduleColor)) %>% # unique modules
                                          dplyr::filter(moduleColor %in% c('red', 
                                                                           'greenyellow', 
                                                                           'brown', 
                                                                           'pink', 
                                                                           'magenta')) %>% # sig modules
                                          dplyr::mutate(Day = "Day2_hightemp")


# Day 18  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::; #

d18_Annot_ModuleMembership     <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv") %>% 
                                    dplyr::select(c('TranscriptID','geneSymbol','Protein_name','moduleColor')) %>%  
                                    na.omit() %>% 
                                    dplyr::mutate(Day = "Day18")

d18ModCols                     <- data.frame(moduleColor = unique(d18_Annot_ModuleMembership$moduleColor)) %>% # unique module colors 
                                    dplyr:: filter(moduleColor %in% c('tan', 
                                                                      'red', 
                                                                      'turquoise', 
                                                                      'salmon', 
                                                                      'blue',
                                                                      'green',
                                                                      'greenyellow')) %>% # sig modules
                                    dplyr::mutate(Day = "Day18")




# master WGCNA module data for for loops! 
WGCNA_MasterModData   <-  merge( (as.data.frame(rbind(d2_Annot_ModuleMembership, 
                                                      d2_Annot_ModuleMembership_hightemp, 
                                                      d18_Annot_ModuleMembership)) %>% 
                                    dplyr::rename(Cvirginica_TranscriptID=TranscriptID)), 
                                 Master_ref, by="Cvirginica_TranscriptID")

WGCNA_ColorList       <-  rbind(d2ModCols, d2ModCols_hightemp, d18ModCols) # master WGCNA color list - use this to loop all the analysis 



# ============================================================================= #
# Co-expression Pattern: High Salinity == High Expression
# ============================================================================= #

d2_brown <- WGCNA_MasterModData %>% 
                dplyr::filter(case_when(Day == 'Day2_hightemp' ~ moduleColor %in% 'brown'))

d18blue  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'blue')) 

#merge
modules_HighSalHighExp <- rbind(d2_brown, d18blue)

#do genes overlap?
View(modules_HighSalHighExp %>% 
  group_by(Cvirginica_TranscriptID,Protein_name) %>% 
  dplyr::summarise(n = n()) %>%  
  filter(!n %in% 1))


# ============================================================================= #
# Co-expression Pattern: Low Salinity == High Expression
# ============================================================================= #

d2red <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day2_hightemp' ~ moduleColor %in% 'red'))

d18red.salmon  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% c('red', 'salmon'))) 

#merge
modules_LowSalHighExp <- rbind(d2red, d18red.salmon)

#do genes overlap?
View(modules_LowSalHighExp %>% 
       group_by(Cvirginica_TranscriptID,Protein_name) %>% 
       dplyr::summarise(n = n()) %>%  
       filter(!n %in% 1))


# ============================================================================= #
# Co-expression Pattern: Low Aragonite == Low Expression
# ============================================================================= #

d2magenta <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day2_hightemp' ~ moduleColor %in% 'magenta'))

d18green  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'green')) 

#merge
modules_LowAragoniteLowExp <- rbind(d2magenta, d18green)

#do genes overlap?
View(modules_LowAragoniteLowExp %>% 
       group_by(Cvirginica_TranscriptID,Protein_name) %>% 
       #group_by(Protein_name) %>% 
       dplyr::summarise(n = n()) %>%  
       filter(!n %in% 1))


# ============================================================================= #
# Co-expression Pattern: Low Aragonite == High Expression
# ============================================================================= #

d2red <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day2_hightemp' ~ moduleColor %in% 'red'))

d18turquoise  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'turquoise')) 

#merge
modules_LowAragoniteHighExp <- rbind(d2red, d18turquoise)

#do genes overlap?
View(modules_LowAragoniteHighExp %>% 
       group_by(Cvirginica_TranscriptID,Protein_name) %>% 
       #group_by(Protein_name) %>% 
       dplyr::summarise(n = n()) %>%  
       filter(!n %in% 1))

