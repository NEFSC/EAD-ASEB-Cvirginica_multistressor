---
  # title: "WGCNA overlapped genes"
  # author: "Samuel Gurr"
  # date: "Oct 5, 2022"
---

#  LOAD LIBRARIES
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

# Day 2 - all including low and high temperature :::::::::::::::::::::::::::::::::::::::::::::::; #

d2_Annot_ModuleMembership      <- read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv") 
# for loop to call the threshold or module mdembership
d2ModCols                      <- data.frame(moduleColor = unique(d2_Annot_ModuleMembership$moduleColor)) %>% # unique module colors 
                                    dplyr::filter(moduleColor %in% c('pink', 
                                                                     'blue', 
                                                                     'turquoise', 
                                                                     'brown', 
                                                                     'black', 
                                                                     'red')) %>% # sig modules
                                    dplyr::mutate(Day = "Day2")
df_total              <- data.frame() # start dataframe 
day2_mod_RR           <- data.frame(matrix(nrow = 1, ncol = 5)) # create dataframe to save cumunalitively during for loop
colnames(day2_mod_RR) <- c('Day', 'moduleColor', 'TranscriptID','geneSymbol','Protein_name') # names for comuns in the for loop
for (i in 1:nrow(d2ModCols)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor         <- d2ModCols[i,1]
  loopmodColor_cor <- paste("MM.", modColor, sep = '') # column name for mod color - PEarsons correlation value 
  loopmodColor_p   <- paste("p.MM.", modColor, sep = '') # column name for mod color - Students asymptotic p value 
  
  Mod_loop_d2            <- d2_Annot_ModuleMembership %>% 
    dplyr::filter(moduleColor %in% modColor) %>% 
    dplyr::select(c('TranscriptID','geneSymbol','Protein_name','moduleColor', loopmodColor_cor,loopmodColor_p))
  Mod_Loop_d2_RRcutoff   <- as.data.frame(Mod_loop_d2 %>% 
                                            dplyr::filter(!(.[[5]] <  0.6 & .[[6]] > 0.05)) %>% 
                           dplyr::mutate(Day = "Day2")) %>% 
                           dplyr::select(!c(5,6))
  
  
  # write csv file for the data reduced mod mem 
  
  loopdf       <- data.frame(Mod_Loop_d2_RRcutoff) # name dataframe for this single row
  day2_mod_RR  <- rbind(day2_mod_RR,loopdf) #bind to a cumulative list dataframe
  print(day2_mod_RR) # print to monitor progress
  
}
View(day2_mod_RR)

# Day 22  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::; #

d22_Annot_ModuleMembership     <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv") 
d22ModCols                     <- data.frame(moduleColor = unique(d22_Annot_ModuleMembership$moduleColor)) %>% # unique module colors 
                                    dplyr:: filter(moduleColor %in% c('tan', 
                                                                      'red', 
                                                                      'turquoise', 
                                                                      'salmon', 
                                                                      'blue',
                                                                      'green',
                                                                      'greenyellow')) %>% # sig modules
                                    dplyr::mutate(Day = "Day22")
df_total              <- data.frame() # start dataframe 
day22_mod_RR          <- data.frame(matrix(nrow = 1, ncol = 5)) # create dataframe to save cumunalitively during for loop
colnames(day22_mod_RR) <- c('Day', 'moduleColor', 'TranscriptID','geneSymbol','Protein_name') # names for comuns in the for loop
for (i in 1:nrow(d22ModCols)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor         <- d22ModCols[i,1]
  loopmodColor_cor <- paste("MM.", modColor, sep = '') # column name for mod color - PEarsons correlation value 
  loopmodColor_p   <- paste("p.MM.", modColor, sep = '') # column name for mod color - Students asymptotic p value 
  
  Mod_loop_d22            <- d22_Annot_ModuleMembership %>% 
    dplyr::filter(moduleColor %in% modColor) %>% 
    dplyr::select(c('TranscriptID','geneSymbol','Protein_name','moduleColor', loopmodColor_cor,loopmodColor_p))
  Mod_Loop_d22_RRcutoff   <- as.data.frame(Mod_loop_d22 %>% 
                                            dplyr::filter(!(.[[5]] <  0.6 & .[[6]] > 0.05)) %>% 
                                            dplyr::mutate(Day = "Day22")) %>% 
    dplyr::select(!c(5,6))
  
  loopdf       <- data.frame(Mod_Loop_d22_RRcutoff) # name dataframe for this single row
  day22_mod_RR  <- rbind(day22_mod_RR,loopdf) #bind to a cumulative list dataframe
  print(day22_mod_RR) # print to monitor progress
  
}

View(day22_mod_RR)

# master WGCNA module data for for loops! 
WGCNA_MasterModData   <-  merge( (as.data.frame(rbind(day2_mod_RR, 
                                                     # d2_Annot_ModuleMembership_hightemp, 
                                                     day22_mod_RR)) %>% 
                                    dplyr::rename(Cvirginica_TranscriptID=TranscriptID)), 
                                 Master_ref, by="Cvirginica_TranscriptID")

#WGCNA_ColorList       <-  rbind(d2ModCols, d18ModCols) # master WGCNA color list - use this to loop all the analysis 



# ============================================================================= #
# Co-expression Pattern: High Salinity == High Expression
# ============================================================================= #

d2blue <- WGCNA_MasterModData %>% 
                dplyr::filter(case_when(Day == 'Day2' ~ moduleColor %in% 'blue'))

d22blue  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'blue')) 

#merge
modules_HighSalHighExp <- rbind(d2blue, d22blue)

#do genes overlap?
modules_HighSalHighExp_filtered <- as.data.frame(modules_HighSalHighExp %>% 
                                                  group_by(Cvirginica_TranscriptID,Protein_name) %>% 
                                                  dplyr::summarise(n = n()) %>%  
                                                  filter(!n %in% 1) %>% 
                                                  select(!n))
write.csv(modules_HighSalHighExp_filtered, 
          "C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/WGCNA/Shared_Day2blue.Day22blue_HighSal_HighExp.csv")# write


# ============================================================================= #
# Co-expression Pattern: Low Salinity == High Expression - focus on high expression under low aragonite
# ============================================================================= #

d22red <- WGCNA_MasterModData %>% # low salinity and mid aragonite == high expression
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'turquoise'))

d2turquoise <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day2' ~ moduleColor %in% 'turquoise')) 

#merge
modules_LowSalHighExp <- rbind(d22red, d2turquoise)

#do genes overlap?
modules_LowSalHighExp_filtered <- as.data.frame(modules_LowSalHighExp %>% 
       group_by(Cvirginica_TranscriptID,Protein_name) %>% 
       dplyr::summarise(n = n()) %>%  
       filter(!n %in% 1) %>% 
       select(!n))

write.csv(modules_LowSalHighExp_filtered, 
          "C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/WGCNA/Shared_Day2turquoise.Day22turquoise_LowArag_HighExp.csv")# write


# ============================================================================= #
# Co-expression Pattern: Low Aragonite == Low Expression
# ============================================================================= #


d2blue <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day2' ~ moduleColor %in% 'blue'))

d22green  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'green')) 

#merge
modules_LowAragoniteLowExp <- rbind(d2blue, d22green)

#do genes overlap?
modules_LowAragoniteLowExp_filtered <- as.data.frame(modules_LowAragoniteLowExp %>% 
                                                   group_by(Cvirginica_TranscriptID,Protein_name) %>% 
                                                   dplyr::summarise(n = n()) %>%  
                                                   filter(!n %in% 1) %>% 
                                                   select(!n))
write.csv(modules_LowAragoniteLowExp_filtered, 
          "C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/WGCNA/Shared_Day2blue.Day22green_LowArag_LowExp.csv")# write



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

