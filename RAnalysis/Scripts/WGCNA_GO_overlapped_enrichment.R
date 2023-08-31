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


# RR cutoff GO Enrichment results 

day2_GO_RR      <- read.csv("Output/WGCNA/day2_larvae/GO_analysis/RR_cutoff/Day2_GO_Results_Master.csv") 

day22_GO_RR      <- read.csv("Output/WGCNA/day18_spat/GO_analysis/RR_cutoff/Day18_GO_Results_Master.csv") 

WGCNA_MasterGO   <-  merge( (as.data.frame(rbind(day2_GO_RR, 
                                                 day22_GO_RR)) %>% 
                               dplyr::rename(Cvirginica_TranscriptID=Transcript_ID)), 
                            Master_ref, by="Cvirginica_TranscriptID")

# ============================================================================= #
# Co-expression Pattern HIGH from low salnity and LOW from elevated pCO2 
# ============================================================================= #
# About: 
# here we will look at genes in D22 module red that are not present in D22 tan
# D22 red = Mod_arag_low_sal > mod_arag_high_pCO2
# D22 tan = mog_arag has low expression (regardless of low salinity or elevated pCO2 in isolattion)
# these represent unique genes we know are only elevated by moderate aragonite caused by low salinity (reduced by pCO2)
# therefore genes uniquely present in module red and not in tan filters those without overlapped function 

d22tan <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'tan'))

d22red <- WGCNA_MasterGO %>% 
  dplyr::filter(Day == 'Day18' & moduleColor %in% 'red')

# create a vector to filter with 
d22tan_filter          <- d22tan %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22tan_proteinNames    <- paste(substr(d22tan_filter$Protein_name,1,15), collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

d22red_filter          <- d22red %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22red_proteinNames    <- paste(substr(d22red_filter$Protein_name,1,15), collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

# subsets
d22red_SUBSET  <- subset(d22red_filter, !grepl(d22tan_proteinNames, Protein_name)) 
d22red_SUBSET$Protein_name
nrow(d22red_filter) - nrow(d22red_SUBSET) # only 10 genes omitted in this subsetting


# ============================================================================= #
# Co-expression Pattern HIGH  Expression under  LOW Aragonite (high salinity + low pCO2)
# ============================================================================= #
# About: 
# here we will look at overlapped genes in D2 Turquoise, D22 salmon + D22 turquoise
d2turquoise <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day2' ~ moduleColor %in% 'turquoise'))

d22turquoise  <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'turquoise')) 

d22salmon  <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'salmon')) 

# create a vector to filter with 
d2turquoise_filter         <- d2turquoise %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d2turquoise_proteinNames   <- paste(d2turquoise_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

d22turquoise_filter        <- d22turquoise %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22turqupoise_proteinNames <- paste(d22turquoise_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

d22salmon_filter           <- d22salmon %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22salmon_proteinNames     <- paste(d22salmon_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 


# genes enriched SHARED between both d22 and t2 turquoise
d22turquoise_SUBSET  <- subset(d22turquoise, grepl(d2turquoise_proteinNames, Protein_name)) 
sort(unique(d22turquoise_SUBSET$Protein_name))
# [1] " calmodulin-like"                                           " caveolin-1-like"                                          
# [3] " enolase-phosphatase E1-like"                               " eukaryotic translation initiation factor 3 subunit E-like"
# [5] " eukaryotic translation initiation factor 3 subunit G-like" " eukaryotic translation initiation factor 3 subunit L-like"
# [7] " group XIIA secretory phospholipase A2-like"                " grpE protein homolog 1, mitochondrial-like"               
# [9] " histone H4"                                                " nucleoside diphosphate kinase 7-like"                     
# [11] " prefoldin subunit 1-like"                                  " proteasome subunit alpha type-6-like"                     
# [13] " ruvB-like 1"                                               " signal peptidase complex subunit 1-like"                  
# [15] " spliceosome RNA helicase DDX39B"                           " T-complex protein 1 subunit alpha-like"                   
# [17] " T-complex protein 1 subunit beta-like"                     " T-complex protein 1 subunit delta-like"                   
# [19] " T-complex protein 1 subunit zeta-like"                     " transcription initiation factor TFIID subunit 13-like"

# genes unique to d2 turquioise - ommit those that are shared with d22 turquoise 
d2turquoise_SUBSET  <- d2turquoise %>% filter(!(Protein_name %in% unique(d22turquoise_SUBSET$Protein_name))) 
sort(d2turquoise_SUBSET$Protein_name)
# zinc finger,, carnosine synthase, EIF3, glutathione peroxidase 2, 
# 

# ============================================================================= #
# Co-expression Pattern: LOW Expression under Low Aragonite
# ============================================================================= #
# About: 
# here we are looking for the following 
# overlap between D1 blue and D22 green - each have significnat low expression from aragonite low)

d22green <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'green'))

d22blue  <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day18' ~ moduleColor %in% 'blue')) 

d2blue  <- WGCNA_MasterGO %>% 
  dplyr::filter(case_when(Day == 'Day2' ~ moduleColor %in% 'blue')) 

# create a vector to filter with 
d22green_filter            <- d22green %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22green_proteinNames      <- paste(d22green_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

d22blue_filter             <- d22blue %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22blue_proteinNames       <- paste(d22blue_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

d2blue_filter             <- d2blue %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d2blue_proteinNames       <- paste(d2blue_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

# subsets
d22green_SUBSET  <- subset(d22green_filter, grepl(d2blue_proteinNames, Protein_name)) 
d22green_SUBSET$Protein_name # all of green D22 that also contains blue D1

d22blue_SUBSET  <- subset(d22blue_filter, grepl(d2blue_proteinNames, Protein_name)) 
sort(d22blue_SUBSET$Protein_name) # all of green D22 that also contains blue D1


# genes unique to d2 blue - ommit those that are shared with d22 blue 
d2blue_SUBSET  <- d2blue %>% filter(!(Protein_name %in% unique(d22blue_SUBSET$Protein_name))) 
sort(d2blue_SUBSET$Protein_name)
# 


# (2) Take these UNIQUE module D22 Green genes - which overlap with module D22 Tan?

d22tan <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'tan'))

# create a vector to filter with 
d22tan_filter              <- d22tan %>% dplyr::filter(!Protein_name %in% ' uncharacterized') # removed "uncharacterized from protein names
d22tan_proteinNames        <- paste(d22tan_filter$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

d22geenUNIQUE_proteinNames <- paste(d22green_SUBSET$Protein_name, collapse="|") # separate all protein names by | to use grep to subset anoter dataframe 

# subsets
d22tan_SUBSET  <- subset(d22tan_filter, grepl(d22geenUNIQUE_proteinNames, Protein_name)) 
d22tan_SUBSET$Protein_name

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

d2turquoise <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day2' ~ moduleColor %in% 'turquoise')) 

d22turquoise <- WGCNA_MasterModData %>% # low salinity and mid aragonite == high expression
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'turquoise'))

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

d22blue <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'blue')) %>% 
  dplyr::select(c('Protein_name', 'moduleColor')) %>% 
  unique()


d22green  <- WGCNA_MasterModData %>% 
  dplyr::filter(case_when(Day == 'Day22' ~ moduleColor %in% 'green')) %>% 
  dplyr::select(c('Protein_name', 'moduleColor')) %>% 
  unique()

#merge
modules_LowAragoniteLowExp_Day22 <- rbind(d22blue, d22green)
View(modules_LowAragoniteLowExp_Day22)
#do genes overlap?
View(modules_LowAragoniteLowExp_Day22 %>% 
       group_by(Protein_name) %>% 
       #group_by(Protein_name) %>% 
       dplyr::summarise(n = n()) %>%  
       filter(!n %in% 1))

