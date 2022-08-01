---
# title: "Overlapped genes WGCNa DESeq2"
# author: "Samuel Gurr"
# date: "January 8, 2021"
---

  
# OBJECTIVE: 
# When interpretting the WGCNA resutls, we have a multitude of genes associated with enriched GO terms
# If we focus on main-effect WGCNA modules (i.e. higher expression ue to single treamtent factor like salinity), 
# we can narrow the focus of these by investigating the following questions: 
#   (1) What proportion of main-effect DEGs (via DESeq2) are shown in corresponding modules representing the same expression-treatment association?
#   i.e. DEGs showing upregualtion under high v. low salinity and a WGCAN module showing higher expressio under high salinity 
#   (2) Which genes are present as BOTH DEGs and WGCNA (following 1) and within enriched terms? 

# the genes here can give us an idea of gens with strong differential expressio due to MAIN TREATMENT EFFECTS! 



# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")
# load in the count matrix with rownames as sample ID and colnames as gene ID

AllDEGs         <- read.csv(file="Output/DESeq2/DE_All_Annotation.csv", sep=',', header=TRUE) 

# WGNA 
Day2All_WGNCA   <- read.csv(file="Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv", sep=',', header=TRUE) %>%  dplyr::rename(Cvirginica_TranscriptID = X)
Day2GO_WGCNA    <- read.csv(file="Output/WGCNA/day2_larvae/GO_analysis/Day2_GO_Results_Master.csv", sep=',', header=TRUE) %>%  dplyr::rename(Cvirginica_TranscriptID = Transcript_ID)

Day18All_WGNCA  <- read.csv(file="Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv", sep=',', header=TRUE) %>%  dplyr::rename(Cvirginica_TranscriptID = X)
Day18GO_WGCNA   <- read.csv(file="Output/WGCNA/day18_spat/GO_analysis/Day18_GO_Results_Master.csv", sep=',', header=TRUE) %>%  dplyr::rename(Cvirginica_TranscriptID = Transcript_ID)


# ------------------------------------------------------------------------------------------------------------------ #
# filter the WGCNA and  DESeq2 datasets by the experimen-treatment patterns (putative overlap datasets!) 
# ------------------------------------------------------------------------------------------------------------------ #

d2DEGs       <- AllDEGs %>% 
                  dplyr::filter(time %in% 'd2_larva')

d18DEGs      <- AllDEGs %>% 
                  dplyr::filter(time %in% 'd18_spat')



# what effects shall we follow through wit??

unique(d2DEGs$Effect)  # "Day2_Salinity"    "Day2_Temperature"  "Day2_OA" - each with upregulated (High expression under High __ ) and downregulated (High expression under Low __ )
unique(d18DEGs$Effect) # "Day18_OA"       "Day18_Salinity" - each with upregulated (High expression under High __ ) and downregulated (High expression under Low __ )




# Day 2: HIGH SALINITY = HIGH EXPRESSION  ============================================================== #


d2_SalHigh_DEGs     <- d2DEGs %>% 
                         dplyr::filter(Dir %in% 'up' & Effect %in% 'Day2_Salinity') # DEGS 

d2_SalHigh_WGCNA    <- Day2All_WGNCA %>% 
                         dplyr::filter(moduleColor %in% c('brown', 'blue')) # ALL MODULE MEMBERSHUP 

d2_SalHigh_WGCNA_GO <- Day2GO_WGCNA %>% 
                         dplyr::filter(moduleColor %in% c('brown','blue')) # GENES IN ENRICHED GO TERMS (WGCNA MODULE(S))

# Overlap - DEGs WITH Whole module membership 
paste( "Numer of genes shared from DESeq2 =",
       length(intersect(d2_SalHigh_DEGs$Cvirginica_TranscriptID,  d2_SalHigh_WGCNA$Cvirginica_TranscriptID)),
       "; proprotion of total DEGs is ",
       (( length(intersect(d2_SalHigh_DEGs$Cvirginica_TranscriptID,  d2_SalHigh_WGCNA$Cvirginica_TranscriptID))) / nrow(d2_SalHigh_DEGs) * 100)) 
# 90.25  gemes are shared in the enriched GO terms with the same WGCNA expression-treatment pattern 
# 69% of DEGs (220 of 318) are IN MODULE BLUE - 21% (220 of 318) IN MODULE BROWN


# Overlap - DEGs WITH  GO Enrichment 
paste( "Numer of genes shared from DESeq2 =",
    length(intersect(d2_SalHigh_DEGs$Cvirginica_TranscriptID,  d2_SalHigh_WGCNA_GO$Cvirginica_TranscriptID)),
    "; proprotion of total DEGs is ",
   (( length(intersect(d2_SalHigh_DEGs$Cvirginica_TranscriptID,  d2_SalHigh_WGCNA_GO$Cvirginica_TranscriptID))) / nrow(d2_SalHigh_DEGs) * 100)) # 124
# 38.99  gemes are shared in the enriched GO terms with the same WGCNA expression-treatment pattern 


# call the unique IDs to filter out the overlap protein names...
IDs_d2_SalHigh_GO <- as.data.frame(intersect(d2_SalHigh_DEGs$Cvirginica_TranscriptID,  d2_SalHigh_WGCNA_GO$Cvirginica_TranscriptID)) # shared with DEGs (40% of them here)
colnames(IDs_d2_SalHigh_GO) = "Cvirginica_TranscriptID" # change the col name to filter easily below 
d2_HighSal_GO..verlap <-  d2_SalHigh_WGCNA_GO %>%  # Get the genes for the GO overlap!!!!!
                              dplyr::select(c('Cvirginica_TranscriptID', 'GO_term', 'Protein_name')) %>% 
                              dplyr::filter(Cvirginica_TranscriptID %in% IDs_d2_SalHigh_GO$Cvirginica_TranscriptID)
unique(d2_HighSal_GO..verlap$Protein_name) # 90 genes - view the protein names! 


# data for the ven diagram 
nrow(d2_SalHigh_DEGs) # 318 # of DEGs
nrow(d2_SalHigh_WGCNA) # 1460 # of genes in module(s) with putatively same expression-treatment association (in this case higher under High Salinity) 
nrow(d2_SalHigh_WGCNA_filt) # total 287 overlapped





# module membership alue 
d2_SalHigh_WGCNA_blueMM_DEGfilered <- d2_SalHigh_WGCNA %>% 
                                        dplyr::filter(moduleColor %in% 'blue') %>%  
                                        dplyr::filter(Cvirginica_TranscriptID %in% d2_SalHigh_DEGs$Cvirginica_TranscriptID) %>% 
                                        dplyr::select(MM.blue)
d2_SalHigh_WGCNA_blueMM <- d2_SalHigh_WGCNA %>% 
                              dplyr::filter(moduleColor %in% 'blue') %>%  
                              dplyr::select(MM.blue)

# module membership alue 
d2_SalHigh_WGCNA_bluepMM_DEGfilered <- d2_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') %>%  
  dplyr::filter(Cvirginica_TranscriptID %in% d2_SalHigh_DEGs$Cvirginica_TranscriptID) %>% 
  dplyr::select(p.MM.blue)
d2_SalHigh_WGCNA_bluepMM <- d2_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') %>%  
  dplyr::select(p.MM.blue)

par(mfrow = c(2,2))
MASS::truehist(d2_SalHigh_WGCNA_bluepMM_DEGfilered$p.MM.blue, col = 'grey85', nbins = 20, ymax = 500, xlim = c(0,1), main = "p.MM.filtered_mod_blue", prob = FALSE)
MASS::truehist(d2_SalHigh_WGCNA_bluepMM$p.MM.blue, col = 'white', nbins = 25, ymax = 500, xlim = c(0,1), main = "p.MM.all_mod_blue", prob = FALSE)
MASS::truehist(d2_SalHigh_WGCNA_blueMM_DEGfilered$MM.blue, col = 'grey85', nbins = 25, ymax = 100, xlim = c(0,1), main = "MM.filtered_mod_blue", prob = FALSE)
MASS::truehist(d2_SalHigh_WGCNA_blueMM$MM.blue, col = 'white', nbins = 25, ymax = 100, xlim = c(0,1), main = "MM.all_mod_blue", prob = FALSE)


(length(d2_SalHigh_WGCNA_bluepMM_DEGfilered$p.MM.blue < 0.5) + length(d2_SalHigh_WGCNA_brownpMM_DEGfilered$p.MM.brown < 0.5)) /(nrow(d2_SalHigh_DEGs)) * 100 # 90.25157 % of the total DEGs are in WGCNA moduleat p va;ue < 0.05

# ALL DEGs that fell within the WGCNA modeuls of the same effect were present here (~90% or 287 of the DEGs) 



d2_SalHigh_WGCNA_brownpMM_DEGfilered <- d2_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'brown') %>%  
  dplyr::filter(Cvirginica_TranscriptID %in% d2_SalHigh_DEGs$Cvirginica_TranscriptID) %>% 
  dplyr::select(p.MM.brown)


# What is the mean ModuleMembership calue? (Pearsons correlation of eigengen velues)
cor_d2_HighSal <- d2_SalHigh_WGCNA  %>%  dplyr::filter(Cvirginica_TranscriptID %in% d2_SalHigh_DEGs$Cvirginica_TranscriptID & 
                                                          case_when(moduleColor=="blue" ~ p.MM.blue < 0.5))


mean(cor_d2_HighSal$MM.blue) # 0.7325375
sd(cor_d2_HighSal$MM.blue)   # 0.1804349


cor_d2_HighSal <- d2_SalHigh_WGCNA  %>%  dplyr::filter(Cvirginica_TranscriptID %in% d2_SalHigh_DEGs$Cvirginica_TranscriptID & 
                                                          case_when(moduleColor=="brown" ~ p.MM.brown < 0.5))


mean(cor_d2_HighSal$MM.brown) # 0.6400661
sd(cor_d2_HighSal$MM.brown)   # 0.1637451







# Day 18: HIGH SALINITY = HIGH EXPRESSION  ============================================================== #


d18_SalHigh_DEGs     <- d18DEGs %>% 
  dplyr::filter(Dir %in% 'up' & Effect %in% 'Day18_Salinity') # DEGS 

d18_SalHigh_WGCNA    <- Day18All_WGNCA %>% 
  dplyr::filter(moduleColor %in% 'blue') # ALL MODULE MEMBERSHUP 

d18_SalHigh_WGCNA_GO <- Day18GO_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') # GENES IN ENRICHED GO TERMS (WGCNA MODULE(S))

# Overlap - DEGs WITH Whole module membership 
paste( "Numer of genes shared from DESeq18 =",
       length(intersect(d18_SalHigh_DEGs$Cvirginica_TranscriptID,  d18_SalHigh_WGCNA$Cvirginica_TranscriptID)),
       "; proprotion of total DEGs is ",
       (( length(intersect(d18_SalHigh_DEGs$Cvirginica_TranscriptID,  d18_SalHigh_WGCNA$Cvirginica_TranscriptID))) / nrow(d18_SalHigh_DEGs) * 100)) 
# 79  gemes are shared in the enriched GO terms with the same WGCNA expression-treatment pattern 


# Overlap - DEGs WITH  GO Enrichment 
paste( "Numer of genes shared from DESeq18 =",
       length(intersect(d18_SalHigh_DEGs$Cvirginica_TranscriptID,  d18_SalHigh_WGCNA_GO$Cvirginica_TranscriptID)),
       "; proprotion of total DEGs is ",
       (( length(intersect(d18_SalHigh_DEGs$Cvirginica_TranscriptID,  d18_SalHigh_WGCNA_GO$Cvirginica_TranscriptID))) / nrow(d18_SalHigh_DEGs) * 100)) # 1184
# 17 gemes are shared in the enriched GO terms with the same WGCNA expression-treatment pattern 

# call the unique IDs to filter out the overlap protein names...
IDs_d18_SalHigh_GO <- as.data.frame(intersect(d18_SalHigh_DEGs$Cvirginica_TranscriptID,  d18_SalHigh_WGCNA_GO$Cvirginica_TranscriptID)) # shared with DEGs (40% of them here)
colnames(IDs_d18_SalHigh_GO) = "Cvirginica_TranscriptID" # change the col name to filter easily below 
d18_HighSal_GO..verlap <-  d18_SalHigh_WGCNA_GO %>%  # Get the genes for the GO overlap!!!!!
  dplyr::select(c('Cvirginica_TranscriptID', 'GO_term', 'Protein_name')) %>% 
  dplyr::filter(Cvirginica_TranscriptID %in% IDs_d18_SalHigh_GO$Cvirginica_TranscriptID)
unique(d18_HighSal_GO..verlap$Protein_name) # 90 genes - view the protein names! 


# data for the ven diagram 
nrow(d18_SalHigh_DEGs) # 476 # of DEGs
nrow(d18_SalHigh_WGCNA) # 808 # of genes in module(s) with putatively same expression-treatment association (in this case higher under High Salinity) 
nrow(d18_SalHigh_WGCNA_blueMM_DEGfilered) # total 378 overlapped



# module membership alue 
d18_SalHigh_WGCNA_blueMM_DEGfilered <- d18_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') %>%  
  dplyr::filter(Cvirginica_TranscriptID %in% d18_SalHigh_DEGs$Cvirginica_TranscriptID) %>% 
  dplyr::select(MM.blue)
d18_SalHigh_WGCNA_blueMM <- d18_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') %>%  
  dplyr::select(MM.blue)

# module membership alue 
d18_SalHigh_WGCNA_bluepMM_DEGfilered <- d18_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') %>%  
  dplyr::filter(Cvirginica_TranscriptID %in% d18_SalHigh_DEGs$Cvirginica_TranscriptID) %>% 
  dplyr::select(p.MM.blue)
d18_SalHigh_WGCNA_bluepMM <- d18_SalHigh_WGCNA %>% 
  dplyr::filter(moduleColor %in% 'blue') %>%  
  dplyr::select(p.MM.blue)

par(mfrow = c(2,2))
MASS::truehist(d18_SalHigh_WGCNA_bluepMM_DEGfilered$p.MM.blue, nbins = 10, ymax = 600, xlim = c(0,1), main = "p.MM.filtered_mod_blue", prob = FALSE)
MASS::truehist(d18_SalHigh_WGCNA_bluepMM$p.MM.blue, nbins = 10, ymax = 600, xlim = c(0,1), main = "p.MM.all_mod_blue", prob = FALSE)
MASS::truehist(d18_SalHigh_WGCNA_blueMM_DEGfilered$MM.blue, nbins = 25, ymax = 150, xlim = c(0,1), main = "MM.filtered_mod_blue", prob = FALSE)
MASS::truehist(d18_SalHigh_WGCNA_blueMM$MM.blue, nbins = 25, ymax = 150, xlim = c(0,1), main = "MM.all_mod_blue", prob = FALSE)


(length(d18_SalHigh_WGCNA_bluepMM_DEGfilered$p.MM.blue < 0.5)) /(nrow(d18_SalHigh_DEGs)) * 100 # 79.41176 % of the total DEGs are in WGCNA moduleat p va;ue < 0.05
# 378 genes and 79.4 % this is all of the DEGS -Summary, module membership < 0.05 contains ALL DEGs that overlapped with the module 


# What is the mean ModuleMembership calue? (Pearsons correlation of eigengen velues)
cor_d18_HighSal <- d18_SalHigh_WGCNA  %>%  dplyr::filter(Cvirginica_TranscriptID %in% d18_SalHigh_DEGs$Cvirginica_TranscriptID & 
                                                           moduleColor %in% 'blue' & 
                                                           p.MM.blue < 0.5) %>% 
                                           dplyr::select(MM.blue)

mean(cor_d18_HighSal$MM.blue) # 0.8179135
sd(cor_d18_HighSal$MM.blue)   # 0.1143647

