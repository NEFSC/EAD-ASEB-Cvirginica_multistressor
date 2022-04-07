#Author: Sam Gurr
#Edit by: Sam Gurr
#Date Last Modified: 20210216

rm(list=ls())

# Load packages and pacage version/date/import/depends info
library(dplyr)
library(VennDiagram)
library("ggVennDiagram")
library(ggvenn)
library(gridExtra)
install.packages("ggVennDiagram")
#set working directory--------------------------------------------------------------------------------------------------
# SET WORKING DIRECTORY AND LOAD DATA
setwd("C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/")

# upload data
## main pairwise effect DEGs (salinity, temperature, OA - high vs. low for each)
DE_d2.Temperature <- read.csv(file="C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Data/TagSeq/DESeq2/Day2_larva/Day2.Temperature_DESeq2results.csv", sep=',', header=TRUE)
DE_d2.Salinity    <- read.csv(file="C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Data/TagSeq/DESeq2/Day2_larva/Day2.Salinity_DESeq2results.csv", sep=',', header=TRUE)
DE_d2.OA          <- read.csv(file="C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Data/TagSeq/DESeq2/Day2_larva/Day2.OA_DESeq2results.csv", sep=',', header=TRUE)
colnames(DE_d2.Temperature)[2]  == "TranscriptID" # TRUE
colnames(DE_d2.Salinity)[2]     == "TranscriptID" # TRUE
colnames(DE_d2.OA)[2]           == "TranscriptID" # TRUE

DE_d18.Salinity    <- read.csv(file="C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Data/TagSeq/DESeq2/Day18_spat/Day18.salinity_DESeq2results.csv", sep=',', header=TRUE)
DE_d18.OA          <- read.csv(file="C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Data/TagSeq/DESeq2/Day18_spat/Day18.OA_DESeq2results.csv", sep=',', header=TRUE)
colnames(DE_d18.Salinity)[2]     == "TranscriptID" # TRUE
colnames(DE_d18.OA)[2]           == "TranscriptID" # TRUE

# call only sig DEGs - osolate just the gene names - rename column to day

##############################################################################################################################
# primary Treatmennt A v M   :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
############################################################################################################################## 

# Day2_larva ------------------------------------------------------------------------ #

## TEMPERATURE
DE_d2.Temperature_UP        <- DE_d2.Temperature %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('TranscriptID')
DE_d2.Temperature_UP$Dir    <- "upregulated"
DE_d2.Temperature_UP$time   <- "Day2_larva"
nrow(DE_d2.Temperature_UP) # 359
DE_d2.Temperature_DOWN      <- DE_d2.Temperature %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 1) %>%  dplyr::select('TranscriptID')
DE_d2.Temperature_DOWN$Dir  <- "downregulated"
DE_d2.Temperature_DOWN$time <- "Day2_larva"
nrow(DE_d2.Temperature_DOWN) # 278

## SALINITY
DE_d2.Salinity_UP           <- DE_d2.Salinity %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('TranscriptID')
DE_d2.Salinity_UP$Dir       <- "upregulated"
DE_d2.Salinity_UP$time      <- "Day2_larva"
nrow(DE_d2.Salinity_UP) # 341

DE_d2.Salinity_DOWN         <- DE_d2.Salinity %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 1) %>%  dplyr::select('TranscriptID')
DE_d2.Salinity_DOWN$Dir     <- "downregulated"
DE_d2.Salinity_DOWN$time    <- "Day2_larva"
nrow(DE_d2.Salinity_DOWN) # 378

## OA
DE_d2.OA_UP                 <- DE_d2.OA %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('TranscriptID')
DE_d2.OA_UP$Dir             <- "upregulated"
DE_d2.OA_UP$time            <- "Day2_larva"
nrow(DE_d2.OA_UP) # 82

DE_d2.OA_DOWN               <- DE_d2.OA %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 1) %>%  dplyr::select('TranscriptID')
DE_d2.OA_DOWN$Dir           <- "downregulated"
DE_d2.OA_DOWN$time          <- "Day2_larva"
nrow(DE_d2.OA_DOWN) # 83



# Day18_spat ------------------------------------------------------------------------ #


## SALINITY
DE_d18.Salinity_UP           <- DE_d18.Salinity %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('TranscriptID')
DE_d18.Salinity_UP$Dir       <- "upregulated"
DE_d18.Salinity_UP$time      <- "Day18_spat"
nrow(DE_d18.Salinity_UP) # 483

DE_d18.Salinity_DOWN         <- DE_d18.Salinity %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 1) %>%  dplyr::select('TranscriptID')
DE_d18.Salinity_DOWN$Dir     <- "downregulated"
DE_d18.Salinity_DOWN$time    <- "Day18_spat"
nrow(DE_d18.Salinity_DOWN) # 585

## OA
DE_d18.OA_UP                 <- DE_d18.OA %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 1) %>%  dplyr::select('TranscriptID')
DE_d18.OA_UP$Dir             <- "upregulated"
DE_d18.OA_UP$time            <- "Day18_spat"
nrow(DE_d18.OA_UP) # 148

DE_d18.OA_DOWN               <- DE_d18.OA %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 1) %>%  dplyr::select('TranscriptID')
DE_d18.OA_DOWN$Dir           <- "downregulated"
DE_d18.OA_DOWN$time          <- "Day18_spat"
nrow(DE_d18.OA_DOWN) # 228



# all Primary treatement effects in the GROUP MOD  ------------------------------------------------------------------------ #


# OA time- up and downregulated genes
Day2_larvae_all <- list(
  OA.UP        = DE_d2.OA_UP$TranscriptID, 
  OA.DOWN      = DE_d2.OA_DOWN$TranscriptID, 
  
  Temp.UP      = DE_d2.Temperature_UP$TranscriptID, 
  Temp.DOWN    = DE_d2.Temperature_DOWN$TranscriptID, 
  
  Sal.UP       = DE_d2.Salinity_UP$TranscriptID, 
  Sal.DOWN     = DE_d2.Salinity_DOWN$TranscriptID 
)




# OA time- up and downregulated genes
OA_ALL <- list(
  Day2.UP      = DE_d2.OA_UP$TranscriptID, 
  Day18.UP     = DE_d18.OA_UP$TranscriptID,
  Day2.DOWN    = DE_d2.OA_DOWN$TranscriptID, 
  Day18.DOWN   = DE_d18.OA_DOWN$TranscriptID
)


# Salinity time- up and downregulated genes
Salinity_ALL <- list(
  Day2.UP      = DE_d2.Salinity_UP$TranscriptID, 
  Day18.UP     = DE_d18.Salinity_UP$TranscriptID,
  Day2.DOWN    = DE_d2.Salinity_DOWN$TranscriptID, 
  Day18.DOWN   = DE_d18.Salinity_DOWN$TranscriptID
)


# Venn diagram - OA all time
venn.d_larvae <- ggVennDiagram(Day2_larvae_all, 
                        label_alpha = 2)
venn.d_larvae <- venn.d_larvae + ggtitle("Day2 larvae: all DEGs")




# Venn diagram - OA all time
vennOA <- ggVennDiagram(OA_ALL, 
              label_alpha = 2)
VennOA <- vennOA + ggtitle("OA all DEGs")




# Venn diagram - SALINITY all time
vennSalinity <- ggVennDiagram(Salinity_ALL, 
                        label_alpha = 2)
VennSalinity <- vennSalinity + ggtitle("Salinity all DEGs")




#  grid
grid.arrange(VennOA, VennSalinity, ncol=1, nrow=2, clip="off")

pdf("C:/Users/samuel.gurr/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/DESeq2/Venn_OA_Salinity.pdf")
grid.arrange(VennOA, VennSalinity, ncol=1, nrow=2, clip="off")
dev.off()



# Venn diagram - ggvenn
up.venn2 <- ggvenn(UPREG.primary, 
                   fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                   stroke_size = 0.5, set_name_size = 4)
up.venn2 <- up.venn2 + ggtitle("UPREGULATED: Primary treatment")
down.venn2 <- ggvenn(DWONREG.primary, 
                     fill_color = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     stroke_size = 0.5, set_name_size = 4)
down.venn2 <- down.venn2 + ggtitle("DOWNREGULATED: Primary treatment")
grid.arrange(up.venn2, down.venn2, ncol=2, nrow=1, clip="off")

# with day 0
up.venn2_w0 <- ggvenn(UPREG.primary_withday0, 
                   fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                   stroke_size = 0.5, set_name_size = 4)
up.venn2_w0 <- up.venn2_w0 + ggtitle("UPREGULATED: Primary treatment")
down.venn2_w0 <- ggvenn(DWONREG.primary_withday0, 
                     fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                     stroke_size = 0.5, set_name_size = 4)
down.venn2_w0 <- down.venn2_w0 + ggtitle("DOWNREGULATED: Primary treatment")
Venn_primary <- grid.arrange(up.venn2_w0, down.venn2_w0, ncol=2, nrow=1, clip="off")




# Day 7 sanity check of module design   ------------------------------------------------------------------------ #



UPREG.Day7_primary_allmods <- list(
  Full.Mod_A_v_M                  = full_day7.primaryDE_UP$Gene, 
  GroupTest_AA_v_MA               = grouptest_day7.primaryDE_UP$Gene, 
  Main.Mod_A_v_M                  = main_day7.primaryDE_UP$Gene, 
  Group.Mod_AA.AM.AS_v_MA.MM.MS   = day7.primaryDE_UP$Gene
)
DWNREG.Day7_primary_allmods <- list(
  Full.Mod_A_v_M                  = full_day7.primaryDE_DOWN$Gene, 
  GroupTest_AA_v_MA               = grouptest_day7.primaryDE_DOWN$Gene, 
  Main.Mod_A_v_M                  = main_day7.primaryDE_DOWN$Gene, 
  Group.Mod_AA.AM.AS_v_MA.MM.MS   = day7.primaryDE_DOWN$Gene
)

# Venn Diagram
up.venn2_day7_mods <- ggvenn(UPREG.Day7_primary_allmods, 
                      fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                      stroke_size = 0.5, set_name_size = 4)
up.venn2_day7_mods <- up.venn2_day7_mods + ggtitle("UPREGULATED: Day 7 mod tests (primary treatment)")
down.venn2__day7_mods <- ggvenn(DWNREG.Day7_primary_allmods, 
                        fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                        stroke_size = 0.5, set_name_size = 4)
down.venn2__day7_mods <- down.venn2__day7_mods + ggtitle("DOWNREGULATED: Day 7 mod tests (primary treatment)")
Day7_Venn_primary_tests <- grid.arrange(up.venn2_day7_mods, down.venn2__day7_mods, ncol=2, nrow=1, clip="off")




# Day 14 sanity check of module design   ------------------------------------------------------------------------ #



UPREG.Day14_primary_allmods <- list(
  Full.Mod_A_v_M                  = full_day14.primaryDE_UP$Gene, 
  GroupTest_AA_v_MA               = grouptest_day14.primaryDE_UP$Gene, 
  Main.Mod_A_v_M                  = main_day14.primaryDE_UP$Gene, 
  Group.Mod_AA.AM.AS_v_MA.MM.MS   = day14.primaryDE_UP$Gene
)
DWNREG.Day14_primary_allmods <- list(
  Full.Mod_A_v_M                  = full_day14.primaryDE_DOWN$Gene, 
  GroupTest_AA_v_MA               = grouptest_day14.primaryDE_DOWN$Gene, 
  Main.Mod_A_v_M                  = main_day14.primaryDE_DOWN$Gene, 
  Group.Mod_AA.AM.AS_v_MA.MM.MS   = day14.primaryDE_DOWN$Gene
)

# Venn Diagram
up.venn2_day14_mods <- ggvenn(UPREG.Day14_primary_allmods, 
                             fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                             stroke_size = 0.5, set_name_size = 4)
up.venn2_day14_mods <- up.venn2_day14_mods + ggtitle("UPREGULATED: Day 14 mod tests (primary treatment)")
down.venn2__day14_mods <- ggvenn(DWNREG.Day14_primary_allmods, 
                                fill_color = c("white", "#E69F00", "#56B4E9", "#009E73"),
                                stroke_size = 0.5, set_name_size = 4)
down.venn2__day14_mods <- down.venn2__day14_mods + ggtitle("DOWNREGULATED: Day 14 mod tests (primary treatment)")
Day14_Venn_primary_tests <- grid.arrange(up.venn2_day14_mods, down.venn2__day14_mods, ncol=2, nrow=1, clip="off")





# output ------------------------------------------------------------------------ #




pdf("Analysis/Output/DESeq2/Venn_DE.pdf")
grid.arrange(up.venn2_w0, down.venn2_w0, ncol=2, nrow=1, clip="off")
dev.off()

pdf("Analysis/Output/DESeq2/Day7/D7_res_model_tests/Venn_DE_Day7_test.pdf")
grid.arrange(up.venn2_day7_mods, down.venn2__day7_mods, ncol=1, nrow=2, clip="off")
dev.off()

pdf("Analysis/Output/DESeq2/Day14/D14_res_model_tests/Venn_DE_Day14_test.pdf")
grid.arrange(up.venn2_day14_mods, down.venn2__day14_mods, ncol=1, nrow=2, clip="off")
dev.off()


##############################################################################################################################
# Day 7 Second treatment ( with MA v. AM pairwise )   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
############################################################################################################################## 

# prep data - filter log fold change and call sig DEGs + build new columns (not necessary for the venns done below but whatevs)
day7.secondDE_AvM_UP <- day7.secondDE_AvM %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 0) %>%  dplyr::select('Gene')
day7.secondDE_AvM_UP$Dir <- "upregulated"
day7.secondDE_AvM_UP$Effect <- "Second_AvM"
day7.secondDE_AvM_DOWN <- day7.secondDE_AvM %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 0)  %>%  dplyr::select('Gene')
day7.secondDE_AvM_DOWN$Dir <- "downregulated"
day7.secondDE_AvM_DOWN$Effect <- "Second_AvM"
# group test of the full mod
day7.fullDE_MA_AM_UP <- day7.fullDE_MA_AM %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange > 0) %>%  dplyr::select('Gene')
day7.fullDE_MA_AM_UP$Dir <- "upregulated"
day7.fullDE_MA_AM_UP$Effect <- "MA_vs_AM"
day7.fullDE_MA_AM_DOWN <- day7.fullDE_MA_AM %>% dplyr::filter(padj < 0.05) %>% dplyr::filter(log2FoldChange < 0)  %>%  dplyr::select('Gene')
day7.fullDE_MA_AM_DOWN$Dir <- "downregulated"
day7.fullDE_MA_AM_DOWN$Effect <- "MA_vs_AM"


# create list for Venn to call common and unique gene names assocaited with the up and down reg of the two models...
ListUPREG_Day7 <- list(
  Second_AvM    = day7.secondDE_AvM_UP$Gene, 
  MA_vs_AM      = day7.fullDE_MA_AM_UP$Gene
)
ListDOWNREG_Day7 <- list(
  Second_AvM  = day7.secondDE_AvM_DOWN$Gene, 
  MA_vs_AM    = day7.fullDE_MA_AM_DOWN$Gene
)


# Venn Diagram
Ureg_Day7_main_v_full <- ggvenn(ListUPREG_Day7, 
                              fill_color = c("white", "white"),
                              stroke_size = 0.5, set_name_size = 4)
Ureg_Day7_main_v_full <- Ureg_Day7_main_v_full + ggtitle("UPREGULATED: Day 7 Second treatment mods")
Ureg_Day7_main_v_full
Downreg_Day7_main_v_full <- ggvenn(ListDOWNREG_Day7, 
                                 fill_color = c("white", "white"),
                                 stroke_size = 0.5, set_name_size = 4)
Downreg_Day7_main_v_full <- Downreg_Day7_main_v_full + ggtitle("DOWNREGULATED: Day 7 Second treatment mods")
Day7_main_v_full <- grid.arrange(Ureg_Day7_main_v_full, Downreg_Day7_main_v_full, ncol=2, nrow=1, clip="off")


