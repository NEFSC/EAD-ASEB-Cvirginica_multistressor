---
  # title: "KEGG analysis"
  # author: "Samuel Gurr"
  # date: "May 2, 2022"
---
  
# INFORMATION FOR KEGG IN R FOUND HERE: (http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test)

# LOAD PACKAGES
library(reactome.db)
library(clusterProfiler)
library(KEGGREST)
library(tidyr)
library(stringr)
library(forcats)
library(ggplot2)
library(scales)
library(ape)
library(data.table)
library(tidyverse)
library(fBasics)
library(dplyr)
library(KEGGREST)
library(ggplot2)
# SET WORKING DIRECTORY   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #

# WGCNA results (all treatments)
d2_WGCNA     <- read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv")
d18_WGCNA    <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv")

Ref_Master   <- read.csv(file = "Data/TagSeq/Seq_details/Seq_Reference_Master.csv",header = T) %>% 
                    dplyr::rename('TranscriptID' = 'Cvirginica_TranscriptID')

###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   WGCNA DATA: KEGG ANALYSIS OF ONLY GENES THAT ALIGNED WITH OYSTER KEGG DATABASE (BLASTED AGINST THE OYSTER)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################

# USING KEGGPROFILE -- NOTE: this is outdated and marked as so... 

module_blue           <- d2_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% 'blue' & p.MM.blue < 0.05)
KEGG_vector_mod.blue  <- as.vector(gsub(".*:","",module_blue$Cgigas_KEGGID)) # ommit the 'crg:' before the actual terms
KEGG_mod.blue         <- enrichKEGG(gene = KEGG_vector_mod.blue, 
                          organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                          pvalueCutoff = 0.05) 

module_brown           <- d2_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% 'brown' & p.MM.brown < 0.05)
KEGG_vector_mod.brown  <- as.vector(gsub(".*:","",module_brown$Cgigas_KEGGID)) # ommit the 'crg:' before the actual terms
KEGG_mod.brown         <- enrichKEGG(gene = KEGG_vector_mod.brown, 
                                    organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                                    pvalueCutoff = 0.05) 

# if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...

  # creat dateframe and write the csv file out 
  df <- as.data.frame(head(KEGG_mod.brown))
  rownames(df) <- c()
  KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
  KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
  KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
  
  write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_treatments_all/Day7_frontloaded_severe_KEGG_allgenes.csv", sep ='')) 
  
  # Plot
  plot <- KEGGoutput %>%  
            ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
            geom_point( aes(col=qvalue, size=Count)) +   # Draw points
            geom_segment(aes(x=Description, 
                             xend=Description, 
                             y=min(Rich_Factor), 
                             yend=max(Rich_Factor)),  
                         linetype=NA, 
                         size=0) +   # Draw dashed lines
            labs(title="Day 2 larvae", 
                 x = "Pathway",
                 y = "Rich Factor",
                 subtitle="module blue") +
            theme_bw() +
            coord_flip()
  #pdf(paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day0_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
  print(plot)
  dev.off()
  
  
# stringsplit and unnest for a data set of genes and IDs associated with each pathway 
  df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
  df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
  colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
  df_3 <- unnest(df_2, gene_IDs)
  df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
  Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
  df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
#write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day0_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 

  
  
  
  




# Day 7 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day7_WGCNA_sigmodules <- as.data.frame(c('brown','yellow', 'green'))
for (i in 1:nrow(Day7_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day7_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  module_without_filter   <- d7_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module        <- length(unique(module_without_filter$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module   <- length(na.omit(module_without_filter$HGNC)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  ModuleLoop_blasthit            <- d7_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module_blasthit      <- na.omit(module_without_filter)  %>% # ommit  genes without gene name annotation 
                                    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
                                    dplyr::filter(geneSymbol %in% d7_WGCNA_crgKEGGhits$geneSymbol) %>%  # call only genes that have a blast hit      
                                    nrow() # 
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  # calaculate the percent mapped and print this...
  print(paste("Day7", modColor, " ", 
              genes_per_module, "genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",ModuleLoop$crg_KO)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day7_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 7", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
     pdf(paste("Analysis/Output/KEGG/WGCNA/Day7_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
     print(plot)
     dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day7_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day7 module = ", modColor, sep = " "))
}





# Day 14 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day14_WGCNA_sigmodules <- as.data.frame(c('brown','black', 'magenta', 'pink'))
for (i in 1:nrow(Day14_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day14_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  module_without_filter   <- d14_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module        <- length(unique(module_without_filter$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module   <- length(na.omit(module_without_filter$HGNC)) # nrow(ModuleLoop) # use this for the looped print out 
  
  ModuleLoop_blasthit            <- d14_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module_blasthit      <- na.omit(module_without_filter)  %>% # ommit  genes without gene name annotation 
    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
    dplyr::filter(geneSymbol %in% d14_WGCNA_crgKEGGhits$geneSymbol) %>%  # call only genes that have a blast hit      
    nrow() # 
  
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  
  # calaculate the percent mapped and print this...
  print(paste("Day14", modColor, " ", 
              genes_per_module, "genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",ModuleLoop_blasthit$crg_KO)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05)
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day14_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 14", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day14_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day14_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day14 module = ", modColor, sep = " "))
}




# Day 21 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day21_WGCNA_sigmodules <- as.data.frame(c('blue','magenta', 'yellow', 'black', 'pink', 'red', 'turquoise'))
for (i in 1:nrow(Day21_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day21_WGCNA_sigmodules[3,1]
  
  # call the module color in the Day 7 data
  module_without_filter   <- d21_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module        <- length(unique(module_without_filter$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module   <- length(na.omit(module_without_filter$HGNC)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  ModuleLoop_blasthit            <- d21_WGCNA_crgKEGGhits %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module_blasthit      <- na.omit(module_without_filter)  %>% # ommit  genes without gene name annotation 
    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
    dplyr::filter(geneSymbol %in% d21_WGCNA_crgKEGGhits$geneSymbol) %>%  # call only genes that have a blast hit      
    nrow() # 
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  
  # calaculate the percent mapped and print this...
  print(paste("Day21", modColor, " ", 
              genes_per_module, "genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",ModuleLoop_blasthit$crg_KO)) # omit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05)
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 21", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    
    
    pdf(paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/subseq_treatments_all/WGCNA/Day21_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day21 module = ", modColor, sep = " "))
}
























###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   BETTER PLOTS FOR THE KEGG ENRICHMENT ANALYSIS (rich factor plots) 
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################


# LOAD DATA    :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
# AMBIENT EFFECT MODUELS 
d7_KEGG_brown          <- read.csv("Analysis/Output/KEGG/WGCNA/Day7_brown_KEGG_allgenes.csv")
d7_KEGG_brown$ModDay   <- "Day7_brown"
d14_KEGG_brown         <- read.csv("Analysis/Output/KEGG/WGCNA/Day14_brown_KEGG_allgenes.csv")
d14_KEGG_brown$ModDay  <- "Day14_brown"
d21_KEGG_blue          <- read.csv("Analysis/Output/KEGG/WGCNA/Day21_blue_KEGG_allgenes.csv")
d21_KEGG_blue$ModDay   <- "Day21_blue"
# SET-UP AND PLOT  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_ambienteffect_mods        <- rbind(d7_KEGG_brown, d14_KEGG_brown, d21_KEGG_blue)
KEGG_ambienteffect_mods$ModDay <- factor(KEGG_ambienteffect_mods$ModDay , levels = c("Day7_brown", "Day14_brown", "Day21_blue")) # for the correct order of facets in the plot below
Ambient_Eff_Mods_WGCNA <- KEGG_ambienteffect_mods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Ambient-Effect Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
Ambient_Eff_Mods_WGCNA


# AMBIENT EFFECT MODUELS 
d7_KEGG_yelow            <- read.csv("Analysis/Output/KEGG/WGCNA/Day7_yellow_KEGG_allgenes.csv")
d7_KEGG_yelow$ModDay     <- "Day7_yellow"
d14_KEGG_black           <- read.csv("Analysis/Output/KEGG/WGCNA/Day14_black_KEGG_allgenes.csv")
d14_KEGG_black$ModDay    <- "Day14_black"
d21_KEGG_yellow          <- read.csv("Analysis/Output/KEGG/WGCNA/Day21_yellow_KEGG_allgenes.csv")
d21_KEGG_yellow$ModDay   <- "Day21_yellow"
# SET-UP AND PLOT  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_moderateeffect_mods        <- rbind(d7_KEGG_yelow, d14_KEGG_black, d21_KEGG_yellow)
KEGG_moderateeffect_mods$ModDay <- factor(KEGG_moderateeffect_mods$ModDay , levels = c("Day7_yellow", "Day14_black", "Day21_yellow")) # for the correct order of facets in the plot below
Moderate_Eff_Mods_WGCNA <- KEGG_moderateeffect_mods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Moderate-Effect Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
Moderate_Eff_Mods_WGCNA


# ALL MODULES PRIMARY EFFECT :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
KEGG_ALLprimaryeffectmods      <- rbind(d7_KEGG_brown, d14_KEGG_brown, d21_KEGG_blue,d7_KEGG_yelow, d14_KEGG_black, d21_KEGG_yellow)
KEGG_ALLprimaryeffectmods$ModDay <- factor(KEGG_ALLprimaryeffectmods$ModDay , levels = c("Day7_brown", "Day14_brown", "Day21_blue","Day7_yellow", "Day14_black", "Day21_yellow")) # for the correct order of facets in the plot below
PrimaryEff_KEGG_WGCNA <- KEGG_ALLprimaryeffectmods  %>%  
  ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
  geom_point( aes(col=qvalue, size=Count)) +   # Draw points
  scale_colour_gradient(low = "#D55E00", high = "#56B4E9", limits = c(0,0.05)) +
  geom_segment(aes(x=Description, 
                   xend=Description, 
                   y=min(Rich_Factor), 
                   yend=max(Rich_Factor)),  
               linetype=NA, 
               size=0) +   # Draw dashed lines
  labs(title="KEGG pathway enrichment analysis", 
       x = "Pathway",
       y = "Rich Factor",
       subtitle="Primaryt-Effect Modules (WGCNA)") +
  theme_bw() +
  coord_flip() +
  facet_wrap(~ModDay,scales="free_y")
PrimaryEff_KEGG_WGCNA


pdf(paste("Analysis/Output/KEGG/WGCNA/PrimaryEffectModules_RichFactor.pdf", sep =''), width=12, height=10)
print(PrimaryEff_KEGG_WGCNA)
dev.off()








###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   DESeq2 DATA: KEGG ANALYSIS OF ALL GENES (upregualted adn downregulated) FROM PAIRWISE DIFFERENTIAL GENE EXPRESSION ANALYSIS
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################


# For loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
days <- c('Day7_Primary_AvM','Day7_Primary_AvM', 'Day14_Primary_AvM', 'Day14_Primary_AvM', 'Day21_Primary_AvM', 'Day21_Primary_AvM')
DE_dir <- c('up','down', 'up','down', 'up','down')
PrimEff_loop <- data.frame(days, DE_dir)
PrimEff_loop

for (i in 1:nrow(PrimEff_loop)) {
  day_trmt <- PrimEff_loop[i,1]
  dir      <- PrimEff_loop[i,2]
  
  DEGs_loop <- DESeq2_PrimaryEffects_master %>%  dplyr::filter(Day_Trmt == day_trmt, DE_dir ==dir) # note this dataset was already sorted by descending LFC
  total_DEGs <- nrow(DEGs_loop)
  colnames(DEGs_loop)[2] <- 'Pgenerosa_Gene_IDs'
 
  # merge with the best blast hits to get  KO identifiers 
  DEGs_loop_2   <- merge(DEGs_loop, crgKEGG_PgenREF_besthits, by.x=c('Pgenerosa_Gene_IDs'),  by.y=c('geneSymbol'))
  KO_per_effect <- nrow(DEGs_loop_2) # use this for the looped print out 
  
  # calaculate the percent mapped and print this...
  Percent_mapped <- (KO_per_effect / total_DEGs) * 100
  print(paste(day_trmt, " ", dir, ":", Percent_mapped, "percent mapped to Cgigas", "(or ", KO_per_effect, " genes) ", total_DEGs, "total DEGs",  sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",DEGs_loop_2$crg_KO)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio <- gsub("/"," of ", KEGGoutput$GeneRatio)
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/DESeq2/",unlist(strsplit(day_trmt, "_"))[1],"_", dir,"regulated_KEGG.csv", sep ='')) 
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/DESeq2/",unlist(strsplit(day_trmt, "_"))[1],"_", dir,"regulated_KEGG_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste(unlist(strsplit(day_trmt, "_"))[1], "_",dir,"regulated - COMPLETED!", sep = ""))
}






###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   MERGE BY GENE NAME AND RUN THE KEGG FOR WGCNA DATA: - NOTE: this is NOT the bet approach, we did NOT use this for the MS!!!!!!!!!!!!!!!
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################



# PREPARE THE KEGG  Crassostra gigas (Pacific oyster) genome 'Crass_gigas_genome_dataframe' will be called in the KEGG analysis for loops

# condence the ref annotatation file for fewer collumns and edit what we need
annot.condenced <- Geoduck_annotation[,c(1,7)] # load just the PGEN ID and the putative gene terms
annot.condenced$Gene_term    <- sub(" \\(EC.*", "", annot.condenced$V7)  # call the gene term BEFORe the EC is listed
Pgen_reference <- na.omit(annot.condenced[c(1,3)]) # ommit unannotated genes
names(Pgen_reference)  <- c('Pgenerosa_Gene_IDs', 'Gene_terms') # rename the columns to merge with the modules below...



# to cal the gene ID We can use gsub. 
# Match the pattern of one or more characters that are not a - ([^-]+) from the start (^) of the string followed by a - or (| a - followed by characters (.*) and replace it with blank ("")
ref_update_20210602$Pgenerosa_Gene_IDs  <- substr(ref_update_20210602$attributes, 4, 18)
ref_update_20210602$gene_notes  <-  sub(".*Notes//", "", ref_update_20210602$attributes)
ref_update_condenced <- ref_update_20210602[,c(9,10)]
test <- merge(Pgen_reference, ref_update_condenced, by = 'Pgenerosa_Gene_IDs')
View(test) # 34940


# Notes on this pseudo-cluster....
# this cluster was writteni in May  2021 to tackle the KEGG analysis of ALL genes in significant WGCNA modules

# Day 7 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day7_WGCNA_sigmodules <- as.data.frame(c('brown','yellow', 'green'))
for (i in 1:nrow(Day7_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day7_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- d7_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module <- nrow(ModuleLoop) # use this for the looped print out 
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(ModuleLoop$geneSymbol) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- length(unique(EnrichedGenes_2$Gene_terms)) # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset (unique genes)
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(EnrichedGenes_2$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 243 calls
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # calaculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste("Day7", modColor, " ", genes_per_module, " genes per module -->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/WGCNA/Day7_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 7", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Analysis/Output/KEGG/WGCNA/Day7_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/WGCNA/Day7_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day7 module = ", modColor, sep = " "))
}





# Day 14 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day14_WGCNA_sigmodules <- as.data.frame(c('brown','black', 'magenta', 'pink'))
for (i in 1:nrow(Day14_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day14_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- d14_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  colnames(ModuleLoop)
  genes_per_module <- nrow(ModuleLoop) # use this for the looped print out 
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(ModuleLoop$geneSymbol) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- length(unique(EnrichedGenes_2$Gene_terms))  # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  vector                 <- as.vector(EnrichedGenes_2$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 243 calls
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # calculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste("Day14", modColor, " ", genes_per_module, " genes per module -->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  
  
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/WGCNA/Day14_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 14", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Analysis/Output/KEGG/WGCNA/Day14_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/WGCNA/Day14_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day14 module = ", modColor, sep = " "))
}




# Day 21 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
Day21_WGCNA_sigmodules <- as.data.frame(c('blue','magenta', 'yellow', 'black', 'pink', 'red', 'turquoise'))
for (i in 1:nrow(Day21_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day21_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- d21_WGCNA_all %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module <- nrow(ModuleLoop) # use this for the looped print out 
  
  # Prep the enriched genes of interest wihtin each module...
  #NOTE: the master file contains PGEN gene IDs separated by ';' - each row of PGEN IDs contains all genes associated with signifcant GO terms following  GO enrichment (in 'goseq')
  # thus, these genes are called as important due to their GO annotation and significant calls  from functional enrichment analysis 
  EnrichedGenes <- as.data.frame(ModuleLoop$geneSymbol) # split the genes for the significant goseq enrichment 
  colnames(EnrichedGenes)[1] <- 'Pgenerosa_Gene_IDs' # rename the single column in this call
  EnrichedGenes_2 <- merge(Pgen_reference, EnrichedGenes, by ='Pgenerosa_Gene_IDs') # merge with the Pgen reference (Gene terms before the EC is stated!)
  EnrichedGenes_2$Gene_terms    <- sub(" \\(.*", "", EnrichedGenes_2$Gene_terms) # removes all of gene term before the parenthesis 
  EnrichedGenes_2$Gene_terms    <- sub(" \\[.*", "", EnrichedGenes_2$Gene_terms) # calls al of the string BEFORE an open bracket (in some cases the sttring has just an open bracket without the closed 
  EnrichedGenes_2$Gene_terms    <- gsub("\\[.*?\\]", "", EnrichedGenes_2$Gene_terms) # now removes all bracketted terms 
  EnrichedGenes_2$Gene_terms    <- tolower(EnrichedGenes_2$Gene_terms) # convert all to lower case to merge with Cgigas
  EnrichedGenes_2 <- unique(EnrichedGenes_2) # call only unique occurances 
  
  Pgenerosa_gene_calls <- length(unique(EnrichedGenes_2$Gene_terms))  # total genes called - use this to learn the percent merged with the C.gigas KEGG dataset
  
  
  # merge and calculate the percent mapping
  # Lets merge with the Cgigas genome! 
  # EnrichedGenes_2$Gene_terms <- gsub("[[]", "", EnrichedGenes_2$Gene_terms) # need to run this for turquoise - an open bracket [ in one of the genes
  vector                 <- as.vector(EnrichedGenes_2$Gene_terms) # call the gene term as a vector 
  vector_with_asserts    <- paste("^", vector, "$", sep='') # assets help to call the exact beginning and end of the gene term in grep1 (exact mtatch!) - otheriwse grep will use an 'in it' like funciton and call any genes that contain the whole term in addition to extra string chracters
  
  EnrichedGenes_Cgigas <- subset(Crass_gigas_genome_dataframe, grepl(paste("^",vector_with_asserts,"$", sep="", collapse= "|"), Gene_terms, ignore.case = TRUE)) # subset the C gigas genome by the common and exact match to gene terms (with sig GO enrichmetn from 'goseq'!) in the module of interest 
  length(unique(EnrichedGenes_Cgigas$Gene_terms)) # 243 calls
  EnrichedGenes_Cgigas_removeduplicates = EnrichedGenes_Cgigas[!duplicated(EnrichedGenes_Cgigas$Gene_terms),]
  
  Cgigas_gene_calls <- nrow(EnrichedGenes_Cgigas_removeduplicates) # total Pgen gens mapped to the Cgigas KEGG genome  - use this to learn the percent merged with the C.gigas KEGG dataset
  
  # calculate the percent mapped and print this...
  Percent_mapped <- (Cgigas_gene_calls / Pgenerosa_gene_calls) * 100
  print(paste("Day21", modColor, " ", genes_per_module, " genes per module -->", Percent_mapped, "percent mapped to Cgigas KEGG of N =", Pgenerosa_gene_calls, "total genes", sep = ' '))
  
  # Run KEGG analysis
  KEGG_vector_Pgen_Cgigas   <- as.vector(gsub(".*:","",EnrichedGenes_Cgigas_removeduplicates$Cgigas_KEGG_IDs)) # ommit the 'crg:' before the actual terms
  
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_Pgen_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
  # if lloop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df <- as.data.frame(head(KEGG_cgigas))
    rownames(df) <- c()
    KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2 <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Analysis/Output/KEGG/WGCNA/Day21_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    theme_set(theme_classic())
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 14", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Analysis/Output/KEGG/WGCNA/Day21_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2) <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    df_final <- merge(df_3, Crass_gigas_genome_dataframe, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Analysis/Output/KEGG/Day21_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day21 module = ", modColor, sep = " "))
}















