---
  # title: "KEGG analysis"
  # author: "Samuel Gurr"
  # date: "May 11, 2022"
---
  
# INFORMATION FOR KEGG IN R FOUND HERE: (http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test)
install.packages("fBasics")
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
d2_WGCNA.data                <- read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership.csv")
d2_WGCNA.data_hightemp       <- read.csv("Output/WGCNA/day2_larvae_hightemp/d2.WGCNA_ModulMembership.csv")
d18_WGCNA.data               <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv")

Ref_Master                   <- read.csv(file = "Data/TagSeq/Seq_details/Seq_Reference_Master.csv",header = T) %>% 
                                        dplyr::rename('TranscriptID' = 'Cvirginica_TranscriptID')

Crass_gigas_genome           <- keggList("crg") # call the C. gigas genome! - notice the csa terms are rownames! (REUIQRES INTERNET)
Crass_gigas_genome_dataframe <- as.data.frame(Crass_gigas_genome) %>%  rownames_to_column() # with will allow us to merge 
colnames(Crass_gigas_genome_dataframe) <- c('sseqid', 'Gene_name')

###############################################################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#   WGCNA DATA: KEGG ANALYSIS OF ONLY GENES THAT ALIGNED WITH OYSTER KEGG DATABASE (BLASTED AGINST THE OYSTER)
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
###############################################################################################################################################

# USING KEGGPROFILE 

#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;:::::::::::::: #
# Day 2 for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;:::::::::::::: #

Day2_WGCNA_sigmodules <- as.data.frame(c('pink', 
                                         'blue', 
                                         'turquoise', 
                                         'brown', 
                                         'black', 
                                         'red', 
                                         'yellow'))
# note: added ModMemership cutoff (Pearson cor and pvalue) 
# when calling 'ModuleLoop_KEGGIDs' - saved to a new directory
for (i in 1:nrow(Day2_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day2_WGCNA_sigmodules[i,1]
  loopmodColor_cor <- paste("MM.", modColor, sep = '')
  loopmodColor_p   <- paste("p.MM.", modColor, sep = '')
  
  # call the module color in the Day 7 data
  ModuleLoop                     <- d2_WGCNA.data %>% dplyr::filter(moduleColor %in% modColor)  %>% 
                                      dplyr::select(c('TranscriptID', 'KEGG_ID', 'geneSymbol', moduleColor, loopmodColor_p, loopmodColor_cor)) 
  genes_per_module               <- length(unique(ModuleLoop$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module          <- nrow(ModuleLoop %>% dplyr::filter(!KEGG_ID %in% NA)) # nrow(ModuleLoop) # use this for the looped print out 
  

  genes_per_module_blasthit      <- na.omit(ModuleLoop)  %>% # ommit  genes without gene name annotation 
                                    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
                                    dplyr::filter(TranscriptID %in% Ref_Master$TranscriptID) %>%  # call only genes that have a blast hit      
                                    nrow() # 
  
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  # calaculate the percent mapped and print this...
  print(paste("Day2", modColor, " ", 
              genes_per_module, "unique genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  # Run KEGG analysis
  ModuleLoop_KEGGIDs       <- as.data.frame(ModuleLoop %>% 
                                              dplyr::filter(.[[5]] < 0.05 & .[[6]] > 0.6) %>% # filter based on thresholds set in the WGCNA + DESeq2 overlap w/PCA 
                                              dplyr::select(c('TranscriptID', 'KEGG_ID')) %>% 
                                              na.omit() %>% 
                                              dplyr::mutate(KEGG_ID = str_split(KEGG_ID,";")) %>% 
                                              unnest(KEGG_ID) %>% 
                                              dplyr::select(-'TranscriptID'))
  
  KEGG_vector_cvirg_Cgigas <- as.vector(gsub(".*:","", ModuleLoop_KEGGIDs$KEGG_ID)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas                     <- enrichKEGG(gene = KEGG_vector_cvirg_Cgigas, 
                                                organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                                                pvalueCutoff = 0.05) 
        # if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
        if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
          # creat dateframe and write the csv file out 
          df                      <- as.data.frame(head(KEGG_cgigas))
          rownames(df)            <- c()
          KEGGoutput              <- as.data.frame(do.call(cbind.data.frame, df))
          KEGGoutput$GeneRatio_2  <- gsub("/"," of ", KEGGoutput$GeneRatio)
          KEGGoutput$Rich_Factor  <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
          
          write.csv(KEGGoutput, file = paste("Output/WGCNA/day2_larvae/KEGG/MM_cutoff/Day2_",modColor,"_KEGG_allgenes.csv", sep ='')) 
          
          # Plot
          plot<- KEGGoutput %>%  
            ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
            geom_point( aes(col=qvalue, size=Count)) +   # Draw points
            geom_segment(aes(x=Description, 
                             xend=Description, 
                             y=min(Rich_Factor), 
                             yend=max(Rich_Factor)),  
                         linetype=NA, 
                         size=0) +   # Draw dashed lines
            labs(title="Day 2", 
                 x = "Pathway",
                 y = "Rich Factor",
                 subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
            coord_flip()
           pdf(paste("Output/WGCNA/day2_larvae/KEGG/MM_cutoff/Day2_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
           print(plot)
           dev.off()
    
    
            # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
            df_2                 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
            df_2$gene_IDs        <- as.vector(strsplit(as.character(df_2$geneID), "/"))
            colnames(df_2)       <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
            df_3                 <- unnest(df_2, gene_IDs)
            df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
            Crass_gigas_ref      <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
            df_final             <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
            write.csv(df_final, file = paste("Output/WGCNA/day2_larvae/KEGG/MM_cutoff/Day2_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
            } else {}
            
    print(paste("Finished! Day2 module = ", modColor, sep = " "))
}




#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;:::::::::::::: #
# Day 2 for loop HIGH TEMP ONLY :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;; #
#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;:::::::::::::: #
View(d2_WGCNA.data_hightemp)
Day2_hightemp_WGCNA_sigmodules <- as.data.frame(c('red', 
                                                  'greenyellow', 
                                                  'brown', 
                                                  'pink', 
                                                  'magenta'))
for (i in 1:nrow(Day2_hightemp_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day2_hightemp_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop                <- d2_WGCNA.data_hightemp %>% dplyr::filter(moduleColor %in% modColor)
  genes_per_module          <- length(unique(ModuleLoop$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module     <- nrow(ModuleLoop %>% dplyr::filter(!KEGG_ID %in% NA)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  genes_per_module_blasthit <- na.omit(ModuleLoop)  %>% # ommit  genes without gene name annotation 
                                  dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
                                  dplyr::filter(TranscriptID %in% Ref_Master$TranscriptID) %>%  # call only genes that have a blast hit      
                                  nrow() # 
  
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  # calaculate the percent mapped and print this...
  print(paste("Day2 (high temp only)", modColor, " ", 
              genes_per_module, "unique genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  # Run KEGG analysis
  ModuleLoop_KEGGIDs       <- as.data.frame(ModuleLoop %>% 
                                            dplyr::select(c('TranscriptID', 'KEGG_ID')) %>% 
                                            na.omit() %>% 
                                            dplyr::mutate(KEGG_ID = str_split(KEGG_ID,";")) %>% 
                                            unnest(KEGG_ID) %>% 
                                            dplyr::select(-'TranscriptID'))
  
  KEGG_vector_cvirg_Cgigas <- as.vector(gsub(".*:","", ModuleLoop_KEGGIDs$KEGG_ID)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas              <- enrichKEGG(gene = KEGG_vector_cvirg_Cgigas, 
                                        organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                                        pvalueCutoff = 0.05) 
  # if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
  if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
    # creat dateframe and write the csv file out 
    df                      <- as.data.frame(head(KEGG_cgigas))
    rownames(df)            <- c()
    KEGGoutput              <- as.data.frame(do.call(cbind.data.frame, df))
    KEGGoutput$GeneRatio_2  <- gsub("/"," of ", KEGGoutput$GeneRatio)
    KEGGoutput$Rich_Factor  <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
    
    write.csv(KEGGoutput, file = paste("Output/WGCNA/day2_larvae_hightemp/KEGG/Day2_",modColor,"_KEGG_allgenes.csv", sep ='')) 
    
    # Plot
    plot<- KEGGoutput %>%  
      ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
      geom_point( aes(col=qvalue, size=Count)) +   # Draw points
      geom_segment(aes(x=Description, 
                       xend=Description, 
                       y=min(Rich_Factor), 
                       yend=max(Rich_Factor)),  
                   linetype=NA, 
                   size=0) +   # Draw dashed lines
      labs(title="Day 2 (high temp only)", 
           x = "Pathway",
           y = "Rich Factor",
           subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
      coord_flip()
    pdf(paste("Output/WGCNA/day2_larvae_hightemp/KEGG/Day2_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
    print(plot)
    dev.off()
    
    
    # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
    df_2                 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
    df_2$gene_IDs        <- as.vector(strsplit(as.character(df_2$geneID), "/"))
    colnames(df_2)       <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
    df_3                 <- unnest(df_2, gene_IDs)
    df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
    Crass_gigas_ref      <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
    df_final             <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
    write.csv(df_final, file = paste("Output/WGCNA/day2_larvae_hightemp/KEGG/Day2_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
    
  } else {}
  
  print(paste("Finished! Day2 (high temp only) module = ", modColor, sep = " "))
}








#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;:::::::::::::: #
# Day 18  for loop ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::; #
#  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;;;:::::::::::::: #

Day18_WGCNA_sigmodules <- as.data.frame(c('tan', 
                                          'red', 
                                          'turquoise', 
                                          'salmon', 
                                          'blue',
                                          'green',
                                          'greenyellow'))

for (i in 1:nrow(Day18_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day18_WGCNA_sigmodules[i,1]
  loopmodColor_cor <- paste("MM.", modColor, sep = '')
  loopmodColor_p   <- paste("p.MM.", modColor, sep = '')
  
  # call the module color in the Day 7 data
  ModuleLoop                     <- d18_WGCNA.data %>% dplyr::filter(moduleColor %in% modColor)  %>% 
    dplyr::select(c('TranscriptID', 'KEGG_ID', 'geneSymbol', moduleColor, loopmodColor_p, loopmodColor_cor)) 
  genes_per_module               <- length(unique(ModuleLoop$geneSymbol)) # nrow(ModuleLoop) # use this for the looped print out 
  annotgenes_per_module          <- nrow(ModuleLoop %>% dplyr::filter(!KEGG_ID %in% NA)) # nrow(ModuleLoop) # use this for the looped print out 
  
  
  genes_per_module_blasthit      <- na.omit(ModuleLoop)  %>% # ommit  genes without gene name annotation 
    dplyr::filter(moduleColor %in% modColor) %>%   # filter for the module loop
    dplyr::filter(TranscriptID %in% Ref_Master$TranscriptID) %>%  # call only genes that have a blast hit      
    nrow() # 
  
  perc_annot_genes_with_blasthit <- ( genes_per_module_blasthit / annotgenes_per_module ) *100
  # calaculate the percent mapped and print this...
  print(paste("Day18", modColor, " ", 
              genes_per_module, "unique genes per module", 
              annotgenes_per_module, " annotated; ", 
              genes_per_module_blasthit, " or", perc_annot_genes_with_blasthit,"% annotated genes with blasthit", sep = ' '))
  
  # Run KEGG analysis
  ModuleLoop_KEGGIDs       <- as.data.frame(ModuleLoop %>% 
                                              dplyr::filter(.[[5]] < 0.05 & .[[6]] > 0.6) %>%
                                              dplyr::select(c('TranscriptID', 'KEGG_ID')) %>% 
                                              na.omit() %>% 
                                              dplyr::mutate(KEGG_ID = str_split(KEGG_ID,";")) %>% 
                                              unnest(KEGG_ID) %>% 
                                              dplyr::select(-'TranscriptID'))
  
  KEGG_vector_cvirg_Cgigas <- as.vector(gsub(".*:","", ModuleLoop_KEGGIDs$KEGG_ID)) # ommit the 'crg:' before the actual terms
  KEGG_cgigas <- enrichKEGG(gene = KEGG_vector_cvirg_Cgigas, 
                            organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                            pvalueCutoff = 0.05) 
          # if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
          if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
            # creat dateframe and write the csv file out 
            df <- as.data.frame(head(KEGG_cgigas))
            rownames(df) <- c()
            KEGGoutput <- as.data.frame(do.call(cbind.data.frame, df))
            KEGGoutput$GeneRatio_18 <- gsub("/"," of ", KEGGoutput$GeneRatio)
            KEGGoutput$Rich_Factor <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
            
            write.csv(KEGGoutput, file = paste("Output/WGCNA/day18_spat/KEGG/RR_cutoff/Day18_",modColor,"_KEGG_allgenes.csv", sep ='')) 
            
            # Plot
            plot<- KEGGoutput %>%  
              ggplot(aes(x=reorder(Description, Rich_Factor), y= Rich_Factor)) + 
              geom_point( aes(col=qvalue, size=Count)) +   # Draw points
              geom_segment(aes(x=Description, 
                               xend=Description, 
                               y=min(Rich_Factor), 
                               yend=max(Rich_Factor)),  
                           linetype=NA, 
                           size=0) +   # Draw dashed lines
              labs(title="Day 18", 
                   x = "Pathway",
                   y = "Rich Factor",
                   subtitle=paste("WGCNA Module:", modColor, sep =' ')) +
              coord_flip()
            pdf(paste("Output/WGCNA/day18_spat/KEGG/RR_cutoff/Day18_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
            print(plot)
            dev.off()
            
            
            # stringsplit and unnest for a data set of genes and IDs associated with each pathway 
            df_18                 <- as.data.frame(KEGG_cgigas)[c(1:2,8:9)]
            df_18$gene_IDs        <- as.vector(strsplit(as.character(df_18$geneID), "/"))
            colnames(df_18)       <- c("Cgigas_PathwayID", "Pathway_Description", "Pathway_gene_list", "Pathway_gene_count", "gene_IDs")
            df_3                 <- unnest(df_18, gene_IDs)
            df_3$Cgigas_KEGG_IDs <- paste("crg:", df_3$gene_IDs, sep='')
            Crass_gigas_ref      <- Crass_gigas_genome_dataframe %>% mutate(Cgigas_KEGG_IDs = Crass_gigas_genome_dataframe$sseqid) %>% select(c('Cgigas_KEGG_IDs','Gene_name'))
            df_final             <- merge(df_3, Crass_gigas_ref, by='Cgigas_KEGG_IDs')
            write.csv(df_final, file = paste("Output/WGCNA/day18_spat/KEGG/RR_cutoff/Day18_",modColor,"_KEGG_allgenes_unlisted.csv", sep ='')) 
            
          } else {}
  
  print(paste("Finished! Day18 module = ", modColor, sep = " "))
}




