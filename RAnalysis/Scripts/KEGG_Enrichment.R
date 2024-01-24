---
  # title: "KEGG analysis"
  # author: "Samuel Gurr"
  # date: "May 11, 2022"
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
d2_WGCNA.data                <- read.csv("Output/WGCNA/day2_larvae/d2.WGCNA_ModulMembership_RRcutoff.csv")
d2_WGCNA.data_hightemp       <- read.csv("Output/WGCNA/day2_larvae_hightemp/d2.WGCNA_ModulMembership.csv")
d18_WGCNA.data               <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership.csv")
d18_WGCNA.data               <- read.csv("Output/WGCNA/day18_spat/d18.WGCNA_ModulMembership_RRcutoff.csv")


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

# Using KEGGREST

# KEGGREST prep 
pathways.list <- keggList("pathway", "crg")
head(pathways.list)
# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list)) 
genes.by.pathway <- sapply(pathway.codes,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(TRUE,FALSE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           }
)
head(genes.by.pathway)

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
# prep loop for cumulative output table 
df_total            <- data.frame() # start dataframe 
KEGG.Day2           <- data.frame(matrix(nrow = 1, ncol = 9)) # create dataframe to save cumunalitively during for loop
colnames(KEGG.Day2) <- c('Day', 'modColor', 'KEGGID_pathway', 'pathway.name' , 
                         'Num.genes.all', 'Num.genes.exp', 'Gene.IDs', 'p.value', 'log10_pvalue') # names for comuns in the for loop

pathtable <- as.data.frame(pathways.list) %>% 
                dplyr::mutate(pathname = sapply(strsplit(pathways.list, " - Crassostrea"), "[",1)) %>% 
                tibble::rownames_to_column("crg_code")

for (i in 1:nrow(Day2_WGCNA_sigmodules)) {
    modColor   <- Day2_WGCNA_sigmodules[i,1]
    ModuleLoop <- as.data.frame(d2_WGCNA.data %>% 
                                  dplyr::filter(moduleColor %in% modColor)  %>% 
                                  dplyr::select(c('TranscriptID', 'KEGG_ID', 'geneSymbol', 'MM.p', 'MM.cor')) %>% 
                                  dplyr::filter(.[[4]] < 0.05 & .[[5]] > 0.6) %>% # filter based on thresholds set in the WGCNA + DESeq2 overlap w/PCA 
                                  dplyr::select(c('TranscriptID', 'KEGG_ID', 'MM.p')) %>% 
                                  na.omit() %>% 
                                  dplyr::mutate(KEGG_ID = gsub(".*:","",KEGG_ID)) %>% 
                                  unnest(KEGG_ID))
    geneList <- ModuleLoop[,3]
    names(geneList) <- ModuleLoop$KEGG_ID
    # Wilcoxon test for each pathway
    # length(genes.by.pathway$crg04814)
    pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                                 function(pathway) {
                                   pathway.genes         <- genes.by.pathway[[pathway]]
                                   list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                                   list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                                   scores.in.pathway <- geneList[list.genes.in.pathway]
                                   scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                                   if (length(scores.in.pathway) > 0){
                                     p.value <- wilcox.test(scores.in.pathway, 
                                                            (length(pathway.genes) - scores.in.pathway),
                                                            #scores.not.in.pathway,
                                                            alternative = "less")$p.value
                                   } else{
                                     p.value <- NA
                                   }
                                   return(c(p.value = p.value, Num_total_in_pathway = length(pathway.genes), Annotated = length(list.genes.in.pathway), GeneIDs = list(list.genes.in.pathway) ))
                                 }
    ))
    
    # Assemble output table
    outdat <- data.frame(pathway.code = rownames(pVals.by.pathway)) %>% 
                dplyr::mutate(Num.genes.all = as.data.frame(pVals.by.pathway)$Num_total_in_pathway) %>%
                dplyr::mutate(Num.genes.exp = as.data.frame(pVals.by.pathway)$Annotated) %>% 
                dplyr::mutate(ene.IDs = as.data.frame(pVals.by.pathway)$GeneIDsG) %>% 
                dplyr::rename(KEGGID_pathway = pathway.code) %>% 
                dplyr::mutate(pathway.name =  (pathtable %>% dplyr::filter(crg_code == KEGGID_pathway))$pathname ) %>% 
                # dplyr::mutate(pathway.name = gsub(" -.*","",pathway.name))  %>% 
                dplyr::mutate(p.value =  pVals.by.pathway[,"p.value"])  %>% 
                dplyr::mutate(Day = 'Day2')  %>% 
                dplyr::mutate(modColor = modColor)  %>% 
                dplyr::filter(p.value < 0.05) %>% 
                #na.omit() %>% 
                dplyr::mutate(log10_pvalue = -log10(as.numeric(p.value))) 
    
    KEGG.Day2 <- rbind(KEGG.Day2,outdat) #bind to a cumulative list dataframe
    print(KEGG.Day2) # print to monitor progress
  }
View(KEGG.Day2)
KEGG.Day2OM <- na.omit(KEGG.Day2) %>% 
   dplyr::arrange()
KEGG.Day2OM
library(tidytext)


# View(KEGG.Day2OM)
Day2_KEGGSegmentPlot <- KEGG.Day2 %>%  
                          dplyr::filter(modColor %in% c('brown', 'blue', 'turquoise', 'red', 'pink', 'black')) %>% 
                          dplyr::mutate(modColor = factor(modColor , levels = c('brown', 'blue', 'turquoise', 'red', 'pink', 'black'))) %>%  # for the correct order of facets in the plot below
                          dplyr::mutate(pathway.name = reorder_within(pathway.name, log10_pvalue, modColor)) %>% 
                          ggplot(aes(x=reorder(pathway.name, log10_pvalue), y= log10_pvalue, group = modColor)) + 
                          geom_segment(aes(x=pathway.name, xend=pathway.name, 
                                           y=1, yend=log10_pvalue, color=modColor,
                                       #aes(x=pathway.name, xend=pathway.name, y=min(log10_pvalue), yend=max(log10_pvalue)),  
                                       #linetype=NA, 
                                       size=3)) +   # Draw dashed lines
                          geom_point( aes(col=modColor, size=Num.genes), shape =21,  fill = "white") +   # Draw points
                          scale_color_manual(values = c('brown', 'blue', 'turquoise', 'red', 'pink', 'black')) +
                          #ylim(0,8) +
                          labs(title="Day 2", 
                               x = "Pathway",
                               y = "-Log(pvalue)",
                               subtitle="KEGGREST") +
                          theme_classic() + 
                          theme(
                            panel.grid.minor.y = element_blank(),
                            panel.grid.major.y = element_blank(),
                            legend.position="bottom"
                          ) +
                          xlab("") +
                          ylab("") +
                          ggtitle("Day 2 KEGGREST") +
                          theme(panel.border = element_blank(), # Set border
                                panel.grid.major = element_blank(), #Set major gridlines
                                panel.grid.minor = element_blank(), #Set minor gridlines
                                axis.line = element_line(colour = "black"), #Set axes color
                                plot.background=element_blank()) + #Set the plot background #set title attributes
                          coord_flip() +
                          facet_wrap(modColor ~., 
                                     scales="free_y", 
                                     ncol= 1, 
                                     strip.position="right", 
                                     shrink = T)
#pdf(paste("Analysis/Output/KEGG/WGCNA/Day7_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
print(Day2_KEGGSegmentPlot)



  
  
  
  
  
# day 22 :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Day18_WGCNA_sigmodules <- as.data.frame(c('tan', 
                                            'red', 
                                            'turquoise', 
                                            'salmon', 
                                            'blue',
                                            'green'))


# prep loop for cumulative output table 
df_total             <- data.frame() # start dataframe 
KEGG.Day22           <- data.frame(matrix(nrow = 1, ncol = 9)) # create dataframe to save cumunalitively during for loop
colnames(KEGG.Day22) <- c('Day', 'modColor', 'KEGGID_pathway', 'pathway.name' , 
                         'Num.genes.all', 'Num.genes.exp', 'Gene.IDs', 'p.value', 'log10_pvalue') # names for comuns in the for loop
for (i in 1:nrow(Day18_WGCNA_sigmodules)) {
        modColor <- Day18_WGCNA_sigmodules[i,1]
        # loopmodColor_cor <- paste("MM.", modColor, sep = '')
        # loopmodColor_p   <- paste("p.MM.", modColor, sep = '')
        
        ModuleLoop <- as.data.frame(d18_WGCNA.data %>% dplyr::filter(moduleColor %in% modColor)  %>% 
                                      dplyr::select(c('TranscriptID', 'KEGG_ID', 'geneSymbol', 'MM.p', 'MM.cor')) %>% 
                                      dplyr::filter(.[[4]] < 0.05 & .[[5]] > 0.6) %>% # filter based on thresholds set in the WGCNA + DESeq2 overlap w/PCA 
                                      dplyr::select(c('TranscriptID', 'KEGG_ID', 'MM.p')) %>% 
                                      na.omit() %>% 
                                      dplyr::mutate(KEGG_ID = gsub(".*:","",KEGG_ID)) %>% 
                                      unnest(KEGG_ID))
        geneList <- ModuleLoop[,3]
        names(geneList) <- ModuleLoop$KEGG_ID
        # Wilcoxon test for each pathway
        
        pVals.by.pathway <- t(sapply(names(genes.by.pathway),
                                     function(pathway) {
                                       pathway.genes         <- genes.by.pathway[[pathway]]
                                       list.genes.in.pathway <- intersect(names(geneList), pathway.genes)
                                       list.genes.not.in.pathway <- setdiff(names(geneList), list.genes.in.pathway)
                                       scores.in.pathway <- geneList[list.genes.in.pathway]
                                       scores.not.in.pathway <- geneList[list.genes.not.in.pathway]
                                       if (length(scores.in.pathway) > 0){
                                         p.value <- wilcox.test(scores.in.pathway, 
                                                                (length(pathway.genes) - scores.in.pathway),
                                                                #scores.not.in.pathway,
                                                                alternative = "less")$p.value
                                       } else{
                                         p.value <- NA
                                       }
                                       return(c(p.value = p.value, Num_total_in_pathway = length(pathway.genes), Annotated = length(list.genes.in.pathway), GeneIDs = list(list.genes.in.pathway) ))
                                     }
        ))
        # Assemble output table
        outdat <- data.frame(pathway.code = rownames(pVals.by.pathway)) %>% 
                  dplyr::mutate(Num.genes.all = as.data.frame(pVals.by.pathway)$Num_total_in_pathway) %>%
                  dplyr::mutate(Num.genes.exp = as.data.frame(pVals.by.pathway)$Annotated) %>% 
                  dplyr::mutate(Gene.IDs = as.data.frame(pVals.by.pathway)$GeneIDs) %>% 
                  dplyr::rename(KEGGID_pathway = pathway.code) %>% 
                  dplyr::mutate(pathway.name =  (pathtable %>% dplyr::filter(crg_code == KEGGID_pathway))$pathname ) %>% 
                  # dplyr::mutate(pathway.name = gsub(" -.*","",pathway.name))  %>% 
                  dplyr::mutate(p.value =  pVals.by.pathway[,"p.value"])  %>% 
                  dplyr::mutate(Day = 'Day22')  %>% 
                  dplyr::mutate(modColor = modColor)  %>% 
                  dplyr::filter(p.value < 0.05) %>% 
                  #na.omit() %>% 
                  dplyr::mutate(log10_pvalue = -log10(as.numeric(p.value))) 
        
        KEGG.Day22 <- rbind(KEGG.Day22,outdat) #bind to a cumulative list dataframe
        # print(KEGG.Day22) # print to monitor progress
}
View(KEGG.Day22)

KEGG.Day22OM <-  na.omit(KEGG.Day22) %>% 
  dplyr::mutate(Num.genes = as.numeric(Num.genes))

plot<- KEGG.Day22OM %>%  
  dplyr::filter(modColor %in% c('blue', 
                                'red', 
                                'salmon',
                                #'tan',
                                'green',
                                'turquoise')) %>% 
  dplyr::mutate(modColor = factor(modColor , levels = c('blue', 
                                                        'red', 
                                                        'salmon',
                                                        #'tan',
                                                        'green',
                                                        'turquoise'))) %>%  # for the correct order of facets in the plot below
  dplyr::mutate(pathway.name = reorder_within(pathway.name, log10_pvalue, modColor)) %>% 
  ggplot(aes(x=reorder(pathway.name, log10_pvalue), y= log10_pvalue, group = modColor)) + 
  geom_segment(aes(x=pathway.name, xend=pathway.name, y=1, yend=log10_pvalue, color=modColor,
                   #aes(x=pathway.name, xend=pathway.name, y=min(log10_pvalue), yend=max(log10_pvalue)),  
                   #linetype=NA, 
                   size=3)) +   # Draw dashed lines
  geom_point( aes(col=modColor, size=Num.genes), shape =21,  fill = "white") +   # Draw points
  scale_color_manual(values = c('blue', 
                                'red', 
                                'salmon',
                                #'tan',
                                'green',
                                'turquoise')) +
  #ylim(0,8) +
  labs(title="Day 22", 
       x = "Pathway",
       y = "-Log(pvalue)",
       subtitle="KEGGREST") +
  theme_classic() + 
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="bottom"
  ) +
  xlab("") +
  ylab("") +
  ggtitle("Day 22 KEGGREST") +
  theme(panel.border = element_blank(), # Set border
        panel.grid.major = element_blank(), #Set major gridlines
        panel.grid.minor = element_blank(), #Set minor gridlines
        axis.line = element_line(colour = "black"), #Set axes color
        plot.background=element_blank()) + #Set the plot background #set title attributes
  coord_flip() +
  facet_wrap(modColor ~., 
             scales="free_y", 
             ncol= 1, 
             strip.position="right", 
             shrink = T)
#pdf(paste("Analysis/Output/KEGG/WGCNA/Day7_",modColor,"_RichFactorPlot.pdf", sep =''), width=5, height=6)
print(plot)


  
  
  
  
  
  
  
  
  
  

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

?enrichKEGG
# note: added ModMemership cutoff (Pearson cor and pvalue) 
# when calling 'ModuleLoop_KEGGIDs' - saved to a new directory
for (i in 1:nrow(Day2_WGCNA_sigmodules)) {
  # start with loop by calling the row value common with the 'Master_KEGG_BPTerms' data frind from rbind above 
  modColor <- Day2_WGCNA_sigmodules[i,1]
  
  # call the module color in the Day 7 data
  ModuleLoop <- as.data.frame(d2_WGCNA.data %>% 
                                dplyr::filter(moduleColor %in% modColor)  %>% 
                                dplyr::select(c('TranscriptID', 'KEGG_ID', 'geneSymbol', 'MM.p', 'MM.cor')) %>% 
                                dplyr::filter(.[[4]] < 0.05 & .[[5]] > 0.6) %>% 
                                dplyr::select(c('TranscriptID', 'KEGG_ID', 'MM.p')) %>% 
                                na.omit() %>% 
                                dplyr::mutate(KEGG_ID = gsub(".*:","",KEGG_ID)) %>% 
                                unnest(KEGG_ID))
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
  options(download.file.method = "auto")
  KEGG_vector_cvirg_Cgigas <- as.vector(paste0('LOC', ModuleLoop$KEGG_ID)) # ommit the 'crg:' before the actual terms
  x <- as.vector(c("105348584", "105317413", "105321007", "105318604", "105347464", "105333341", "105318673", "105333480"))
  x2 <- as.vector(as.numeric(ModuleLoop$KEGG_ID))
  data(geneList, package='DOSE')
  de <- names(geneList)[1:100]
  KEGG_cgigas                     <- enrichKEGG(gene = de, 
                                                # organism  = 'crg', # 'hsa' is human 'crg' is pacific oyster 
                                                pvalueCutoff = 0.05) 
  data(geneList, package="DOSE")
  install.packages(c('cli', 'data.table', 'digest', 'dplyr', 'fansi', 'glue', 'Matrix', 'plyr', 'purrr', 'rlang', 'tidyr', 'utf8', 'vctrs'))
  .libPaths()
  install.packages('rlang')
  devtools::install_github("YuLab-SMU/clusterProfiler")
  BiocManager::install("DOSE")
  kk2 <- gseKEGG(geneList     = geneList,
                 organism     = 'hsa',
                 minGSSize    = 120,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
  head(geneList)
        # if loop to output the KEGG enrichment analysis ONLY if genes were successfully mapped...
        if (  nrow(as.data.frame(head(KEGG_cgigas))) > 0 ) {
          # creat dateframe and write the csv file out 
          df                      <- as.data.frame(head(KEGG_cgigas))
          rownames(df)            <- c()
          KEGGoutput              <- as.data.frame(do.call(cbind.data.frame, df))
          KEGGoutput$GeneRatio_2  <- gsub("/"," of ", KEGGoutput$GeneRatio)
          KEGGoutput$Rich_Factor  <- (  (as.numeric(sub("/.*", "", KEGGoutput$GeneRatio))) / (as.numeric(sub("/.*", "", KEGGoutput$BgRatio)))  ) 
          
          write.csv(KEGGoutput, file = paste("Output/WGCNA/day2_larvae/KEGG/RR_cutoff/Day2_",modColor,"_KEGG_allgenes.csv", sep ='')) 
          
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




