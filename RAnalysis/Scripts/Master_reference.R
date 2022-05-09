---
  # title: "KEGG analysis"
  # author: "Samuel Gurr"
  # date: "May 2, 2022"
---
  
# INFORMATION FOR KEGG IN R FOUND HERE: (http://yulab-smu.top/clusterProfiler-book/chapter6.html#kegg-over-representation-test)
  
# LOAD PACKAGES
library(dplyr)
# SET WORKING DIRECTORY   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #


Cvirginica_GOterms  <- unique(read.csv(file="Data/TagSeq/Seq_details/Cviginiva_GOterms.csv", sep=',', header=TRUE) %>%  # Cvirginica GO terms (from collaborator - shared from K. McFarland)
                             dplyr::select(c('GeneID' ,'Annotation_GO_ID', 'Length')) %>% 
                             dplyr::group_by(GeneID) %>% 
                             dplyr::filter(Length == max(Length)))



Cvirginica_annot_reference  <- read.csv(file="Data/TagSeq/Seq_details/seq_id_master.csv", sep=',', header=TRUE) %>% 
                                  dplyr::select(c('TranscriptID','Function','GeneID')) %>% 
                                  dplyr::mutate(TranscriptID = gsub(" ", "", TranscriptID)) %>% # remove the space at the end of each transcript ID
                                  dplyr::mutate(Protein_name = gsub("\\s\\(LOC.*|\\sLOC111.*", "", perl=TRUE, Function)) %>% 
                                  dplyr::select(!Function)
nrow(Cvirginica_annot_reference) #66625 including all genes
nrow(Cvirginica_annot_reference %>% dplyr::filter(grepl('uncharacterized', Protein_name))) #27164 EXCLUDING protein names labeled  as 'uncharacterized'

# nrow(Cvirginica_GOterms) # 34608 - contains gene ID as LOC111... gene length and GO ternms as GO...; GO...; (semicolon sep)

Cgigas_GOterms      <- read.csv(file="Data/TagSeq/Seq_details/Cgigas_refs/Cgigas_GOterms_uniprot.csv", sep=',', header=TRUE) %>%  # Cvirginica GO terms (from uniprot online)
                             dplyr::select(c('Entry.name','Protein.names','Gene.ontology.IDs'))

Cgigas_KEGGIDs      <- as.data.frame(read.delim2(file = "C:/Users/samjg/Documents/Bioinformatics/refs/Cgigas/T03920_(2021_06_23 21_09_09 UTC).nuc", header=FALSE)) %>%  # C gigas KEGG terms from KEGG database (pain subscription allows access to these databases)
                            dplyr::filter(grepl(">crg:", V1)) %>% 
                            dplyr::mutate(Cgigas_KEGGID = gsub("\\s.*|>", "",V1)) %>% 
                            dplyr::mutate(Cgigas_Protein_name = gsub("\\sLOC105.*|*.;", "", substring(V1, 16))) %>% 
                            dplyr::select(!V1)

blastx_CgigasCvirg  <- read.delim2(file="Data/TagSeq/Seq_details/Cgigas_refs/crgKEGG_diamond_out.txt",header=FALSE) %>% # Blast results of Cvirginica to Cgigas using diamond - upload the ouput file (raw!)
                            dplyr::rename(qseqid = V1,sseqid = V2,pident = V3,length = V4,mismatch = V5,gapopen = V6,qstart = V7,qend = V8,sstart = V9,send = V10,evalue = V11,bitscore = V12)  %>% # rename columns 
                            dplyr::group_by(qseqid)
length(unique(blastx_CgigasCvirg$qseqid)) # 61147
# FILTER FOR BEST BLAST HIT (get the KEGIDss merged to Cvrignica) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
blastx_CgigasCvirg_besthits   <- blastx_CgigasCvirg %>% 
                                    dplyr::rename(Cvirginica_TranscriptID = qseqid, Cgigas_KEGGID = sseqid, Cgigas_length = length) %>% # rename these columns for clarity
                                    dplyr::group_by(Cvirginica_TranscriptID) %>% # for each unique transcript ID Cvirginica....
                                    dplyr::mutate(bitscore = as.numeric(bitscore)) %>%  # convert bitscore to umeric - this helps grab the correct value
                                    dplyr::mutate(evalue = as.numeric(evalue)) %>%  # convert bitscore to umeric - this helps grab the correct value
                                    dplyr::filter(bitscore == max(bitscore) | evalue == min(evalue)) %>% # call the  best biscore (highes) and evalue (lowest) from blastx(diamond) against the C gigas genome for KEGG IDs
                                    dplyr::select(c('Cvirginica_TranscriptID','Cgigas_KEGGID', 'Cgigas_length', 'bitscore', 'pident'))  %>%  # call the two target columns - the query (C virginica transcript IDs) and the database of the blast (Cgigas KEGG IDs)
                                    dplyr::mutate(Cgigas_KEGGID= paste0(Cgigas_KEGGID, collapse = ";")) %>%  # concatenate duplicate KEGG IDs to the same row with ; delimiter
                                    dplyr::distinct(Cvirginica_TranscriptID, Cgigas_KEGGID, .keep_all = TRUE) # the line above does not omit duplicate rows - this line does...
length(unique(blastx_CgigasCvirg_besthits$Cvirginica_TranscriptID)) # 61147

# mean sd bitscore
mean(blastx_CgigasCvirg_besthits$bitscore) # 781.6509
sd(blastx_CgigasCvirg_besthits$bitscore) # 927.4943
# mean sd percent identity
mean(as.numeric(blastx_CgigasCvirg_besthits$pident)) # 70.04658
sd(as.numeric(blastx_CgigasCvirg_besthits$pident)) # 18.19463

# check for duplicate rows (sould be zero yielded from theduplyr pipeline above!)
n_occur <- data.frame(table(blastx_CgigasCvirg_besthits$Cvirginica_TranscriptID))
n_occur[n_occur$Freq > 1,]
blastx_CgigasCvirg_besthits[blastx_CgigasCvirg_besthits$Cvirginica_TranscriptID %in% n_occur$Var1[n_occur$Freq > 1],] # should be none!


# check the status of this pipeline...
nrow(blastx_CgigasCvirg_besthits) # 35124 total rows 
nrow(blastx_CgigasCvirg_besthits) == length(unique(blastx_CgigasCvirg_besthits$Cvirginica_TranscriptID)) # must be TRUE - this means there are NO duplicate rows/TranscriptIDs
length(unique(Cvirginica_annot_reference$TranscriptID)) # 66625
length(unique(Cvirginica_GOterms$GeneID)) # 34604

# percent of hits to total Cvirginica transcriptome (mRNA only) 
nrows_Cvirg_mRNA_reference <- nrow(Cvirginica_annot_reference %>% dplyr::filter(!grepl('XR_', TranscriptID))) # 60201
nrows_Cvirg_mRNA_blasthits <- nrow(blastx_CgigasCvirg_besthits %>% dplyr::filter(!grepl('XR_', Cvirginica_TranscriptID ))) # 57662
( nrows_Cvirg_mRNA_blasthits / nrows_Cvirg_mRNA_reference ) * 100 # 95.78246 % of the transcriptome! (excluding ncRNAs!)


# MERGE CVRINGICA GOTERMS WITH THE REFERENCE (ref contains protein name) ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #

# merge the Cvirginica annotated reference with the Cvirnginca GO terms (contains the gene length!) 
Cvirginica_GOterms_MERGE <- merge(Cvirginica_annot_reference, Cvirginica_GOterms, by= 'GeneID') %>% 
                                dplyr::rename(Cvirginica_GeneID = GeneID, Cvirginica_length = Length, Cvirginica_TranscriptID = TranscriptID, Cvirginica_Protein_name = Protein_name)
# nrow(Cvirginica_GOterms_MERGE) # 61891


# MASTER REFERENCE ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::;: #
# we now have 'Cvirginica_GOterms_MERGE' containing the protein name, transcript and gene ID, and GO term annotation - we now want the KEGGIDs from the blastx hits!

Master_ref    <- (as.data.frame(merge(blastx_CgigasCvirg_besthits, Cvirginica_GOterms_MERGE, by ='Cvirginica_TranscriptID')))[,c(1,6,2,7,8,9)] # create a referecne from the Geoduck annotation file to merge with the DE sumamry lists 
nrow(Master_ref) # 59089

ncRNAs_NOT.merged <- blastx_CgigasCvirg_besthits %>%  dplyr::filter(!Cvirginica_TranscriptID %in% Cvirginica_GOterms_MERGE$Cvirginica_TranscriptID)
nrow(ncRNAs_NOT.merged) # 2104 IDs in the Cvirginica_GOterms_MERGE file that do no have the same Transcript ID as the blast hits!


# write master reference
write.csv(Master_ref, file = "Data/TagSeq/Seq_details/Seq_Reference_Master.csv", row.names = FALSE)
