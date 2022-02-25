



# Setup: 
# Load libraries

library(dplyr)
library(edgeR)
library(cowplot)
library(ggplot2)
library(reshape2)
library(stringr)

### Set working directory
print(getwd())  # working directory is the scipts folder 
# path for output ting all .csv filtered count files

path = 'C:/Users/samjg/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis/Data/TagSeq/Filtered_Counts/' # personnal computer
# path = 'C:/Users/samuel.gurr/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis/Data/TagSeq/Filtered_Counts' # work computer

raw_counts_mtx <- read.csv("C:/Users/samjg/Documents/Github_repositories/Cvriginica_multistressor/HPC_Analysis/Output/Cvirginica_transcript_count_matrix.csv", sep=',', header=TRUE) # read the output count matrix from sedna
head(raw_counts_mtx) # notice that there are some instances of NAs - prepDE.py3 appeared to skips instances when 
ncol(raw_counts_mtx) -1 # 70 total columns without the transcript ID; 35 total samples (duplicated for sequencing lanes 1 and 2)
View(raw_counts_mtx)

# Edit raw count matrix ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# two major alterations before starting read filtering...
# (1) convert NAs to 0 
# (2) read counts are parsed between sequencing lanes (L001 and L002) - sum within sample ID

# (1) 
raw_counts_mtx[is.na(raw_counts_mtx)] <- 0 # all instances of NA are now 0 

# (2) 
colnames(raw_counts_mtx) # sequencing lanes can be ommitted from the header so samples can be merged by common ID (sum of gene counts) 
raw_counts_melt <- raw_counts_mtx %>%  # convert the matrix from wide to long format 
                          reshape2::melt(id.vars = "transcript_id",
                                        variable.name = "sample_id", 
                                        value.name = "transcript_count")
raw_counts_melt$sample_id <- gsub("_L00.*","",raw_counts_melt$sample_id) # omit the lane identifier from sample_id
sum_raw_counts <- raw_counts_melt %>% # call the long form data
                          mutate(sample_id = as.factor(sub('.[^.]*$', '',raw_counts_melt$sample_id)) ) %>% # grab all of the sample name before the last occurance of a period (i.e. now A1.spat, a2.spat, etc...)
                          mutate(transcript_sample_IDs = paste(transcript_id, sample_id, sep='~')) %>%  # merge the transcript ID and sample Id with unique deliminator '~' 
                          dplyr::group_by(transcript_sample_IDs) %>%  # group by the merged IDs 
                          dplyr::summarise(sum_count = sum(transcript_count)) %>% # sum together the transcript counts by this grouped variable 
                          dplyr::mutate(transcript_id = as.factor(gsub("~.*","",transcript_sample_IDs)) ) %>%  # mutate transcript ID (back from the pasted merge character to group_by)
                          dplyr::mutate(sample_id = as.factor(gsub(".*~","",transcript_sample_IDs)) ) # mutate sample ID (back from the pasted merge character to group_by)
sum_raw_counts <- subset(sum_raw_counts, select = -c(transcript_sample_IDs) ) # ommit the merged column 
# dcast back to the wide format
sum_raw_counts_matrix <- sum_raw_counts %>%
                          dcast(transcript_id  ~ sample_id, value.var = "sum_count")

# Write new raw read file
write.csv(sum_raw_counts_matrix,paste(path,"Sum_raw_count_matrix.csv"))


# Next step (edit here)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Next step (edit here)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Next step (edit here)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Next step (edit here)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Next step (edit here)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Next step (edit here)::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

