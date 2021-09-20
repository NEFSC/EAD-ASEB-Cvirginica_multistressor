# Purpose: Oyster project 2021 - Respiration rate data 
# analysis of respiration rate data

# Written by: Sam J Gurr (last edit 9/14/2021)

# LOAD PACKAGES :::::::::::::::::::::::::::::::::::::::::::::::::::::::

library(dplyr)
library(ggplot2)
library(forcats)

# SET WORKING DIRECTORY :::::::::::::::::::::::::::::::::::::::::::::::

setwd("C:/Users/samjg/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis")

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

exp_metadata <- read.csv(file="Data/ExperimentMetadata.csv", header=T) # treatment assignments to 'Chamber_Tank'
resp.ref     <- read.csv(file="Data/Respiration/Reference_master.csv", header=T) # reference for the respirometry data - contains the 'Chamber_Tank' for each sensor channel (whether an animal or a blank)
resp.data    <- read.csv(file="Output/Respiration/Cumulative_resp_alpha0.4.csv", header=T) # read the calculate raw rates from 'resp_LoLin' script - contains the calculated rate (not normalized for blanks) for each sensor-channel

# merge the exp_metadata with the resp.data
resp.ref_merged                 <- merge(exp_metadata, resp.ref, by = 'Chamber_tank', all=TRUE) # all TRUE allows us to keep the blanks
resp.data_merged                <- merge(resp.data, resp.ref_merged, by = c('Date', 'Channel','Notes','Filename')) # out master file moving forward....
resp.data_merged$TempCarbSal    <- paste(resp.data_merged$Temp, resp.data_merged$pCO2, resp.data_merged$Salinity, sep ='')

# CALCULATE RESPIRATION RATES :::::::::::::::::::::::::::::::::::::::::::::::

# ---------------------- (1) --------------------------------#
# get a summary table of blanks to normalize respiration rates

dates.runs <- resp.data_merged %>%  # call table
  dplyr::distinct(Date, Notes) # call all unique values for date run and sw condition
dates.runs <- na.omit(dates.runs)
# call dataframe and build table to rbind in for loop
blanks_total <- data.frame() # start dataframe 
blanks.table <- data.frame(matrix(nrow = 1,ncol = 5)) # make a table template
colnames(blanks.table)<-c('Date', 'Notes', 'BLANK.mean.Lpc', 'BLANK.mean.Leq' , 'BLANK.mean.Lz') # names for comuns in the for loop

# for loop. objective = obtian a mean value for all blanks specific to date, run #, seawater treatment
for(i in 1:nrow(dates.runs)) {
  data <- resp.data_merged %>% 
    dplyr::select(Date, Chamber_tank, Notes, Lpc,  Leq, Lz) %>% 
    dplyr::filter(!is.na(Lpc)) %>% # ommits empty resp channels (if any)
    dplyr::filter(Notes == dates.runs[i,2])
  
  blanks <- data %>%
    dplyr::filter(Chamber_tank == "blank") %>% 
    dplyr::summarise(mean_Lpc = mean(abs(Lpc)),
                     mean_Leq = mean(abs(Leq)), 
                     mean_Lz = mean(abs(Lz)))
  
  blanks.table$Date           <- dates.runs[i,1] # all files have date in the form of yyyymmdd at the start of each csv name
  #blanks.table$RUN <- dates.runs[i,2] # assign the run to the number in the title for the trials completed that day
  blanks.table$Notes      <- dates.runs[i,2]
  blanks.table$BLANK.mean.Lpc <- blanks[1,1]
  blanks.table$BLANK.mean.Leq <- blanks[1,2]
  blanks.table$BLANK.mean.Lz  <- blanks[1,3]
  # blanks.table$alpha <- data[1,9] # set at start of script - reresents the proportion of data for final estimate of slopes (Lpc, Leq, Lz)
  
  df <- data.frame(blanks.table) # name dataframe for this singl e row
  blanks_total <- rbind(blanks_total,df) #bind to a cumulative list dataframe
  print(blanks_total) # print to monitor progress
}

blanks_total # view blanks table


# ---------------------- (2) --------------------------------#
# merge blanks with the summary table and calculate the normalized rates 

Resp.Master <- merge(resp.data_merged, blanks_total, by=c("Date", "Notes")) %>% # NOTE: this repeats for every distinct length value
  dplyr::filter(!Chamber_tank =='blank') %>% 
  dplyr::mutate(resp_norm = Lpc - BLANK.mean.Lpc) # ommits respiration rate values showing an increase in O2 over time 

Resp.Master_OM <- Resp.Master %>% dplyr::filter(!resp_norm > 0) 
# NOTE: look at the following table to troubleshoot if needed
Resp.outliers <- Resp.Master %>% dplyr::filter(resp_norm > 0) 

# calculate resp rates
vial.vol <- 0.08 # milliliters (ml)
Resp.Master_OM$resp_ug_L_hr <- ( (abs(Resp.Master_OM$resp_norm)*1000)*(vial.vol/1000)*(60)) # correct for Liters and hours (currently in minutes fro the resp.lolin.R script)

# mean sd rates
mean(Resp.Master_OM$resp_ug_L_hr)
sd(Resp.Master_OM$resp_ug_L_hr)


# ANALYSIS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# model effect of treatment on resp rate 20210507
Resp_0430 <- Resp.Master_OM %>% dplyr::filter(Date %in% '4/30/2021') %>% dplyr::filter(!resp_ug_L_hr >1)

mod <- aov(lm(resp_ug_L_hr~Temp*pCO2*Salinity,data=Resp_0430))
summary(mod)

Resp_0430_FIGS <- Resp_0430 %>% dplyr::select(c('resp_ug_L_hr', 'Temp', 'pCO2', 'Salinity'))
Resp_0430_melt <- melt(Resp_0430_FIGS, id=c("resp_ug_L_hr"))

ggplot(Resp_0430_melt, aes(value , resp_ug_L_hr , fill = factor(value ))) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y")

Plot <- Resp_0430_FIGS %>%
  dplyr::mutate(full.treatment = (paste(Temp, pCO2, Salinity,sep=''))) %>%
  dplyr::mutate(full.treatment = fct_relevel(full.treatment,
                            "HHH", "HLH", "LHH",'LLH',
                            "HHL", "HLL", "LHL", 'LLL')) %>%
  ggplot(aes(full.treatment, resp_ug_L_hr , fill = factor(full.treatment))) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  theme_bw()


# model effect of treatment on resp rate 20210507
Resp_0507 <- Resp.Master_OM %>% dplyr::filter(Date %in% '5/7/2021') %>% dplyr::filter(!resp_ug_L_hr >1)

mod <- aov(lm(resp_ug_L_hr~Temp*pCO2*Salinity,data=Resp_0507))
summary(mod)

Resp_0507_FIGS <- Resp_0507 %>% dplyr::select(c('resp_ug_L_hr', 'Temp', 'pCO2', 'Salinity'))
Resp_0507_melt <- melt(Resp_0507_FIGS, id=c("resp_ug_L_hr"))

ggplot(Resp_0507_melt, aes(value , resp_ug_L_hr , fill = factor(value ))) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y")

# Plot <- Resp_0507_FIGS %>% 
#   dplyr::mutate(full.treatment = (paste(Temp, pCO2, Salinity,sep=''))) %>% 
#   dplyr::mutate(full.treatment = fct_relevel(full.treatment, 
#                             "HHH", "HLH", "LHH",'LLH', 
#                             "HHL", "LHL",'LLL')) %>% 
#   ggplot(aes(full.treatment, resp_ug_L_hr , fill = factor(full.treatment))) +
#   scale_color_manual(values=c("#56B4E9","#D55E00")) +
#   geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
#   geom_boxplot(size=0.2, alpha=0.1) +
#   theme_bw() 

