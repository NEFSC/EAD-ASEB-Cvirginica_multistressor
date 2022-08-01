# Purpose: Oyster project 2021 - Respiration rate data 
# analysis of respiration rate data
# Written by: Sam J Gurr (last edit 9/14/2021)

# LOAD PACKAGES :::::::::::::::::::::::::::::::::::::::::::::::::::::::

library(dplyr)
library(ggplot2)
library(reshape2)
library(lme4)
library(lmerTest)
library(performance) # check check_model QC 
library(see)
library(patchwork)
library(forcats)
library(lawstat)
library(car)
library(bestNormalize)

# SET WORKING DIRECTORY :::::::::::::::::::::::::::::::::::::::::::::::

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/") # personal computer
setwd("C:/Users/samuel.gurr/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis") # Work computer

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

exp_metadata <- read.csv(file="Data/ExperimentMetadata.csv", header=T) # treatment assignments to 'Chamber_Tank'
counts_resp  <- read.csv(file="Data/Counts_resp.csv", header=T) # reference for the respirometry data - contains the 'Chamber_Tank' for each sensor channel (whether an animal or a blank)
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

Resp.blanks.merge <- merge(resp.data_merged, blanks_total, by=c("Date", "Notes")) %>% # NOTE: this repeats for every distinct length value
  dplyr::filter(!Chamber_tank =='blank') %>% 
  dplyr::mutate(resp_norm = Lpc - BLANK.mean.Lpc) # calc resp norm - note Lpc is still raw data  - thus a positive resp norm means blank > sample resp - bad data! 

Resp.blanks.merge_OM <- Resp.blanks.merge %>% dplyr::filter(!resp_norm > 0)  # omits respiration rate values showing an increase in O2 over time 
Resp.outliers <- Resp.blanks.merge_OM %>% dplyr::filter(resp_norm > 0)  # lets lopok at this data - possible troubleshootiung to ID these outliers...

# calculate resp rates
vial.vol <- 0.08 # milliliters (ml)
Resp.Master <- merge(Resp.blanks.merge_OM[,c(1,3,5,9:13,17)], counts_resp[,c(1,6:7)], by = c('Date','Chamber_tank'))
Resp.Master$resp_ng_L_indiv_hr <- ( 
  ( ( (abs(Resp.Master$resp_norm)*1000000) * # call absolute value of resp in mg per minute - convert to ng min-1
      (vial.vol/1000) ) / # correct ng minute-1 to ng liter-1 by multiplying by the resp vial in liters
      Resp.Master$Counts ) * # normalize by individual or larvae count - as to ng L-1 individual-1
    (60)) # correct for the time; final value is ng Liter-1 individual-1 hour-1

# mean sd rates
mean(Resp.Master$resp_ng_L_indiv_hr) # mean = 0.4994656
sd(Resp.Master$resp_ng_L_indiv_hr) # sd= 0.8197345





# write master Resp file ::::::::::::::::::::::::::::::::::;
write.csv(Resp.Master, 
          "Output/Respiration/RespirationMaster.csv")# write







# ANALYSIS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# model effect of treatment on resp rate 20210507
Resp_APRIL <- Resp.Master %>% 
  dplyr::filter(Date %in% '4/30/2021') %>% 
  dplyr::filter(!resp_ng_L_indiv_hr >1)


# Stats
# Factorial ANOVA: has two or more categorical independent variables 
# (either with or without the interactions) and a single normally distributed interval dependent variable
LMmod.APRIL   <- aov(lm(resp_ng_L_indiv_hr~Temp*pCO2*Salinity,data=Resp_APRIL))
summary(LMmod.APRIL)
check_model(LMmod.APRIL) # observe the diagnostics of the model
shapiro.test(residuals(LMmod.APRIL)) # non normal
MASS::truehist(Resp_APRIL$resp_ng_L_indiv_hr, nbins = 12) # left skew
leveneTest(LMmod.APRIL) # good

# transformation, check for normality of historgram 
library('bestNormalize')
(boxcox_April      <- boxcox(Resp_APRIL$resp_ng_L_indiv_hr))
(archsign_April    <- arcsinh_x(Resp_APRIL$resp_ng_L_indiv_hr))
(yeojohnson_April  <- yeojohnson(Resp_APRIL$resp_ng_L_indiv_hr))
(orderNorm_April   <- orderNorm(Resp_APRIL$resp_ng_L_indiv_hr))
par(mfrow = c(2,2))
MASS::truehist(boxcox_April$x.t, main = "YBoxcox transformation", nbins = 12) #looks good
MASS::truehist(archsign_April$x.t, main = "Archsign transformation", nbins = 12)
MASS::truehist(yeojohnson_April$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(orderNorm_April$x.t, main = "orderNorm transformation", nbins = 12)

Resp_APRIL$resp_ng_L_indiv_hr_T <- boxcox_April$x.t
LMmod.APRILboxcox   <- aov(lm(resp_ng_L_indiv_hr_T~Temp*pCO2*Salinity,data=Resp_APRIL))
summary(LMmod.APRILboxcox)
check_model(LMmod.APRILboxcox) # observe the diagnostics of the model
shapiro.test(residuals(LMmod.APRILboxcox)) # 0.09613 normal
leveneTest(LMmod.APRILboxcox) # 0.7158 good




MEmod.APRIL    <- lmer(resp_ng_L_indiv_hr_T~Temp*pCO2*Salinity + (1|Chamber_tank),REML=TRUE, data=Resp_APRIL)
summary(MEmod.APRIL)
check_model(MEmod.APRIL)
leveneTest(MEmod.APRIL) # good

# Post-hoc tests and exploration of sig effects
TukeyHSD(AOV_APRIL)

# Figures
Resp_APRIL_select  <- Resp_APRIL %>% dplyr::select(c('resp_ng_L_indiv_hr', 'Temp', 'pCO2', 'Salinity'))
Resp_APRIL_melt    <- tidyr::gather(Resp_APRIL_select, variable, value, -resp_ng_L_indiv_hr)

ggplot(Resp_APRIL_melt, aes(value , resp_ng_L_indiv_hr , fill = factor(value ))) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y")

Plot <- Resp_APRIL_select %>%
  dplyr::mutate(full.treatment = (paste(Temp, pCO2, Salinity,sep=''))) %>%
  dplyr::mutate(full.treatment = fct_relevel(full.treatment,
                            "HHH", "HLH", "LHH",'LLH',
                            "HHL", "HLL", "LHL", 'LLL')) %>%
 
  ggplot(aes(full.treatment, resp_ng_L_indiv_hr , fill = factor(full.treatment))) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  scale_fill_manual(values=c("#D55E00","#56B4E9", "#D55E00","#56B4E9", 
                             "#D55E00","#56B4E9", "#D55E00","#56B4E9")) +
  labs(title = "C virginica, April 30th larvae", 
       y = expression(Respiration~rate~"("~ng~L^{-1}~O[2]%.%indiv^{-1}%.% hr^{-1}~")"), 
       x = "Full treatment (Temp/pCO2/Sal)") + 
  annotate("text", x=6.5, y=.75, label = "Low Salinity") +
  annotate("rect", xmin = 4.5, xmax = 9, ymin = 0, ymax = 1,
           alpha = .2) +
  theme_bw() 
Plot









# model effect of treatment on resp rate 20210507
Resp_MAY <- Resp.Master %>% 
  dplyr::filter(Date %in% '5/7/2021') %>% 
  dplyr::filter(!resp_ng_L_indiv_hr >1)


# Stats
# Factorial ANOVA: has two or more categorical independent variables 
# (either with or without the interactions) and a single normally distributed interval dependent variable
LMmod.MAY   <- aov(lm(resp_ng_L_indiv_hr~Temp*pCO2*Salinity,data=Resp_MAY))
summary(LMmod.MAY)
check_model(LMmod.MAY)
shapiro.test(residuals(LMmod.MAY)) # 0.0008179
leveneTest(LMmod.MAY) # 0.2802 good

# transformation, check for normality of historgram using 'bestNormalize
(boxcox_May      <- boxcox(Resp_MAY$resp_ng_L_indiv_hr))
(archsign_May    <- arcsinh_x(Resp_MAY$resp_ng_L_indiv_hr))
(yeojohnson_May  <- yeojohnson(Resp_MAY$resp_ng_L_indiv_hr))
(orderNorm_May   <- orderNorm(Resp_MAY$resp_ng_L_indiv_hr))
par(mfrow = c(2,2))
MASS::truehist(boxcox_May$x.t, main = "YBoxcox transformation", nbins = 12) #looks good
MASS::truehist(archsign_May$x.t, main = "Archsign transformation", nbins = 12)
MASS::truehist(yeojohnson_May$x.t, main = "Yeo-Johnson transformation", nbins = 12)
MASS::truehist(orderNorm_May$x.t, main = "orderNorm transformation", nbins = 12)

Resp_MAY$resp_ng_L_indiv_hr_T <- boxcox_May$x.t
LMmod.Mayboxcox   <- aov(lm(resp_ng_L_indiv_hr_T~Temp*pCO2*Salinity,data=Resp_MAY))
summary(LMmod.Mayboxcox)
check_model(LMmod.Mayboxcox) # observe the diagnostics of the model
shapiro.test(residuals(LMmod.Mayboxcox)) # 0.5787 normal
leveneTest(LMmod.Mayboxcox) # 0.7022 good
plot(LMmod.Mayboxcox)



# Figures
Resp_MAY_select  <- Resp_MAY %>% dplyr::select(c('resp_ng_L_indiv_hr', 'Temp', 'pCO2', 'Salinity'))
Resp_MAY_melt    <- tidyr::gather(Resp_MAY_select, variable, value, -resp_ng_L_indiv_hr)

ggplot(Resp_MAY_melt, aes(value , resp_ng_L_indiv_hr , fill = factor(value ))) +
  theme(panel.grid=element_blank()) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_y")

Plot <- Resp_MAY_select %>%
  
  dplyr::mutate(full.treatment = (paste(Temp, pCO2, Salinity,sep=''))) %>%
  # dplyr::mutate(full.treatment = (paste(Temp, Salinity,sep=''))) %>%
  # dplyr::mutate(full.treatment = (paste(Temp, pCO2,sep=''))) %>%
  
  # dplyr::mutate(full.treatment = fct_relevel(full.treatment,
  #                           "HHH", "HLH", "LHH",'LLH',
  #                           "HHL", "HLL", "LHL", 'LLL')) %>%
  
  ggplot(aes(full.treatment, resp_ng_L_indiv_hr , fill = factor(full.treatment))) +
  scale_color_manual(values=c("#56B4E9","#D55E00")) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  theme_bw()
Plot
