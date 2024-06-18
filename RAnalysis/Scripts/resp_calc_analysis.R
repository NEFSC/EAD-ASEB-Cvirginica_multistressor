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
library(car)


# SET WORKING DIRECTORY :::::::::::::::::::::::::::::::::::::::::::::::

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/") # personal computer

# setwd("C:/Users/samuel.gurr/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis") # Work computer

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

exp_metadata <- read.csv(file="Data/ExperimentMetadata.csv", header=T) # treatment assignments to 'Chamber_Tank'
counts_resp  <- read.csv(file="Data/Counts_resp.csv", header=T) # reference for the respirometry data - contains the 'Chamber_Tank' for each sensor channel (whether an animal or a blank)
resp.ref     <- read.csv(file="Data/Respiration/Reference_master.csv", header=T) # reference for the respirometry data - contains the 'Chamber_Tank' for each sensor channel (whether an animal or a blank)
resp.data    <- read.csv(file="Output/Respiration/Cumulative_resp_alpha0.4.csv", header=T) %>% # read the calculate raw rates from 'resp_LoLin' script - contains the calculated rate (not normalized for blanks) for each sensor-channel
                          dplyr::filter(!Filename %in% '1_3_19_21_raw.txt') # use 1_3_19_21_new_sensor_for_7_raw.txt
resp.data[2,3] <- -0.0046 # 1_3_19_21_raw - CH2
resp.data[3,3] <- -0.0072 # 1_3_19_21_raw - CH3

# merge the exp_metadata with the resp.data
resp.ref_merged                 <- merge(exp_metadata, resp.ref, by = 'Chamber_tank', all=TRUE) # all TRUE allows us to keep the blanks
resp.data_merged                <- merge(resp.data, resp.ref_merged, by = c('Date', 'Channel','Notes','Filename')) # out master file moving forward....
resp.data_merged$TempCarbSal    <- paste(resp.data_merged$Temp, resp.data_merged$pCO2, resp.data_merged$Salinity, sep ='')




# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# CALCULATE RESPIRATION RATES ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



# ------------------------------------------- (1) -------------------------------------------------------#
# get a summary table of blanks to normalize respiration rates

dates.runs <- resp.data_merged %>%  # call table
  dplyr::distinct(Date, Notes) # call all unique values for date run and sw condition
dates.runs <- na.omit(dates.runs)


# call dataframe and build table to rbind in for loop
blanks_total <- data.frame() # start dataframe 
# for loop. objective = obtian a mean value for all blanks specific to date, run #, seawater treatment
for(i in 1:nrow(dates.runs)) {
  data <- resp.data_merged %>% 
    dplyr::select(Date,Channel, Chamber_tank, Notes, Lpc,  Leq, Lz) %>% 
    dplyr::filter(!is.na(Lpc)) %>% # ommits empty resp channels (if any)
    dplyr::filter(Notes == dates.runs[i,2])
  
  blanks <- data %>%
    dplyr::filter(Chamber_tank == "blank") %>% 
    dplyr::select(!Chamber_tank) %>% 
    dplyr::mutate(filename = dates.runs[i,2])
  
  blanks.table <- data.frame(matrix(nrow = nrow(blanks),ncol = 7)) # make a table template
  colnames(blanks.table)<-c('Date', 'Channel', 'BLANK.Lpc', 'BLANK.Leq' , 'BLANK.Lz', 'filename', 'Notes') # names for comuns in the for loop
  
  blanks.table$Date      <- blanks$Date
  blanks.table$Channel   <- blanks$Channel
  blanks.table$Notes     <- blanks$Notes
  blanks.table$BLANK.Lpc <- blanks$Lpc
  blanks.table$BLANK.Leq <- blanks$Leq
  blanks.table$BLANK.Lz  <- blanks$Lz
  blanks.table$filename  <- blanks$filename
  # blanks.table$alpha <- data[1,9] # set at start of script - reresents the proportion of data for final estimate of slopes (Lpc, Leq, Lz)
  
  df <- data.frame(blanks.table) # name dataframe for this singl e row
  blanks_total <- rbind(blanks_total,df) #bind to a cumulative list dataframe
  print(blanks_total) # print to monitor progress
}

blanks_total2 <- blanks_total %>% 
  #dplyr::filter(!(filename =='20210319_new_sensor_7' & Channel == 'CH8')) %>% # positive blank rate - omit, use other blank in file
  dplyr::filter(!(filename =='20210319' & Channel == 'CH4')) %>% # slightly positive rate - omit, use other blank in the file 
  dplyr::filter(!(filename =='20210430_LOWtemp_HIGHsal' & Channel == 'CH8')) %>% # abnormally higher than other blank, error  -omit, use other blank on file
  dplyr::filter(!(filename =='20210507_LOWtemp_HIGHsal' & Channel == 'CH4')) %>% # abnormally higher than other blank, error -omit, use other blank on file
  mutate(across(everything(), ~ ifelse(. < 0, 0, .))) # all other positive rates hug zero, (i.e. 0.00089) - thus make zero 
  
blanks_means <- blanks_total2 %>% 
  group_by(Date, Notes) %>% 
  dplyr::summarise(BLANK.mean.Lpc = mean(abs(BLANK.Lpc)),
                   BLANK.mean.Leq = mean(abs(BLANK.Leq)), 
                   BLANK.mean.Lz = mean(abs(BLANK.Lz)))




# -------------------------------------------------- (2) -----------------------------------------------------------------#
# merge blanks with the summary table and calculate the normalized rates 

Resp.blanks.merge <- merge(resp.data_merged, blanks_means, by=c("Date", "Notes")) %>% # NOTE: this repeats for every distinct length value
  dplyr::filter(!Chamber_tank =='blank') %>% 
  dplyr::filter(!Lpc > 0) %>% # 14 with positive rate raw value - omitted
  dplyr::mutate(resp_norm = abs(Lpc) - BLANK.mean.Lpc) # calc resp norm - note Lpc is still raw data  - thus a positive resp norm means blank > sample resp - bad data! 



# filter for postive rates - look at the outliers in which the animal rates were < the blank! 
Resp.blanks.merge_OM <- Resp.blanks.merge %>% dplyr::filter(!resp_norm < 0)  # omits respiration rate values showing an increase in O2 over time 
Resp.outliers <- Resp.blanks.merge %>% dplyr::filter(resp_norm < 0)  # only 3 rates were ommitted

# calculate resp rates
# vial.vol <- 0.08 # milliliters (ml)
vial.vol <- 2.2 # milliliters (ml) - small Loligo chambers
Resp.Master <- merge(Resp.blanks.merge_OM[,c(1,3,5,9:13,17)], counts_resp[,c(1,6:7)], by = c('Date','Chamber_tank'))
Resp.Master$resp_ng_L_indiv_hr <- ( 
  ( ( (abs(Resp.Master$resp_norm)*1000000) * # call absolute value of resp in mg per minute - convert to ng min-1 (note negatives should aredy be omitted!_)
      (vial.vol/1000) ) / # correct ng minute-1 to ng liter-1 by multiplying by the resp vial in liters
      Resp.Master$Counts ) * # normalize by individual or larvae count - as to ng L-1 individual-1
    (60)) # correct for the time; final value is ng Liter-1 individual-1 hour-1

# mean sd rates
mean(Resp.Master$resp_ng_L_indiv_hr) # mean = 1.423157
sd(Resp.Master$resp_ng_L_indiv_hr) # sd= 2.027203





# write master Resp file ------------------------------------------------------------------------------------- #
write.csv(Resp.Master, 
          "Output/Respiration/RespirationMaster.csv")# write




# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ANALYSIS :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

library(emmeans)
library(car)
library(ggpubr)
library(rcompanion)
# APRIL 24 hours post ferilization ------------------------------------------------ #

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/") # personal computer
Resp.Master  <- read.csv(file="Output/Respiration/RespirationMaster.csv", header=T) # reference for the respirometry data - contains the 'Chamber_Tank' for each sensor channel (whether an animal or a blank)

# model effect of treatment on resp rate 20210430
Resp_APRIL <- Resp.Master %>% 
  dplyr::filter(Date %in% '4/30/2021') # %>% 

# Three way ANOVA 
LMmod.APRIL   <- aov(lm(resp_ng_L_indiv_hr~Temp*pCO2*Salinity,data=Resp_APRIL))
shapiro.test(residuals(LMmod.APRIL)) # 0.009733 - non normal
# LOG transform
LMmod.APRIL_T   <- aov(lm(log(resp_ng_L_indiv_hr)~Temp*pCO2*Salinity,data=Resp_APRIL))
shapiro.test(residuals(LMmod.APRIL_T)) # 0.6472
leveneTest(LMmod.APRIL_T) # 0.6339
summary(LMmod.APRIL_T)
# Df Sum Sq Mean Sq F value  Pr(>F)   
# Temp                1  4.622   4.622   5.769 0.02881 * 
# pCO2                1  0.303   0.303   0.378 0.54710   
# Salinity            1  1.523   1.523   1.901 0.18693   
# Temp:pCO2           1  8.930   8.930  11.145 0.00417 **
# Temp:Salinity       1  1.591   1.591   1.986 0.17792   
# pCO2:Salinity       1  2.690   2.690   3.357 0.08560 . 
# Temp:pCO2:Salinity  1  0.251   0.251   0.314 0.58321   
# Residuals          16 12.821   0.801 
# post hoc tests 
library(emmeans)
posthoc<-emmeans(LMmod.APRIL_T, pairwise~Temp:pCO2, adjust="tukey")
multcomp::cld(posthoc$emmeans,alpha = 0.5, Letters = letters)
# Temp pCO2 emmean    SE df lower.CL upper.CL .group
# H    L    -1.819 0.365 16   -2.593  -1.0439  a    
# L    H    -0.716 0.365 16   -1.491   0.0586   b   
# H    H    -0.374 0.365 16   -1.149   0.4009   bc  
# L    L     0.279 0.365 16   -0.496   1.0538    c  


# Figures
Resp_APRIL_select  <- Resp_APRIL %>% 
  dplyr::select(c('resp_ng_L_indiv_hr', 'Temp', 'pCO2', 'Salinity')) %>% 
  dplyr::mutate(Age = '24hrs')
# Resp_APRIL_melt <- tidyr::gather(Resp_APRIL_select, variable, value, -resp_ng_L_indiv_hr)

APRIL_all <- Resp_APRIL_select %>%
  dplyr::mutate(Temp_pCO2 = paste(Temp,pCO2, sep = '_')) %>% 
  dplyr::mutate(full.treatment = (paste(Salinity, pCO2, Temp,sep=''))) %>%
  dplyr::mutate(full.treatment = fct_relevel(full.treatment,
                            "LHL", "LHH", "LLL",'LLH',
                            "HHL", "HHH", "HLL", 'HLH')) %>%
  ggplot(aes(Temp_pCO2, resp_ng_L_indiv_hr , fill = factor(Salinity))) +
  geom_boxplot(size=0.2, alpha=0.1) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  scale_fill_manual(values=c("white", "grey40")) +
  labs(title = "C virginica, 24 hr larvae", 
       y = expression(Respiration~rate~"("~ng~L^{-1}~O[2]%.%indiv^{-1}%.% hr^{-1}~")"), 
       x = "Temp_pCO2") + 
  # annotate("text", x=2, y=5.8, label = "Low Salinity") +
  #annotate("rect", xmin = 0, xmax = 4.5, ymin = 0, ymax = 6.5,alpha = .2) +
  theme_classic() 
APRIL_all


APRIL_TemppCO2 <- Resp_APRIL_select %>%
  dplyr::mutate(Temp_pCO2 = paste(Temp,pCO2, sep = '_')) %>% 
  ggplot(aes(Temp_pCO2, resp_ng_L_indiv_hr , fill = factor(Temp_pCO2))) +
  geom_boxplot(size=0.2, alpha=0.1) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  # scale_fill_manual(values=c("white", "grey40")) +
  labs(title = "C virginica, 24 hr larvae", 
       y = expression(Respiration~rate~"("~ng~L^{-1}~O[2]%.%indiv^{-1}%.% hr^{-1}~")"), 
       x = "Temp_pCO2") + 
  annotate("text", x=2.2, y=0.8, label = "a", size  =5) +
  annotate("text", x=1.2, y=1.2, label = "ab", size =5) +
  annotate("text", x=3.2, y=1.2, label = "b", size  =5) +
  annotate("text", x=4.2, y=3.2, label = "c", size  =5) +
  #annotate("rect", xmin = 0, xmax = 4.5, ymin = 0, ymax = 6.5,alpha = .2) +
  theme_classic() 
APRIL_TemppCO2


setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/")
pdf("Respiration/Day1_RR.pdf", width=8, height=12)
ggarrange(APRIL_all,APRIL_TemppCO2, nrow=2)
dev.off()


Resp_means_APRIL <- Resp_APRIL_select %>% 
  na.omit() %>% 
  dplyr::group_by(Temp, pCO2, Salinity, Age) %>% 
  dplyr::summarise(mean_RR = mean(resp_ng_L_indiv_hr), 
                   n       = n(),
                   sd_RR   = sd(resp_ng_L_indiv_hr),
                   se_RR   = sd_RR/(sqrt(n)))


APRIL_MeanSE <- Resp_means_APRIL %>%
  dplyr::mutate(Temp_pCO2 = paste(Temp,pCO2, sep = '_')) %>% 
  ggplot(aes(Temp_pCO2, mean_RR , fill = factor(Salinity))) +
  geom_errorbar(aes(ymin = mean_RR - se_RR, 
                    ymax = mean_RR + se_RR), 
                width = 0.5, 
                position= "dodge2") +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.1))+
  scale_fill_manual(values=c("white", "grey40")) +
  labs(title = "C virginica, 24 hr larvae", 
       y = expression(Respiration~rate~"("~ng~L^{-1}~O[2]%.%indiv^{-1}%.% hr^{-1}~")"), 
       x = "Temp_pCO2") + 
  theme_classic() 
APRIL_MeanSE

Heatplot_RR_APRIL <- Resp_means_APRIL %>% 
                        dplyr::mutate(OA_Sal = (paste(pCO2, Salinity,sep=''))) %>% 
                        ggplot(aes(x = as.factor(Age),
                                   y = mean_RR)) +
                                   #shape = stage)) + 
                        geom_rect(aes(fill = mean_RR), 
                                  xmin = -Inf, 
                                  xmax = Inf, 
                                  ymin = -Inf, 
                                  ymax = Inf, 
                                  alpha = 0.3) +
                        geom_point(color = 'black') +
                        geom_errorbar(aes(ymin = mean_RR - se_RR, 
                                          ymax = mean_RR + se_RR), 
                                      width = 0.5, 
                                      position= "dodge2") +
                        theme_bw() +  
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank()) +
                        facet_grid(vars(Temp), 
                                   vars(fct_relevel(OA_Sal, c("HL", "LL", "HH", "LH")))) +
                        scale_fill_gradient(low = "orange", 
                                            high = "forestgreen", 
                                            aesthetics = "fill")


setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/")
pdf("Respiration/Day1_RR_heatplot.pdf", width=8, height=6)
print(Heatplot_RR_APRIL)
dev.off()


setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/")
pdf("Respiration/Day1_RR_posthoc.pdf", width=8, height=6)
print(APRIL_MeanSE)
dev.off()
# MAY - 8 days post fertilization---------------------------------------------------------- #





# model effect of treatment on resp rate 20210507
Resp_MAY <- Resp.Master %>% 
  dplyr::filter(Date %in% '5/7/2021')  %>% 
  dplyr::mutate(Age = '8days') %>% 
  # dplyr::filter(Temp %in% 'L') %>%  # call only the high temp
  dplyr::mutate()


# Three way ANOVA 
LMmod.MAY   <- aov(lm(resp_ng_L_indiv_hr~Temp*pCO2*Salinity,data=Resp_MAY))
shapiro.test(residuals(LMmod.MAY)) # 0.006028
# log transformation
LMmod.MAY_T   <- aov(lm(log(resp_ng_L_indiv_hr)~Temp*pCO2*Salinity,data=Resp_MAY))
shapiro.test(residuals(LMmod.MAY_T)) # 0.7608
leveneTest(LMmod.MAY_T) # 0.8419 good
summary(LMmod.MAY_T)
# Df Sum Sq Mean Sq F value Pr(>F)  
# Temp           1  3.355   3.355   3.544 0.0924 .
# pCO2           1  0.972   0.972   1.027 0.3373  
# Salinity       1  3.664   3.664   3.870 0.0807 .
# Temp:pCO2      1  3.794   3.794   4.008 0.0763 .
# Temp:Salinity  1  0.022   0.022   0.023 0.8825  
# pCO2:Salinity  1  0.117   0.117   0.123 0.7334  
# Residuals      9  8.521   0.947

# Figures
Resp_MAY_select  <- Resp_MAY %>% 
  dplyr::select(c('resp_ng_L_indiv_hr', 'Temp', 'pCO2', 'Salinity')) %>% 
  dplyr::mutate(Age = '8days')


MAY_all <- Resp_MAY_select %>%
  dplyr::mutate(full.treatment = (paste(Salinity, pCO2, Temp,sep=''))) %>%
  dplyr::mutate(full.treatment = fct_relevel(full.treatment,
                                             "LHL", "LHH", "LLL",#'LLH',
                                             "HHL", "HHH", "HLL", 'HLH')) %>%
  
  ggplot(aes(full.treatment, resp_ng_L_indiv_hr , fill = factor(full.treatment))) +
  geom_point(shape = 21, size = 2, position = position_jitterdodge(jitter.width = 0.5))+
  geom_boxplot(size=0.2, alpha=0.1) +
  scale_fill_manual(values=c("#56B4E9", "#D55E00","#56B4E9", #"#D55E00",
                             "#56B4E9", "#D55E00","#56B4E9", "#D55E00")) +
  labs(title = "C virginica, 8 day larvae", 
       y = expression(Respiration~rate~"("~ng~L^{-1}~O[2]%.%indiv^{-1}%.% hr^{-1}~")"), 
       x = "Full treatment (Sal/pCO2/Temp)") + 
  annotate("text", x=2, y=8, label = "Low Salinity") +
  annotate("rect", xmin = 0, xmax = 3.5, ymin = 0, ymax = 10,
           alpha = .2) +
  theme_bw() 
MAY_all

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/")
pdf("Respiration/Day8_RR.pdf", width=8, height=6)
print(MAY_all)
dev.off()


Resp_means_MAY <- Resp_MAY_select %>% 
  na.omit() %>% 
  dplyr::group_by(Temp, pCO2, Salinity, Age) %>% 
  dplyr::summarise(mean_RR = mean(resp_ng_L_indiv_hr), 
                   n       = n(),
                   sd_RR   = sd(resp_ng_L_indiv_hr),
                   se_RR   = sd_RR/(sqrt(n)))


Heatplot_RR_MAY <- Resp_means_MAY %>% 
  dplyr::mutate(OA_Sal = (paste(pCO2, Salinity,sep=''))) %>% 
  ggplot(aes(x = as.factor(Age),
             y = mean_RR)) +
  #shape = stage)) + 
  geom_rect(aes(fill = mean_RR), 
            xmin = -Inf, 
            xmax = Inf, 
            ymin = -Inf, 
            ymax = Inf, 
            alpha = 0.3) +
  geom_point(color = 'black') +
  geom_errorbar(aes(ymin = mean_RR - se_RR, 
                    ymax = mean_RR + se_RR), 
                width = 0.5, 
                position= "dodge2") +
  theme_bw() +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_grid(vars(Temp), 
             vars(fct_relevel(OA_Sal, c("HL", "LL", "HH", "LH")))) +
  scale_fill_gradient(low = "orange", 
                      high = "forestgreen", 
                      aesthetics = "fill")

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/")
pdf("Respiration/Day8_RR_heatplot.pdf", width=8, height=6)
print(Heatplot_RR_MAY)
dev.off()


