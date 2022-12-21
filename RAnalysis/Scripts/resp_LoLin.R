# Purpose: Oyster project 2021 - Respiration rate data 
# measure respiration rate from raw Loligo output data 
# using Lolin.R (Olito etal. 201?) for reproducible and non-bias calculation of respiration rates

# Written by: Sam J Gurr (last edit 9/14/2021)

# LOAD PACKAGES :::::::::::::::::::::::::::::::::::::::::::::::::::::::

library(devtools) # devtools::install_github # use devtools to instlal github link
library(LoLinR) # install_github('colin-olito/LoLinR') # install LoLinR from github
library(dplyr)
library(lubridate)
library(rMR)

# SET WORKING DIRECTORY :::::::::::::::::::::::::::::::::::::::::::::::

setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")

# CHANGE THE FOLLOWING ..THEN CONTROL A + ENTER ::::::::::::::::::::::

path.p    <- "Data/Respiration" #the location of all your respirometry files 
a         <- 0.4
ouputNAME <- "Output/Respiration/Cumulative_resp_alpha0.4.csv" 

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# bring in the respiration file names 
file.names.full        <- basename(list.files(path = path.p, pattern = "txt$", recursive = TRUE)) #list all csv file names in the folder and subfolders
file.names             <- file.names.full[c(3:4,6:7,10:11, 13,15, 19:20)] # call all target files (excluding tests and empty files) - open in notepad++ to make sure you have the correct files
file.names # check
file.names.table       <- data.frame(file.names)
file.names.table$about <- c('20210319_new_sensor_7', '20210319', '20210430_LOWtemp_HIGHsal', 
                            '20210430_LOWtemp_LOWsal', '20210430_raw', '20210507_HIGHtemp_LOWsal',
                            '20210507_LOWtemp_LOWsal', '20210507_HIGHtemp_HIGHsal', '20210507_LOWtemp_HIGHsal_b', '20210507_LOWtemp_HIGHsal')

# ANALYSIS  :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Objective: use LoLinr to run all respiration rates in a non-bias and autonomous fashion
# Ouputs: there will be two main resources from this script: 
#   (1) cumulative spreasheet of all respiration rate value for each Channel on each day
#   (2) a folder of plots from the LoLinR script to a plots folder - this will allow troubleshooting and sanity checks 

# I. Call the cumulative dataframe that we will write to in the for loop below
df_total             <- data.frame() # start dataframe 
resp.table           <- data.frame(matrix(nrow = 1, ncol = 8)) # create dataframe to save cumunalitively during for loop
colnames(resp.table) <- c('Date', 'Channel', 'Lpc', 'Leq' , 'Lz', 'alpha', 'Notes', 'Filename') # names for comuns in the for loop

# II.A Loop #1 - call filenames one at a time for analysis
for(i in 1:nrow(file.names.table)) { # for every file in list start at the first and run this following function
  Resp.Data           <- read.delim2(file = paste(path.p,'/',file.names.table[i,1], sep=''), header = TRUE,skip = 37) #reads in the data files
  Resp.Data$date      <- paste((sub("2021.*", "", Resp.Data$Date..Time..DD.MM.YYYY.HH.MM.SS.)), '2021', sep='') #  date - use 'sub' to call everything before 2021, add back 2021 using paste
  Resp.Data$time_Sec  <- period_to_seconds(hms(substr((strptime(sub(".*2021/", "", Resp.Data$Date..Time..DD.MM.YYYY.HH.MM.SS.), "%I:%M:%S %p")) , 12,19))) # time - use 'sub' to call target time of the raw date time after 'year/' + strptime' convert to 24 hr clock + 'period_to_seconds' converts the hms to seconds  
  Resp.Data$seconds   <- (Resp.Data$time_Sec - Resp.Data$time_Sec[1])    # secs - calc the sec time series
  Resp.Data$minutes   <- (Resp.Data$time_Sec - Resp.Data$time_Sec[1])/60 # mins - calc the minute time series
  
  temperature_C       <- as.numeric(Resp.Data$CH1.temp...C.[1])
  barromP_kPa         <- as.numeric(Resp.Data$Barometric.pressure..hPa.[1]) / 10
  salinity.pp.thou    <- as.numeric(Resp.Data$Salinity....[1])
  
  Resp.Data           <- Resp.Data %>% # use 'dplyr' 
                          dplyr::filter(!Phase %in% 'Flush') %>% # remove the initial rows labeled flush
                          dplyr::select(c(date, seconds, minutes, contains(".O2...air.sat"))) # call unique column names for the 8 Channels
  colnames(Resp.Data)[c(4:(ncol(Resp.Data)))] <- substr( ( colnames(Resp.Data)[c(4:(ncol(Resp.Data)))] ), 1,3) # clean these column names to make things easier - first 3 characters
  
  # Truncate!
  # the loligo recoreded values every second, this slows the model dramatically with >2000 values for each Channel, call every 30 seconds to speed this up
  # discuss with collaborators on this truncated approach 
  Resp.Data_30sec = Resp.Data[seq(1, nrow(Resp.Data), 30), ]
  
  
      # II.B Loop #2 - within each filename, call each Channel (% O2 sat.) to measure respiration rates, output plots and cumulative spreadsheet
      for(j in 4:(ncol(Resp.Data_30sec))){
        
        Resp_loop         <- na.omit(Resp.Data_30sec[,c(3,j)]) # noticed some random rows have 'NaN' - so I will loop the min and Channels to ommit Nas before proceeding
        Resp_loop$mgL     <- DO.unit.convert(as.numeric(Resp_loop[,2]),  # DO in percent air sat to be converted to mgL
                                             DO.units.in = "pct", DO.units.out ="mg/L", 
                                             bar.units.in = "kPa", bar.press = barromP_kPa, bar.units.out = "kpa",
                                             temp.C = temperature_C, 
                                             salinity.units = "pp.thou", salinity = salinity.pp.thou)
        
          
          
         if (nrow(Resp_loop) < 1) { 
            resp.table$Date                <- Resp.Data_30sec[1,1]
            resp.table$Channel             <- colnames(Resp_loop)[2] 
            resp.table[3:ncol(resp.table)] <- 'NA'
            df       <- data.frame(resp.table) # name dataframe for this single row
            df_total <- rbind(df_total,df) #bind to a cumulative list dataframe
            print(df_total) # print to monitor progress
            
            } else {
                model <- rankLocReg(
                xall    = as.numeric(Resp_loop[, 1]), 
                yall    = as.numeric(Resp_loop[, 3]), # call x as the minute timeseries and y as the % air saturation of the particular Channel
                alpha   = a,  # alpha was assigned earlier as 0.4 by the authors default suggestions - review Olito et al. and their github page for details
                method  = "pc", 
                verbose = TRUE) 
            
                sum.table <- summary(model)
                
                resp.table$Date       <- Resp.Data_30sec[1,1]
                resp.table$Channel    <- colnames(Resp_loop)[2] 
                resp.table$Lpc        <- sum.table$Lcompare[3,6] # Lpc slope - call the column 'b1'
                resp.table$Leq        <- sum.table$Lcompare[2,6] # Leq slope - call the column 'b1'
                resp.table$Lz         <- sum.table$Lcompare[1,6] # Lz slope  - call the column 'b1'
                resp.table$alpha      <- a
                resp.table$Notes      <- file.names.table[i,2]
                resp.table$Filename   <- file.names.table[i,1]
                
                df       <- data.frame(resp.table) # name dataframe for this single row
                df_total <- rbind(df_total,df) #bind to a cumulative list dataframe
                print(df_total) # print to monitor progress
            
                # save plots every inside loop and name by date_run_vialposition
                pdf(paste0("C:/Users/samjg/Documents/Github_repositories/Cvriginica_multistressor/RAnalysis/Output/Respiration/plots_alpha0.4_increm30sec/",file.names.table[i,2],"_",colnames(Resp_loop)[2],"_regression.pdf"))
                plot(model)
                dev.off()
            } # end of if else statement
  } # end of inside for loop
} # end of outside for loop

# merge with the preexisiting table?
cumulative_resp_table <- read.csv(file=ouputNAME, header=T) #call the pre existing cumulative table
new_table             <- rbind(cumulative_resp_table, df_total) # bind the new table from the for loop to the pre exisiting table
write.table(new_table,ouputNAME,sep=",", row.names=FALSE)  # write out to the path names outputNAME





# AFTER VISUAL INSPECTION OF PLOTS....
# we have the following rates that need to be rerun... (check the diagnostic plots and see for yourself!)

# (1) 1_3_19_21_new_sensor_for_7  CH2 - problem is the only the first chunk of data was pulled
# - solution = call more data to represent the O2 consmption

# (2) # (1) 1_3_19_21_new_sensor_for_7   CH3 - problem is the only the first chunk of data was pulled
# - solution = call more data to represent the O2 consmption 



#is the data a txt file? (from Lolin 8 channel
setwd("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis")
resp_rerun           <- read.delim2(file = "Data/Respiration/1_3_19_21_new_sensor_for_7_raw.txt", header = TRUE, skip = 37, fileEncoding = 'latin1')
resp_rerun$date      <- paste((sub("2021.*", "", resp_rerun$Date..Time..DD.MM.YYYY.HH.MM.SS.)), '2021', sep='') #  date - use 'sub' to call everything before 2021, add back 2021 using paste
resp_rerun$time_Sec  <- period_to_seconds(hms(substr((strptime(sub(".*2021/", "", resp_rerun$Date..Time..DD.MM.YYYY.HH.MM.SS.), "%I:%M:%S %p")) , 12,19))) # time - use 'sub' to call target time of the raw date time after 'year/' + strptime' convert to 24 hr clock + 'period_to_seconds' converts the hms to seconds  
resp_rerun$seconds   <- (resp_rerun$time_Sec - resp_rerun$time_Sec[1])    # secs - calc the sec time series
resp_rerun$minutes   <- (resp_rerun$time_Sec - resp_rerun$time_Sec[1])/60 # mins - calc the minute time series

#to calculate mg per L from air saturation....
temperature_C        <- as.numeric(resp_rerun$CH1.temp...C.[1])
barromP_kPa          <- as.numeric(resp_rerun$Barometric.pressure..hPa.[1]) / 10
salinity.pp.thou     <- as.numeric(resp_rerun$Salinity....[1])

resp_rerun_LoLin <- resp_rerun[seq(1, nrow(resp_rerun), 15), ]  %>% # data every 15 seconds to decrease the run time
  dplyr::filter(!colnames(resp_rerun)[2] %in% 'NaN') %>% # Lolin recorede NAs are written as 'Nan' - wonts run unless removed!
  #dplyr::select(c("minutes", "CH2.O2...air.sat..")) %>% # CH2
  dplyr::select(c("minutes", "CH3.O2...air.sat..")) %>% # CH3
  #dplyr::select(c("minutes", "CH3.O2...air.sat..")) %>% # run 2 ch3 20221026
  dplyr::filter(minutes >25) %>% # run 2 ch7 20221026 - before numute 20
  dplyr::mutate(mgL = (DO.unit.convert(as.numeric(CH3.O2...air.sat..),  # DO in percent air sat to be converted to mgL - uses an R package from loligo rMR
                                       DO.units.in = "pct", DO.units.out ="mg/L", 
                                       bar.units.in = "kPa", bar.press = barromP_kPa, bar.units.out = "kpa",
                                       temp.C = temperature_C, 
                                       salinity.units = "pp.thou", salinity = salinity.pp.thou)))

model <- rankLocReg(
  xall    = as.numeric(resp_rerun_LoLin[, 1]), 
  yall    = as.numeric(resp_rerun_LoLin[, 3]), # call x as the minute timeseries and y as the mg L-1 O2 
  alpha   = 0.4,  # alpha was assigned earlier as 0.4 by the authors default suggestions - review Olito et al. and their github page for details
  method  = "pc", 
  verbose = TRUE) 
plot(model) # CH2; Lpc == -0.0046
plot(model) # CH3; Lpc == -0.0072
