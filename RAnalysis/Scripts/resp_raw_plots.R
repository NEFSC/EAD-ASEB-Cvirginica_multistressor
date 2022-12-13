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

# LOAD DATA :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# bring in the respiration file names 
file.names.full        <- basename(list.files(path = path.p, pattern = "txt$", recursive = TRUE)) #list all csv file names in the folder and subfolders
file.names             <- file.names.full[c(3:4,6:7,10:11, 13,15, 19:20)] # call all target files (excluding tests and empty files) - open in notepad++ to make sure you have the correct files
file.names # check
file.names.table       <- data.frame(file.names)
file.names.table$about <- c('20210319_new_sensor_7', '20210319', '20210430_LOWtemp_HIGHsal', 
                            '20210430_LOWtemp_LOWsal', '20210430_raw', '20210507_HIGHtemp_LOWsal',
                            '20210507_LOWtemp_LOWsal', '20210507_HIGHtemp_HIGHsal', '20210507_LOWtemp_HIGHsal_b', '20210507_LOWtemp_HIGHsal')


# II.A Loop #1 - call filenames one at a time for analysis
for(i in 1:nrow(file.names.table)) { # for every file in list start at the first and run this following function
  Resp.Data           <- read.delim2(file = paste(path.p,'/',file.names.table[i,1], sep=''), header = TRUE,skip = 37, fileEncoding="latin1") #reads in the data files
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
  
  plot_title   <- gsub("_raw.*","", file.names.table[i,1])

  
  PLOT <- Resp.Data_30sec %>% 
    dplyr::select(-c('date', 'seconds')) %>%  
    reshape2::melt(id.vars = "minutes",variable.name = "channel", value.name = "air.sat") %>%
    dplyr::filter(!air.sat  %in% 'NaN') %>% 
    dplyr::mutate(mg.L.min =   (DO.unit.convert(as.numeric(air.sat),  # DO in percent air sat to be converted to mgL - uses an R package from loligo rMR
                                                DO.units.in = "pct", DO.units.out ="mg/L", 
                                                bar.units.in = "kPa", bar.press = barromP_kPa, bar.units.out = "kpa",
                                                temp.C = temperature_C, 
                                                salinity.units = "pp.thou", salinity = salinity.pp.thou))) %>% 
    ggplot(aes(x = minutes , y = mg.L.min)) +
    geom_smooth(method = "loess", se=FALSE, color="black", formula = mg.L.min ~ minutes) +
    theme_classic() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none", plot.title = element_text(size=10))+ 
    labs(y = expression(RAW_mg~L^{-1}~O[2]%.%min^{-1})) +
    xlab("minutes") +
    geom_point() + 
    ggtitle(plot_title) +
    facet_wrap(~channel)
  
  #pdf(paste0("C:/Users/samuel.gurr/Documents/Github_repositories/Airradians_multigen_OA/RAnalysis/Output/Respiration/plots_raw/",folder.names.table[i,1],"_", sub("_raw.*","",file.names.table[m,1]),"_regression.pdf"), width=10, height=12)
  pdf(paste0("C:/Users/samjg/Documents/Github_repositories/Cvirginica_multistressor/RAnalysis/Output/Respiration/raw_plots/",gsub("_raw.*","", file.names.table[i,1]),"_regression.pdf"), width=10, height=12)
  print(PLOT)
  dev.off()
}


