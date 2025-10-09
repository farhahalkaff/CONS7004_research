## GET NECESSARY STUFF ==============================================================================
# Set the directory where your files are located
data_directory <- "dataset"

# Get a list of all files with a specific extension (e.g., .csv)
files <- list.files(path = data_directory, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("Helgoland_2002.csv", files)]
files <- files[!grepl("Helgoland_1998.csv", files)]

# Use lapply to read each file into a list of data frames
# Replace read.csv with the appropriate function for your file type (e.g., read.delim, read_excel from 'readxl' package)
list_of_dataframes <- lapply(files, read.csv)

library(dplyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(tidyr)
library(gamlss)
library(gamlss.dist)
library(gamlss.ggplots)
library(modeest)
library(e1071)
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(forecast)        # Construct fourier terms for time series
library(gratia)          # Graceful plotting of smooth terms
library(marginaleffects) # Interrogating regression models
library(janitor)         # Creating clean, tidy variable names
library(cowplot)
library(patchwork)

df <- files %>%
  lapply(read.csv) %>%
  bind_rows()


# add year as a column
df$year <- year(ymd_hm(df$Date))
df <- df %>%
  mutate(Date = as.Date(ymd_hm(Date)))

# add month as a column 
df <- df %>% 
  mutate(month = factor(format(Date, "%m")))

# subset Ammonium 
Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium, year, month)

# subset Nitrite 
Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite, year, month)


library(gamlss.tr)
gen.trun(0,"NO",type="left")
gen.trun(0,"NO2",type="left")
gen.trun(0,"GU",type="left")
gen.trun(0,"LO",type="left")
gen.trun(0,"RG",type="left")
gen.trun(0,"exGAUS",type="left")
gen.trun(0,"NOF",type="left")
gen.trun(0,"PE",type="left")
gen.trun(0,"PE2",type="left")
gen.trun(0,"SN1",type="left")
gen.trun(0,"SN2",type="left")
gen.trun(0,"TF",type="left")
gen.trun(0,"TF2",type="left")
gen.trun(0,"GT",type="left")
gen.trun(0,"JSU",type="left")
gen.trun(0,"JSUo",type="left")
gen.trun(0,"NET",type="left")
gen.trun(0,"SHASH",type="left")
gen.trun(0,"SHASHo",type="left")
gen.trun(0,"SHASHo2",type="left")
gen.trun(0,"SEP1",type="left")
gen.trun(0,"SEP2",type="left")
gen.trun(0,"SEP3",type="left")
gen.trun(0,"SEP4",type="left")
gen.trun(0,"SST",type="left")
gen.trun(0,"ST1",type="left")
gen.trun(0,"ST2",type="left")
gen.trun(0,"ST3",type="left")
gen.trun(0,"ST4",type="left")
gen.trun(0,"ST5",type="left")
gen.trun(0,"EGB2",type="left")
# ========================================================================================================

## RQ2: Are the higher moments changing through time?? 

##############
## Ammonium ##
##############

# continuous distribution scaled (0, infinity) -----------------------------
## one-parameter (EXP)
am_meanonly_EXP <- gamlss(Ammonium ~ year + month, family = EXP(), data = Ammonium) 
## two-parameter (GA)
am_meanonly_GA <- gamlss(Ammonium ~ year + month, family = GA(), data = Ammonium) 
# IG
am_meanonly_IG <- gamlss(Ammonium ~ year + month, family = IG(), data = Ammonium) 
# LOGNO
am_meanonly_LOGNO <- gamlss(Ammonium ~ year + month, family = LOGNO(), data = Ammonium) 
# LOGNO2
am_meanonly_LOGNO2 <- gamlss(Ammonium ~ year + month, family = LOGNO2(), data = Ammonium) 
# PARETI1o
am_meanonly_PARETO1o <- gamlss(Ammonium ~ year + month, family = PARETO1o(), data = Ammonium) 
# PARETI2o
am_meanonly_PARETO2o <- gamlss(Ammonium ~ year + month, family = PARETO2o(), data = Ammonium,
                               #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# PARETO2
am_meanonly_PARETO2 <- gamlss(Ammonium ~ year + month, family = PARETO2(), data = Ammonium,
                               #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# WEI
am_meanonly_WEI <- gamlss(Ammonium ~ year + month, family = WEI(), data = Ammonium) 
# WEI2
am_meanonly_WEI2 <- gamlss(Ammonium ~ year + month, family = WEI2(), data = Ammonium,
                           #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# WEI3
am_meanonly_WEI3 <- gamlss(Ammonium ~ year + month, family = WEI3(), data = Ammonium) 
# three-parameter (BCCGo)
am_meanonly_BCCGo <- gamlss(Ammonium ~ year + month, family = BCCGo(), data = Ammonium) 
# BCCG
am_meanonly_BCCG <- gamlss(Ammonium ~ year + month, family = BCCG(), data = Ammonium) 
# GAF
am_meanonly_GAF <- gamlss(Ammonium ~ year + month, family = GAF(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# GG
am_meanonly_GG <- gamlss(Ammonium ~ year + month, family = GG(), data = Ammonium) 
# GIG
am_meanonly_GIG <- gamlss(Ammonium ~ year + month, family = GIG(), data = Ammonium) 
# four-parameter (BCT)
am_meanonly_BCT <- gamlss(Ammonium ~ year + month, family = BCT(), data = Ammonium) 
# BCTo
am_meanonly_BCTo <- gamlss(Ammonium ~ year + month, family = BCTo(), data = Ammonium) 
# BCPEo
am_meanonly_BCPEo <- gamlss(Ammonium ~ year + month, family = BCPEo(), data = Ammonium) 
# BCPE
am_meanonly_BCPE <- gamlss(Ammonium ~ year + month, family = BCPE(), data = Ammonium) 
# GB2 
am_meanonly_GB2 <- gamlss(Ammonium ~ year + month, family = GB2(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 



# continuous distribution left-truncate (0, infinity) -----------------------------
## two-parameter (NO)
am_meanonly_NOtr <- gamlss(Ammonium ~ year + month, family = NOtr(), data = Ammonium) 
##  NO2
am_meanonly_NO2tr <- gamlss(Ammonium ~ year + month, family = NO2tr(), data = Ammonium) 
## GU
am_meanonly_GUtr <- gamlss(Ammonium ~ year + month, family = GUtr(), data = Ammonium) 
## LO
am_meanonly_LOtr <- gamlss(Ammonium ~ year + month, family = LOtr(), data = Ammonium) 
## RG
am_meanonly_RGtr <- gamlss(Ammonium ~ year + month, family = RGtr(), data = Ammonium) 
## three-parameter (exGAUS)
am_meanonly_exGAUStr <- gamlss(Ammonium ~ year + month, family = exGAUStr(), data = Ammonium,
                               #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## NOF
am_meanonly_NOFtr <- gamlss(Ammonium ~ year + month, family = NOFtr(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## PE
am_meanonly_PEtr <- gamlss(Ammonium ~ year + month, family = PEtr(), data = Ammonium) 
## PE2
am_meanonly_PE2tr <- gamlss(Ammonium ~ year + month, family = PE2tr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SN1
am_meanonly_SN1tr <- gamlss(Ammonium ~ year + month, family = SN1tr(), data = Ammonium)    # took too damn long
## SN2
am_meanonly_SN2tr <- gamlss(Ammonium ~ year + month, family = SN2tr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## TF
am_meanonly_TFtr <- gamlss(Ammonium ~ year + month, family = TFtr(), data = Ammonium) 
## TF2
am_meanonly_TF2tr <- gamlss(Ammonium ~ year + month, family = TF2tr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## four_parameter (EGB2)
am_meanonly_EGB2tr <- gamlss(Ammonium ~ year + month, family = EGB2tr(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## GT
am_meanonly_GTtr <- gamlss(Ammonium ~ year + month, family = GTtr(), data = Ammonium)   # took too damn long
## JSU
am_meanonly_JSUtr <- gamlss(Ammonium ~ year + month, family = JSUtr(), data = Ammonium)   # NA's in the working vector or weights for parameter nu
## JSUo
am_meanonly_JSUotr <- gamlss(Ammonium ~ year + month, family = JSUotr(), data = Ammonium)  # NA's in the working vector or weights for parameter nu
## NET
am_meanonly_NETtr <- gamlss(Ammonium ~ year + month, family = NETtr(), data = Ammonium) 
## SHASH
am_meanonly_SHASHtr <- gamlss(Ammonium ~ year + month, family = SHASHtr(), data = Ammonium,
                              #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SHASHo
am_meanonly_SHASHotr <- gamlss(Ammonium ~ year + month, family = SHASHotr(), data = Ammonium,   # missing value where TRUE/FALSE needed
                               #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SHASHo2
am_meanonly_SHASHo2tr <- gamlss(Ammonium ~ year + month, family = SHASHo2tr(), data = Ammonium) 
## SEP1
am_meanonly_SEP1tr <- gamlss(Ammonium ~ year + month, family = SEP1tr(), data = Ammonium)  # took too damn long
## SEP2
am_meanonly_SEP2tr <- gamlss(Ammonium ~ year + month, family = SEP2tr(), data = Ammonium)  # took too damn long
## SEP3
am_meanonly_SEP3tr <- gamlss(Ammonium ~ year + month, family = SEP3tr(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SEP4
am_meanonly_SEP4tr <- gamlss(Ammonium ~ year + month, family = SEP4tr(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SST
am_meanonly_SSTtr <- gamlss(Ammonium ~ year + month, family = SSTtr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## ST1
am_meanonly_ST1tr <- gamlss(Ammonium ~ year + month, family = ST1tr(), data = Ammonium)  # took too damn long
## ST2
am_meanonly_ST2tr <- gamlss(Ammonium ~ year + month, family = ST2tr(), data = Ammonium,
                            mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))   # took too damn long
## ST3
am_meanonly_ST3tr <- gamlss(Ammonium ~ year + month, family = ST3tr(), data = Ammonium,
                            mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## ST4
am_meanonly_ST4tr <- gamlss(Ammonium ~ year + month, family = ST4tr(), data = Ammonium,
                            mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 
## ST5
am_meanonly_ST5tr <- gamlss(Ammonium ~ year + month, family = ST5tr(), data = Ammonium,
                            mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 




# continuous distribution  (-infinity, infinity) -----------------------------
## two-parameter (NO)
am_meanonly_NO <- gamlss(Ammonium ~ year + month, family = NO(), data = Ammonium) 
##  NO2
am_meanonly_NO2 <- gamlss(Ammonium ~ year + month, family = NO2(), data = Ammonium) 
## GU
am_meanonly_GU <- gamlss(Ammonium ~ year + month, family = GU(), data = Ammonium) 
## LO
am_meanonly_LO <- gamlss(Ammonium ~ year + month, family = LO(), data = Ammonium) 
## RG
am_meanonly_RG <- gamlss(Ammonium ~ year + month, family = RG(), data = Ammonium) 
## three-parameter (exGAUS)
am_meanonly_exGAUS <- gamlss(Ammonium ~ year + month, family = exGAUS(), data = Ammonium,
                               #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## NOF
am_meanonly_NOF <- gamlss(Ammonium ~ year + month, family = NOF(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## PE
am_meanonly_PE <- gamlss(Ammonium ~ year + month, family = PE(), data = Ammonium) 
## PE2
am_meanonly_PE2 <- gamlss(Ammonium ~ year + month, family = PE2(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SN1
am_meanonly_SN1 <- gamlss(Ammonium ~ year + month, family = SN1(), data = Ammonium)    
## SN2
am_meanonly_SN2 <- gamlss(Ammonium ~ year + month, family = SN2(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## TF
am_meanonly_TF <- gamlss(Ammonium ~ year + month, family = TF(), data = Ammonium) 
## TF2
am_meanonly_TF2 <- gamlss(Ammonium ~ year + month, family = TF2(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## four_parameter (EGB2)
am_meanonly_EGB2 <- gamlss(Ammonium ~ year + month, family = EGB2(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## GT
am_meanonly_GT <- gamlss(Ammonium ~ year + month, family = GT(), data = Ammonium,
                         #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))   
## JSU
am_meanonly_JSU <- gamlss(Ammonium ~ year + month, family = JSU(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))   
## JSUo
am_meanonly_JSUo <- gamlss(Ammonium ~ year + month, family = JSUo(), data = Ammonium,
                           #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
## NET
am_meanonly_NET <- gamlss(Ammonium ~ year + month, family = NET(), data = Ammonium) 
## SHASH
am_meanonly_SHASH <- gamlss(Ammonium ~ year + month, family = SHASH(), data = Ammonium,
                              #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SHASHo
am_meanonly_SHASHo <- gamlss(Ammonium ~ year + month, family = SHASHo(), data = Ammonium,   
                               #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SHASHo2
am_meanonly_SHASHo2 <- gamlss(Ammonium ~ year + month, family = SHASHo2(), data = Ammonium) 
## SEP1
am_meanonly_SEP1 <- gamlss(Ammonium ~ year + month, family = SEP1(), data = Ammonium,   # took too damn long
                           #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SEP2
am_meanonly_SEP2 <- gamlss(Ammonium ~ year + month, family = SEP2(), data = Ammonium,
                           #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SEP3
am_meanonly_SEP3 <- gamlss(Ammonium ~ year + month, family = SEP3(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SEP4
am_meanonly_SEP4 <- gamlss(Ammonium ~ year + month, family = SEP4(), data = Ammonium,
                             #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## SST
am_meanonly_SST <- gamlss(Ammonium ~ year + month, family = SST(), data = Ammonium) 
## ST1
am_meanonly_ST1 <- gamlss(Ammonium ~ year + month, family = ST1(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## ST2
am_meanonly_ST2 <- gamlss(Ammonium ~ year + month, family = ST2(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## ST3
am_meanonly_ST3 <- gamlss(Ammonium ~ year + month, family = ST3(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## ST4
am_meanonly_ST4 <- gamlss(Ammonium ~ year + month, family = ST4(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
## ST5
am_meanonly_ST5 <- gamlss(Ammonium ~ year + month, family = ST5(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 




# rank all of them models based on AIC ===============================================================

am_meanonly_AIC_combine <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "EGB2","GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "SEP4" ,"SST", "ST1", "ST2", "ST3", "ST4", "ST5",
             
             "EXP", 
             "GA", "IG", "LOGNO", "LOGNO2", "PARETO1o", "PARENTO2o", "PARETO2", "WEI", "WEI2", "WEI3", 
             "BCCGo", "BCCG","GAF", "GG", "GIG",
             "BCT", "BCTo", "BCPEo", "BCPE","GB2",
             
             "NOtr", "NO2tr", "GUtr", "LOtr", "RGtr",
             "exGAUStr", "NOFtr", "PEtr", "PE2tr","SN2tr","TFtr", "TF2tr", 
             "EGB2tr", "NETtr", "SHASHtr", "SHASHo2",
             "SEP3tr", "SEP4tr" ,"SSTtr", "ST3tr", "ST4tr", "ST5tr"),
  AIC = c(GAIC(am_meanonly_NO, k = 2),GAIC(am_meanonly_NO2, k = 2), GAIC(am_meanonly_GU, k = 2), GAIC(am_meanonly_LO, k = 2), GAIC(am_meanonly_RG, k = 2),
          GAIC(am_meanonly_exGAUS, k = 2), GAIC(am_meanonly_NOF, k = 2), GAIC(am_meanonly_PE, k = 2), GAIC(am_meanonly_PE2, k = 2),GAIC(am_meanonly_SN1, k = 2), 
          GAIC(am_meanonly_SN2, k = 2),GAIC(am_meanonly_TF, k = 2), GAIC(am_meanonly_TF2, k = 2), GAIC(am_meanonly_EGB2, k = 2),
          GAIC(am_meanonly_GT, k = 2),GAIC(am_meanonly_JSU, k = 2), GAIC(am_meanonly_JSUo, k = 2), GAIC(am_meanonly_NET, k = 2),
          GAIC(am_meanonly_SHASH, k = 2),GAIC(am_meanonly_SHASHo, k = 2), GAIC(am_meanonly_SHASHo2, k = 2), GAIC(am_meanonly_SEP2, k = 2), GAIC(am_meanonly_SEP3, k = 2),
          GAIC(am_meanonly_SEP4, k = 2), GAIC(am_meanonly_SST, k = 2), GAIC(am_meanonly_ST1, k = 2),
          GAIC(am_meanonly_ST2, k = 2),GAIC(am_meanonly_ST3, k = 2), GAIC(am_meanonly_ST4, k = 2), GAIC(am_meanonly_ST5, k = 2),
          
          GAIC(am_meanonly_EXP, k = 2),
          GAIC(am_meanonly_GA, k = 2), GAIC(am_meanonly_IG, k = 2), GAIC(am_meanonly_LOGNO, k = 2),GAIC(am_meanonly_LOGNO2, k = 2), GAIC(am_meanonly_PARETO1o, k = 2),
          GAIC(am_meanonly_PARETO2o, k = 2), GAIC(am_meanonly_PARETO2, k = 2), GAIC(am_meanonly_WEI, k = 2), GAIC(am_meanonly_WEI2, k = 2),GAIC(am_meanonly_WEI3, k = 2), 
          GAIC(am_meanonly_BCCGo, k = 2),GAIC(am_meanonly_BCCG, k = 2),GAIC(am_meanonly_GAF, k = 2), GAIC(am_meanonly_GG, k = 2),GAIC(am_meanonly_GIG, k = 2),
          GAIC(am_meanonly_BCT, k = 2),GAIC(am_meanonly_BCTo, k = 2), GAIC(am_meanonly_BCPEo, k = 2), GAIC(am_meanonly_BCPE, k = 2),GAIC(am_meanonly_GB2, k = 2),
          
          GAIC(am_meanonly_NOtr, k = 2),GAIC(am_meanonly_NO2tr, k = 2), GAIC(am_meanonly_GUtr, k = 2), GAIC(am_meanonly_LOtr, k = 2), GAIC(am_meanonly_RGtr, k = 2),
          GAIC(am_meanonly_exGAUStr, k = 2), GAIC(am_meanonly_NOFtr, k = 2), GAIC(am_meanonly_PEtr, k = 2), GAIC(am_meanonly_PE2tr, k = 2),
          GAIC(am_meanonly_SN2tr, k = 2),GAIC(am_meanonly_TFtr, k = 2), GAIC(am_meanonly_TF2tr, k = 2), GAIC(am_meanonly_EGB2tr, k = 2), GAIC(am_meanonly_NETtr, k = 2),
          GAIC(am_meanonly_SHASHtr, k = 2), GAIC(am_meanonly_SHASHo2tr, k = 2),  GAIC(am_meanonly_SEP3tr, k = 2),
          GAIC(am_meanonly_SEP4tr, k = 2), GAIC(am_meanonly_SSTtr, k = 2), 
          GAIC(am_meanonly_ST3tr, k = 2), GAIC(am_meanonly_ST4tr, k = 2), GAIC(am_meanonly_ST5tr, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4","4",
            
            "1",
            "2", "2", "2", "2", "2", "2", "2", "2", "2", "2",
            "3","3","3","3","3",
            "4","4","4","4", "4", 
            
            "2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve", "sym", "sym", "sym", "both","both", "sym", "sym",
               "both","sym", "both", "both", "sym", "both", "both", "both",
               "both","both" ,"both","both", "both", "both", "both", "both", "both",
               
               "+ve", 
               "+ve", "+ve", "+ve", "+ve", "+ve", "+ve", "+ve", "both","both","both",
               "both","both", "+ve", "both","+ve",
               "both", "both", "both","both", "both",
               
               "sym", "sym", "-ve", "sym", "+ve", 
               "+ve", "sym", "sym", "sym", "both","sym", "sym",
               "both", "sym", "both", "both",
               "both","both","both","both","both","both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto", "meso", "both", "both", "meso","meso", "lepto", "lepto",
               "lepto","both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto", "lepto", "lepto", "lepto", "lepto", "lepto",
               
               "lepto",
               "lepto", "lepto", "lepto", "lepto","lepto", "lepto", "lepto","both","both","both",
               "both", "both","lepto", "both", "lepto",
               "lepto", "lepto", "both", "both","both", 
               
               "meso", "meso", "meso", "lepto", "lepto",
               "lepto", "meso", "both", "both", "meso","lepto", "lepto",
               "lepto", "lepto","both", "both",
               "both", "both", "lepto", "lepto","lepto","lepto")
)


pAIC_am_meanonly_AIC_combine <- ggplot(am_meanonly_AIC_combine, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "skew") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(x = "AIC", y = "mean only models (Ammonium)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


# check moment bucket and PIT histogram and worm plot  =======================================================
# of the best model of the three categories of fam distribution 

## moment bucket ########################################################
moment_bucket(am_meanonly_BCT, am_meanonly_EGB2tr, am_meanonly_SHASHo2) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

#------SHASHo2 and BCT does quite well. EGB2tr, not so much 

## PIT histogram #########################################################

# fitted param for (am_meanonly_BCT) -----------
am_mu_hat_BCT    <- predict(am_meanonly_BCT, "mu", type = "response")
am_sigma_hat_BCT <- predict(am_meanonly_BCT, "sigma", type = "response")
am_nu_hat_BCT    <- predict(am_meanonly_BCT, "nu", type = "response")
am_tau_hat_BCT   <- predict(am_meanonly_BCT, "tau", type = "response")

am_pit1 <- pBCT(Ammonium$Ammonium, mu = am_mu_hat_BCT, sigma = am_sigma_hat_BCT, 
              nu = am_nu_hat_BCT, tau = am_tau_hat_BCT)

# Plot PIT histogram
hist(am_pit1, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit1)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


# fitted param for (am_meanonly_EGB2tr) -----------
am_mu_hat_EGB2tr    <- predict(am_meanonly_EGB2tr, "mu", type = "response")
am_sigma_hat_EGB2tr <- predict(am_meanonly_EGB2tr, "sigma", type = "response")
am_nu_hat_EGB2tr    <- predict(am_meanonly_EGB2tr, "nu", type = "response")
am_tau_hat_EGB2tr   <- predict(am_meanonly_EGB2tr, "tau", type = "response")

am_pit2 <- pEGB2tr(Ammonium$Ammonium, mu = am_mu_hat_EGB2tr, sigma = am_sigma_hat_EGB2tr, 
                nu = am_nu_hat_EGB2tr, tau = am_tau_hat_EGB2tr)

# Plot PIT histogram
hist(am_pit2, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit2)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


# fitted param for (am_meanonly_SHASHo2) -----------
am_mu_hat_SHASHo2    <- predict(am_meanonly_SHASHo2, "mu", type = "response")
am_sigma_hat_SHASHo2 <- predict(am_meanonly_SHASHo2, "sigma", type = "response")
am_nu_hat_SHASHo2    <- predict(am_meanonly_SHASHo2, "nu", type = "response")
am_tau_hat_SHASHo2   <- predict(am_meanonly_SHASHo2, "tau", type = "response")

am_pit3 <- pSHASHo2(Ammonium$Ammonium, mu = am_mu_hat_SHASHo2, sigma = am_sigma_hat_SHASHo2, 
                   nu = am_nu_hat_SHASHo2, tau = am_tau_hat_SHASHo2)

# Plot PIT histogram
hist(am_pit3, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit3)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


##### The best looking histogram is BCT #######


# Worm plot ############################################################################

wp(am_meanonly_BCT, ylim.all=1)
wp(am_meanonly_EGB2tr, ylim.all=1)
wp(am_meanonly_SHASHo2, ylim.all=1)

# EGB2tr looks the horrendous, the other looks ok 

## concluding to use BCT as the final mode ##



# time-varying moments  ===============================================================

## mean only 
am_meanonly_BCT <- gamlss(Ammonium ~ year + month, family = BCT(), data = Ammonium) 

## mean sigma 
am_mean_and_sigma_BCT <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, 
                                family = BCT(), data = Ammonium) 

## no skew
am_noskew_BCT <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month, 
                                family = BCT(), data = Ammonium,
                        #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

## no kurt 
am_nokurt_BCT <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, 
                        family = BCT(), data = Ammonium,
                        #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# all param 
am_all_BCT <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = BCT(), data = Ammonium,
                        mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# check AIC 
AIC(am_meanonly_BCT)
AIC(am_mean_and_sigma_BCT)
AIC(am_noskew_BCT)
AIC(am_nokurt_BCT)
AIC(am_all_BCT) # Error in dBCT(y, mu, sigma, nu, tau, log = TRUE) : mu must be positive 


# I assume that this error occurs because the model is assuming the mu to be negative, but because
# BCT is a positive distribution, it can't take in any negative values, and comes out as error 
# ok so according to Robert from stack overflow, negative mu value can happen during iteration (especially
# when there are values of Y that is close to 0), and indeed the lowest Y for ammonium is 0.04, and there 
# are 157 observations that are below the value of 1
# the next best model is BCPE, which actually also have an identity mu link, but let's give it a try 

am_all_BCPE <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCPE(), data = Ammonium,
                     mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ok I got the same error, as expected 
# the next best thing is BCCG which I also expect to give the same error, let's give it a go 
am_all_BCCG <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCCG(), data = Ammonium,
                      mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ok the same error, as expected 
# let's best thing is GB2
am_all_GB2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = GB2(), data = Ammonium,
                      mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ok this works, because the mu is in a log link function
# now let's look at the rest of changing parameters and how it does 

# mean and sigma
am_mean_and_sigma_GB2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month,
                     family = GB2(), data = Ammonium,
                     #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

# no skew
am_noskew_GB2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                     family = GB2(), data = Ammonium,
                     mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# no kurt
am_nokurt_GB2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                     family = GB2(), data = Ammonium,
                     mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# hmm convergence issue for mean and sigma, no skew and no kurt



# now let's try SHASHo2

## mean only 
am_meanonly_SHASHo2 <- gamlss(Ammonium ~ year + month, family = SHASHo2(), data = Ammonium) 

## mean sigma 
am_mean_and_sigma_SHASHo2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, 
                                family = SHASHo2(), data = Ammonium,
                                #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                                method = mixed(10,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

## no skew
am_noskew_SHASHo2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month, 
                        family = SHASHo2(), data = Ammonium,
                        #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

## no kurt 
am_nokurt_SHASHo2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, 
                        family = SHASHo2(), data = Ammonium,
                        #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# all param 
am_all_SHASHo2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SHASHo2(), data = Ammonium,
                     #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# check AIC 
AIC(am_meanonly_SHASHo2)
AIC(am_mean_and_sigma_SHASHo2)
AIC(am_noskew_SHASHo2)
AIC(am_nokurt_SHASHo2)
AIC(am_all_SHASHo2) # tis is the best 



#### let's check moment bucket again 
# compare mean only and all param (time)
moment_bucket(am_all_GB2,am_meanonly_GB2) +
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

moment_bucket(am_all_SHASHo2,am_meanonly_SHASHo2) +
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# hmm why does it seem like mean only is doing better in terms of the skewness and kurtosis 

# lettuce check pit histogram 

# fitted param for (am_all_GB2) -----------
am_mu_hat_GB2    <- predict(am_all_GB2, "mu", type = "response")
am_sigma_hat_GB2 <- predict(am_all_GB2, "sigma", type = "response")
am_nu_hat_GB2    <- predict(am_all_GB2, "nu", type = "response")
am_tau_hat_GB2   <- predict(am_all_GB2, "tau", type = "response")

am_pit4 <- pGB2(Ammonium$Ammonium, mu = am_mu_hat_GB2, sigma = am_sigma_hat_GB2, 
                   nu = am_nu_hat_GB2, tau = am_tau_hat_GB2)

# Plot PIT histogram
hist(am_pit4, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit2)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity

# fitted param for (am_meanonly_GB2) -----------
am_mu_hat_meanonly_GB2    <- predict(am_meanonly_GB2, "mu", type = "response")
am_sigma_hat_meanonly_GB2 <- predict(am_meanonly_GB2, "sigma", type = "response")
am_nu_hat_meanonly_GB2    <- predict(am_meanonly_GB2, "nu", type = "response")
am_tau_hat_meanonly_GB2   <- predict(am_meanonly_GB2, "tau", type = "response")

am_pit4b <- pGB2(Ammonium$Ammonium, mu = am_mu_hat_meanonly_GB2, sigma = am_sigma_hat_meanonly_GB2, 
                nu = am_nu_hat_meanonly_GB2, tau = am_tau_hat_meanonly_GB2)

# Plot PIT histogram
hist(am_pit4b, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit2)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity

# fitted param for (am_all_SHASHo2) -----------
am_mu_hat_all_SHASHo2    <- predict(am_all_SHASHo2, "mu", type = "response")
am_sigma_hat_all_SHASHo2 <- predict(am_all_SHASHo2, "sigma", type = "response")
am_nu_hat_all_SHASHo2    <- predict(am_all_SHASHo2, "nu", type = "response")
am_tau_hat_all_SHASHo2   <- predict(am_all_SHASHo2, "tau", type = "response")

am_pit5 <- pSHASHo2(Ammonium$Ammonium, mu = am_mu_hat_all_SHASHo2, sigma = am_sigma_hat_all_SHASHo2, 
                nu = am_nu_hat_all_SHASHo2, tau = am_tau_hat_all_SHASHo2)

# Plot PIT histogram
hist(am_pit5, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit2)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


# check worm plot 

wp(am_all_SHASHo2, ylim.all=2)
wp(am_meanonly_GB2, ylim.all = 1)
wp(am_all_GB2, ylim.all = 1)

resid_wp(am_all_SHASHo2)
resid_wp(am_meanonly_BCT)
resid_wp(am_meanonly_EGB2tr)


resid_wp(am_all_GB2)

resid_wp(am_all_JSU)

resid_wp(am_all_SSTtr)

?kurtosis

# lemme add polynomial 
am_all_SHASHo2_poly <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                         family = SHASHo2(), data = Ammonium,
                         #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
AIC(am_all_SHASHo2_poly) # tis is better 


# lemme add polynomial 
am_all_GB2_poly <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                              family = GB2(), data = Ammonium,
                              #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
AIC(am_all_GB2_poly) # tis is better but convergence issue



# now lets plot and see how the models predict the mean 
Ammonium$am_mu_hat_all_SHASHo2 <- predict(am_all_SHASHo2, what = "mu", type = "response")
Ammonium$am_mu_hat_meanonly_SHASHo2 <- predict(am_meanonly_SHASHo2, what = "mu", type = "response")
Ammonium$am_mu_hat_all_GB2 <- predict(am_all_GB2, what = "mu", type = "response")
Ammonium$am_mu_hat_meanonly_GB2 <- predict(am_meanonly_GB2, what = "mu", type = "response")


# plot that sucka
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_all_SHASHo2), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_meanonly_SHASHo2), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = am_mu_hat_all_GB2), color = "goldenrod", linewidth = 1) +
  geom_line(aes(y = am_mu_hat_meanonly_GB2), color = "steelblue", linewidth = 1) +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()




## now let's get preidcted param from models 
# create function to get predicted param from model
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- predict(model, what = "mu", type = "response")
  data$sigma_hat <- predict(model, what = "sigma", type = "response")
  data$nu_hat <- predict(model, what = "nu", type = "response")
  data$tau_hat <- predict(model, what = "tau", type = "response")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
    summarise(mean_mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
              mean_nu_hat = mean(nu_hat, na.rm = TRUE),
              mean_tau_hat = mean(tau_hat, na.rm = TRUE))
  
  # get the average per month 
  average_by_month <- data %>%
    group_by(month) %>%
    summarise(mean_mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
              mean_nu_hat = mean(nu_hat, na.rm = TRUE),
              mean_tau_hat = mean(tau_hat, na.rm = TRUE))
  
  #return(list(by_year = average_by_year, by_month = average_by_month))
  return(average_by_year)
  #return(average_by_month)
  
}



best_model_param2 <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- predict(model, what = "mu", type = "link")
  data$sigma_hat <- predict(model, what = "sigma", type = "link")
  data$nu_hat <- predict(model, what = "nu", type = "link")
  data$tau_hat <- predict(model, what = "tau", type = "link")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
    summarise(mean_mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
              mean_nu_hat = mean(nu_hat, na.rm = TRUE),
              mean_tau_hat = mean(tau_hat, na.rm = TRUE))
  
  # get the average per month 
  average_by_month <- data %>%
    group_by(month) %>%
    summarise(mean_mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
              mean_nu_hat = mean(nu_hat, na.rm = TRUE),
              mean_tau_hat = mean(tau_hat, na.rm = TRUE))
  
  #return(list(by_year = average_by_year, by_month = average_by_month))
  return(average_by_year)
  #return(average_by_month)
  
}




# for GB2 all
am_all_GB2_param <- best_model_param(am_all_GB2, Ammonium)
# for GB2 mean only
am_meanonly_GB2_param <- best_model_param(am_meanonly_GB2, Ammonium)

# for SHASHo2 
am_allSHASHo2_param <- best_model_param(am_all_SHASHo2, Ammonium)
am_meanonly_SHASHo2_param <- best_model_param(am_meanonly_SHASHo2, Ammonium)

#### for the link function
# for GB2 all
am_all_GB2_param2 <- best_model_param2(am_all_GB2, Ammonium)
# for GB2 mean only
am_meanonly_GB2_param2 <- best_model_param2(am_meanonly_GB2, Ammonium)




# empirical estimate of moments
am_moments_summary <- Ammonium %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Ammonium, na.rm = TRUE),
    var = var(Ammonium, na.rm = TRUE),
    skew = skewness(Ammonium, na.rm = TRUE),
    kurt = kurtosis(Ammonium, na.rm = TRUE)
  )
print(ni_moments_summary, n=33)



# predicted mean
plot(am_moments_summary$mean ~ am_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(2,15))
lines(am_all_GB2_param$mean_mu_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(2,15))
#plot(am_all_GB2_param$mean_mu_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(2,15))
lines(am_meanonly_GB2_param$mean_mu_hat ~ am_meanonly_GB2_param$year, col = "pink", lwd = 2)
lines(am_allSHASHo2_param$mean_mu_hat ~ am_allSHASHo2_param$year, col = "chocolate", lwd = 2)
lines(am_meanonly_SHASHo2_param$mean_mu_hat ~ am_meanonly_SHASHo2_param$year, col = "goldenrod", lwd = 2)

# predicted variance
plot(am_moments_summary$var ~ am_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(0,22))
lines(am_all_GB2_param$mean_sigma_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(0,7))
#plot(am_all_GB2_param$mean_sigma_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(0,7))
lines(am_meanonly_GB2_param$mean_sigma_hat ~ am_meanonly_GB2_param$year, col = "pink", lwd = 2)
lines(am_allSHASHo2_param$mean_sigma_hat ~ am_allSHASHo2_param$year, col = "chocolate", lwd = 2)
lines(am_meanonly_SHASHo2_param$mean_sigma_hat ~ am_meanonly_SHASHo2_param$year, col = "goldenrod", lwd = 2)

# predicted skewness
plot(am_moments_summary$skew ~ am_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(-1,5))
lines(am_all_GB2_param$mean_nu_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(0,5))
#plot(am_all_GB2_param$mean_nu_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(0,5))
lines(am_meanonly_GB2_param$mean_nu_hat ~ am_meanonly_GB2_param$year, col = "pink", lwd = 2)
lines(am_allSHASHo2_param$mean_nu_hat ~ am_allSHASHo2_param$year, col = "chocolate", lwd = 2)
lines(am_meanonly_SHASHo2_param$mean_nu_hat ~ am_meanonly_SHASHo2_param$year, col = "goldenrod", lwd = 2)

# predicted tails
plot(am_moments_summary$kurt ~ am_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(-1,16))
lines(am_all_GB2_param$mean_tau_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(1,16))
#plot(am_all_GB2_param$mean_tau_hat ~ am_all_GB2_param$year, col = "darkred", lwd = 2, type = "l", ylim = c(1,16))
lines(am_meanonly_GB2_param$mean_tau_hat ~ am_meanonly_GB2_param$year, col = "pink", lwd = 2)
lines(am_allSHASHo2_param$mean_tau_hat ~ am_allSHASHo2_param$year, col = "chocolate", lwd = 2)
lines(am_meanonly_SHASHo2_param$mean_tau_hat ~ am_meanonly_SHASHo2_param$year, col = "goldenrod", lwd = 2)



####### PDF ================================================================================

# pdf of three random years from the data itself 
# 1964
Ammonium_1964 <- Ammonium[Ammonium$year %in% c(1964), ]
Ammonium_1964_pdf <- density(Ammonium_1964$Ammonium)
# 1975
Ammonium_1975 <- Ammonium[Ammonium$year %in% c(1975), ]
Ammonium_1975_pdf <- density(Ammonium_1975$Ammonium)
# 1992
Ammonium_1992 <- Ammonium[Ammonium$year %in% c(1992), ]
Ammonium_1992_pdf <- density(Ammonium_1992$Ammonium)


# get the pdf from models 
# get pdf for am_all_GB2 -----
x <- seq(min(am_all_GB2$y, na.rm = TRUE),
         max(am_all_GB2$y, na.rm = TRUE))

# natural scale (backtransform)
am_all_GB2_1964_pdf  <- dGB2(x, mu = 13.097801, sigma = 5.299032, nu = 1.2419962, tau = 2.469285)
am_all_GB2_1975_pdf  <- dGB2(x, mu = 11.147062, sigma = 3.310689, nu = 1.9729324, tau = 4.767189)
am_all_GB2_1992_pdf  <- dGB2(x, mu = 8.512500, sigma = 1.758140, nu = 4.1428046, tau = 13.094056)


# get pdf for am_meanonly_GB2 -----
x1 <- seq(min(am_meanonly_GB2$y, na.rm = TRUE),
         max(am_meanonly_GB2$y, na.rm = TRUE))

# natural scale (backtransform)
am_meanonly_GB2_1964_pdf  <- dGB2(x1, mu = 13.468603, sigma = 4.460873, nu = 0.5378375, tau = 1.188082)
am_meanonly_GB2_1975_pdf  <- dGB2(x1, mu = 8.394169, sigma = 4.460873, nu = 0.5378375, tau = 1.188082)
am_meanonly_GB2_1992_pdf  <- dGB2(x1, mu = 4.114894, sigma = 4.460873, nu = 0.5378375, tau = 1.188082)



# plot density plots for 1964 ------------------------------
plot(Ammonium_1964_pdf, lwd = 2, lty = 2, col = "black", xlim = c(0,25))
lines(am_all_GB2_1964_pdf, type = "l", lwd = 2, col = "goldenrod") # natural scale
lines(am_meanonly_GB2_1964_pdf, type = "l", lwd = 2, col = "steelblue") # natural scale 


# plot density plots for 1975 ------------------------------
plot(Ammonium_1975_pdf, ylim = c(0,0.18), lwd = 2, lty = 2, col = "black")
lines(am_all_GB2_1975_pdf, type = "l", lwd = 2, col = "goldenrod")
lines(am_meanonly_GB2_1975_pdf, type = "l", lwd = 2, col = "steelblue")


# plot density plots for 1992 ------------------------------
plot(Ammonium_1992_pdf, ylim = c(0,0.4), lwd = 2, lty = 2, col = "black")
lines(am_all_GB2_1992_pdf, type = "l", lwd = 2, col = "goldenrod")
lines(am_meanonly_GB2_1992_pdf, type = "l", lwd = 2, col = "steelblue")


mean(Ammonium_1964$Ammonium)

am_all_GB2_1992_pdf  <- dGB2(x, mu = 8.512500, sigma = 1.758140, nu = 4.1428046, tau = 13.094056)
am_meanonly_GB2_1992_pdf  <- dGB2(x1, mu = 4.114894, sigma = 4.460873, nu = 0.5378375, tau = 1.188082)





############# not lets do SSTtr since it's the next best thing for me use because I know how to get the parameters

## SSTtr
am_meanonly_SSTtr <- gamlss(Ammonium ~ year + month, family = SSTtr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

AIC(am_meanonly_SSTtr)
AIC(am_meanonly_SST)


# Check worm plot
wp(am_meanonly_SSTtr, ylim.all = 1)
wp(am_meanonly_SST, ylim.all = 1)



# time-varying param
# mean only
am_meanonly_SSTtr <- gamlss(Ammonium ~ year + month, family = SSTtr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# mean and sigma
am_mean_and_sigma_SSTtr <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, family = SSTtr(), data = Ammonium,
                            #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# no skew
am_noskew_SSTtr <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month, family = SSTtr(), data = Ammonium, # sigma must be positive
                                  #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                                  method = mixed(10,200),
                                  control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# no kurt
am_nokurt_SSTtr <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, family = SSTtr(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# all param
am_all_SSTtr <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                          family = SSTtr(), data = Ammonium,
                          #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 


AIC(am_meanonly_SSTtr)
AIC(am_mean_and_sigma_SSTtr)
AIC(am_nokurt_SSTtr)
AIC(am_all_SSTtr)


wp(am_all_SSTtr, ylim.all = 1)

moment_bucket(am_all_SSTtr, am_meanonly_SSTtr) +
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")



# all param SST
am_all_JSUo <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = JSUo(), data = Ammonium,
                       #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

AIC(am_all_JSUo)


wp(am_all_SSTtr, ylim.all = 1)
wp(am_all_JSU, ylim.all = 1)
wp(am_meanonly_SSTtr, ylim.all = 1)
wp(am_meanonly_JSU, ylim.all = 1)



# ok we are giong to compare SSTtr and JSU
# mean only
am_meanonly_JSU <- gamlss(Ammonium ~ year + month,
                      family = JSU(), data = Ammonium,
                      #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# mean and sigma
am_mean_and_sigma_JSU <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month,
                           family = JSU(), data = Ammonium,
                           #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# no skew
am_noskew_JSU <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                                 family = JSU(), data = Ammonium,
                                 #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                                 method = mixed(10,200),
                                 control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# no kurt
am_kurt_JSU <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                         family = JSU(), data = Ammonium,
                         #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# all param 
am_all_JSU <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = JSU(), data = Ammonium,
                      #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 


AIC(am_meanonly_JSU)
AIC(am_mean_and_sigma_JSU)
AIC(am_noskew_JSU)
AIC(am_kurt_JSU)
AIC(am_all_JSU)



#### PIT histogram 

# fitted param for (am_all_JSU) -----------
am_mu_hat_all_JSU    <- predict(am_all_JSU, "mu", type = "response")
am_sigma_hat_all_JSU <- predict(am_all_JSU, "sigma", type = "response")
am_nu_hat_all_JSU   <- predict(am_all_JSU, "nu", type = "response")
am_tau_hat_all_JSU   <- predict(am_all_JSU, "tau", type = "response")

am_pit6 <- pJSU(Ammonium$Ammonium, mu = am_mu_hat_all_JSU, sigma = am_sigma_hat_all_JSU, 
                    nu = am_nu_hat_all_JSU, tau = am_tau_hat_all_JSU)

# Plot PIT histogram
hist(am_pit6, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit6)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


# fitted param for (am_all_SSTtr) -----------
am_mu_hat_all_SSTtr    <- predict(am_all_SSTtr, "mu", type = "response")
am_sigma_hat_all_SSTtr <- predict(am_all_SSTtr, "sigma", type = "response")
am_nu_hat_all_SSTtr    <- predict(am_all_SSTtr, "nu", type = "response")
am_tau_hat_all_SSTtr   <- predict(am_all_SSTtr, "tau", type = "response")

am_pit7 <- pSSTtr(Ammonium$Ammonium, mu = am_mu_hat_all_SSTtr, sigma = am_sigma_hat_all_SSTtr, 
                 nu = am_nu_hat_all_SSTtr, tau = am_tau_hat_all_SSTtr)

# Plot PIT histogram
hist(am_pit7, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit7)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity




### how these models predict the mean

Ammonium$am_mu_hat_all_SSTtr <- predict(am_all_SSTtr, what = "mu", type = "response")
Ammonium$am_mu_hat_all_JSU <- predict(am_all_JSU, what = "mu", type = "response")
Ammonium$am_mu_hat_meanonly_SSTtr <- predict(am_meanonly_SSTtr, what = "mu", type = "response")
Ammonium$am_mu_hat_meanonly_JSU <- predict(am_meanonly_JSU, what = "mu", type = "response")

# plot that sucka
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_all_SSTtr), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_meanonly_SSTtr), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = am_mu_hat_all_JSU), color = "goldenrod", linewidth = 1) +
  geom_line(aes(y = am_mu_hat_meanonly_JSU), color = "steelblue", linewidth = 1) +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()





# pdf 
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- predict(model, what = "mu", type = "response")
  data$sigma_hat <- predict(model, what = "sigma", type = "response")
  data$nu_hat <- predict(model, what = "nu", type = "response")
  data$tau_hat <- predict(model, what = "tau", type = "response")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
    summarise(mean_mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
              mean_nu_hat = mean(nu_hat, na.rm = TRUE),
              mean_tau_hat = mean(tau_hat, na.rm = TRUE))
  
  # get the average per month 
  average_by_month <- data %>%
    group_by(month) %>%
    summarise(mean_mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_sigma_hat = mean(sigma_hat, na.rm = TRUE),
              mean_nu_hat = mean(nu_hat, na.rm = TRUE),
              mean_tau_hat = mean(tau_hat, na.rm = TRUE))
  
  #return(list(by_year = average_by_year, by_month = average_by_month))
  return(average_by_year)
  #return(average_by_month)
}



# for SSTtr all
am_all_SSTtr_param <- best_model_param(am_all_SSTtr, Ammonium)
# for SSTtr mean only
am_meanonly_SSTtr_param <- best_model_param(am_meanonly_SSTtr, Ammonium)

# for JSU all
am_all_JSU_param <- best_model_param(am_all_JSU, Ammonium)
# for SSTtr mean only
am_meanonly_JSU_param <- best_model_param(am_meanonly_JSU, Ammonium)


# pdf of three random years from the data itself ----------------
# 1964
Ammonium_1964 <- Ammonium[Ammonium$year %in% c(1964), ]
Ammonium_1964_pdf <- density(Ammonium_1964$Ammonium)
# 1975
Ammonium_1975 <- Ammonium[Ammonium$year %in% c(1975), ]
Ammonium_1975_pdf <- density(Ammonium_1975$Ammonium)
# 1992
Ammonium_1992 <- Ammonium[Ammonium$year %in% c(1992), ]
Ammonium_1992_pdf <- density(Ammonium_1992$Ammonium)


# get the pdf from models 
# get pdf for am_all_SSTtr -----
x <- seq(min(am_all_SSTtr$y, na.rm = TRUE),
         max(am_all_SSTtr$y, na.rm = TRUE))

# natural scale (backtransform)
am_all_SSTtr_1964_pdf  <- dSSTtr(x, mu = 9.343769, sigma = 3.266789, nu = 1.249053, tau = 30.254042)
am_all_SSTtr_1975_pdf  <- dSSTtr(x, mu = 6.769741, sigma = 4.092673, nu = 1.354003, tau = 5.610644)
am_all_SSTtr_1992_pdf  <- dSSTtr(x, mu = 2.997666, sigma = 5.685501, nu = 1.544101, tau = 2.124034)


# get pdf for am_meanonly_SSTtr -----
x1 <- seq(min(am_meanonly_SSTtr$y, na.rm = TRUE),
          max(am_meanonly_SSTtr$y, na.rm = TRUE))

# natural scale (backtransform)
am_meanonly_SSTtr_1964_pdf  <- dSSTtr(x1, mu = 9.437967, sigma = 3.328111, nu = 1.047589, tau = 3.38749)
am_meanonly_SSTtr_1975_pdf  <- dSSTtr(x1, mu = 6.720587, sigma = 3.328111, nu = 1.047589, tau = 3.38749)
am_meanonly_SSTtr_1992_pdf  <- dSSTtr(x1, mu = 2.655498, sigma = 3.328111, nu = 1.047589, tau = 3.38749)


# get pdf for am_all_JSU -----
x2 <- seq(min(am_all_JSU$y, na.rm = TRUE),
          max(am_all_JSU$y, na.rm = TRUE))

# natural scale (backtransform)
am_all_JSU_1964_pdf  <- dJSU(x2, mu = 9.415685, sigma = 3.738264, nu = 1.935820, tau = 6.065708)
am_all_JSU_1975_pdf  <- dJSU(x2, mu = 6.969798, sigma = 2.979428, nu = 2.843136, tau = 4.282139)
am_all_JSU_1992_pdf  <- dJSU(x2, mu = 3.297202, sigma = 2.057176, nu = 3.311038, tau = 1.999770)


# get pdf for am_meanonly_JSU -----
x3 <- seq(min(am_meanonly_JSU$y, na.rm = TRUE),
          max(am_meanonly_JSU$y, na.rm = TRUE))

# natural scale (backtransform)
am_meanonly_JSU_1964_pdf  <- dJSU(x3, mu = 9.137512, sigma = 2.781038, nu = 0.7902295, tau = 1.519528)
am_meanonly_JSU_1975_pdf  <- dJSU(x3, mu = 6.884389, sigma = 2.781038, nu = 0.7902295, tau = 1.519528)
am_meanonly_JSU_1992_pdf  <- dJSU(x3, mu = 3.535727, sigma = 2.781038, nu = 0.7902295, tau = 1.519528)



am_all_JSU_1992_pdf  <- dJSU(x2, mu = 3.297202, sigma = 2.057176, nu = 3.311038, tau = 1.999770)
am_meanonly_JSU_1992_pdf  <- dJSU(x3, mu = 3.535727, sigma = 2.781038, nu = 0.7902295, tau = 1.519528)



# plot density plots for 1964 ------------------------------
plot(Ammonium_1964_pdf, xlim= c(0,25), lwd = 2, lty = 2, col = "black")
#lines(am_all_SSTtr_1964_pdf, type = "l", lwd = 2, col = "goldenrod") # natural scale
#lines(am_meanonly_SSTtr_1964_pdf, type = "l", lwd = 2, col = "chocolate") # natural scale 
lines(am_all_JSU_1964_pdf, type = "l", lwd = 2, col = "darkred") # natural scale
lines(am_meanonly_JSU_1964_pdf, type = "l", lwd = 2, col = "pink") # natural scale
legend("topright", legend = c("empirical", "all; SSTtr", "mean only; SSTtr", 
                              "all; JSU", "mean only; JSU"),
       col = c("black", "goldenrod", "chocolate",
               "darkred", "pink"), cex = 0.7, lty = c(2,1,1,1), bty = "n")

mean(Ammonium_1964$Ammonium)
sd(Ammonium_1964$Ammonium)
skewness(Ammonium_1964$Ammonium)
kurtosis(Ammonium_1964$Ammonium)

# plot density plots for 1975 ------------------------------
plot(Ammonium_1975_pdf, xlim= c(0,25), ylim= c(0,0.19),lwd = 2, lty = 2, col = "black")
#lines(am_all_SSTtr_1975_pdf, type = "l", lwd = 2, col = "goldenrod") # natural scale
#lines(am_meanonly_SSTtr_1975_pdf, type = "l", lwd = 2, col = "chocolate") # natural scale 
lines(am_all_JSU_1975_pdf, type = "l", lwd = 2, col = "darkred") # natural scale
lines(am_meanonly_JSU_1975_pdf, type = "l", lwd = 2, col = "pink") # natural scale
legend("topright", legend = c("empirical", "all; SSTtr", "mean only; SSTtr", 
                              "all; JSU", "mean only; JSU"),
       col = c("black", "goldenrod", "chocolate",
               "darkred", "pink"), cex = 0.7, lty = c(2,1,1,1), bty = "n")


mean(Ammonium_1975$Ammonium)
sd(Ammonium_1975$Ammonium)
skewness(Ammonium_1975$Ammonium)
kurtosis(Ammonium_1975$Ammonium)


# plot density plots for 1992 ------------------------------
plot(Ammonium_1992_pdf, lwd = 2, lty = 2, col = "black")
#lines(am_all_SSTtr_1992_pdf, type = "l", lwd = 2, col = "goldenrod") # natural scale
#lines(am_meanonly_SSTtr_1992_pdf, type = "l", lwd = 2, col = "chocolate") # natural scale 
lines(am_all_JSU_1992_pdf, type = "l", lwd = 2, col = "darkred") # natural scale
lines(am_meanonly_JSU_1992_pdf, type = "l", lwd = 2, col = "pink") # natural scale
legend("topright", legend = c("empirical", "all; SSTtr", "mean only; SSTtr", 
                              "all; JSU", "mean only; JSU"),
       col = c("black", "goldenrod", "chocolate",
               "darkred", "pink"), cex = 0.7, lty = c(2,1,1,1), bty = "n")

mean(Ammonium_1992$Ammonium)
sd(Ammonium_1992$Ammonium)
skewness(Ammonium_1992$Ammonium)
kurtosis(Ammonium_1992$Ammonium)





#### Nitrite ###########################################################################################

# only compare fam distribution that has four parameters and when the mu and mean


# GT
nii_all_GT <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = GT(), data = Nitrite,
                       #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# JSU
nii_all_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = JSU(), data = Nitrite,
                       #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 


# SST
nii_all_SST <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = SST(), data = Nitrite,
                       #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 


AIC(nii_all_GT)
AIC(nii_all_JSU) # tis the best



# now lets check if all param is the best -----------------------------------------

## mean only 
nii_meanonly_JSU <- gamlss(Nitrite ~ year + month,
                      family = JSU(), data = Nitrite,
                      #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# mean and sigma 
nii_mean_and_sigma_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month,
                      family = JSU(), data = Nitrite,
                      #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# no skew 
nii_noskew_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                      family = JSU(), data = Nitrite,
                      mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# no kurt
nii_nokurt_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                      family = JSU(), data = Nitrite,
                      #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# all param
nii_all_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = JSU(), data = Nitrite,
                      #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 

# check AIC 
AIC(nii_meanonly_JSU)
AIC(nii_mean_and_sigma_JSU)
AIC(nii_noskew_JSU) # this has a better AIC (interesting)
AIC(nii_nokurt_JSU)
AIC(nii_all_JSU)


# PIT
# fitted param for (am_all_SSTtr) -----------
nii_mu_hat_all_JSU    <- predict(nii_noskew_JSU, "mu", type = "response")
nii_sigma_hat_all_JSU <- predict(nii_noskew_JSU, "sigma", type = "response")
nii_nu_hat_all_JSU    <- predict(nii_noskew_JSU, "nu", type = "response")
nii_tau_hat_all_JSU   <- predict(nii_noskew_JSU, "tau", type = "response")

am_pit8 <- pJSU(Nitrite$Nitrite, mu = nii_mu_hat_all_JSU, sigma = nii_sigma_hat_all_JSU, 
                  nu = nii_nu_hat_all_JSU, tau = nii_tau_hat_all_JSU)

# Plot PIT histogram
hist(am_pit8, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit7)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity



# empirical estimate of moments
nii_moments_summary <- Nitrite %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Nitrite, na.rm = TRUE),
    var = var(Nitrite, na.rm = TRUE),
    skew = skewness(Nitrite, na.rm = TRUE),
    kurt = kurtosis(Nitrite, na.rm = TRUE)
  )
print(nii_moments_summary, n=33)


plot(nii_moments_summary$year, nii_moments_summary$skew, type = "l")


Nitrite$mu_hat_JSU_best <- predict(nii_noskew_JSU, what = "mu", type = "response")
Nitrite$mu_hat_JSU_2best <- predict(nii_all_JSU, what = "mu", type = "response")


# plot that sucka
ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "grey", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_all_SSTtr), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_meanonly_SSTtr), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = mu_hat_JSU_best), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = mu_hat_JSU_2best), color = "steelblue", linewidth = 1) +
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()


resid_wp(nii_all_JSU)
resid_wp(nii_noskew_JSU)
resid_wp(nii_nokurt_JSU)



# lemme try and see if skewness is only for year 
# no kurt 
nii_nokurt2_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year,
                         family = JSU(), data = Nitrite,
                         mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
AIC(nii_nokurt2_JSU) # better! 

# all with skew only for year
nii_all2_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year, tau.fo = ~ year + month,
                      family = JSU(), data = Nitrite,
                      mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite), nu.start = 2, tau.start = 2,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
AIC(nii_all2_JSU)






#### Nitrate ##############################################################################

# JSU all param
ni_all_JSU <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year, tau.fo = ~ year + month,
                       family = JSU(), data = Nitrate,
                       #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 

# mean only
ni_meanonly_JSU <- gamlss(Nitrate ~ year + month,
                     family = JSU(), data = Nitrate,
                     #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 

# mean and sigma
ni_mean_and_sigma_JSU <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month,
                          family = JSU(), data = Nitrate,
                          #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 

# no skew
ni_noskew_JSU <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                                family = JSU(), data = Nitrate,
                                mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                                method = mixed(10,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 

# no kurt
ni_nokurt_JSU <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                                family = JSU(), data = Nitrate,
                                #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                                method = mixed(10,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 


# check AIC 
AIC(ni_meanonly_JSU)
AIC(ni_mean_and_sigma_JSU)
AIC(ni_noskew_JSU) # convergence issue 
AIC(ni_nokurt_JSU)
AIC(ni_all_JSU)


Nitrate$mu_hat_JSU_best <- predict(ni_all_JSU, what = "mu", type = "response")
Nitrate$mu_hat_SHASHo <- predict(ni_all_SHASHo, what = "mu", type = "response")
Nitrate$mu_hat_GAF <- predict(ni_all_GAF, what = "mu", type = "response")

Nitrate$JSU_sigma <- predict(ni_all_JSU, what = "sigma", type = "response")

ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "grey", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_all_SSTtr), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = am_mu_hat_meanonly_SSTtr), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = mu_hat_JSU_best), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = mu_hat_SHASHo), color = "steelblue", linewidth = 1) +
  geom_line(aes(y = mu_hat_GAF), color = "steelblue", linewidth = 1) +
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()



ni_all_SHASHo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year, tau.fo = ~ year + month,
                     family = SHASHo(), data = Nitrate,
                     mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 



## lemme test BCTo

ni_all_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = BCTo(), data = Nitrate,
                        mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)) 




# GB2 ---------------------------------------------------------------------------------
# mean only 
ni_meanonly_GB2 <- gamlss(Nitrate ~ pb(year) + month, 
                     family = GB2(), data = Nitrate,
                     mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE)).
#Warning message:
  #In regularize.values(x, y, ties, missing(ties)) :
  #collapsing to unique 'x' values


# mean and sigma 
ni_mean_and_sigma_GB2 <- gamlss(Nitrate ~ pb(year) + month, sigma.fo = ~ year + month, 
                          family = GB2(), data = Nitrate,
                          mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

# no skew
ni_noskew_GB2 <- gamlss(Nitrate ~ pb(year) + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                                family = GB2(), data = Nitrate,
                                mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                                method = mixed(10,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

# no kurt
ni_noskew_GB2 <- gamlss(Nitrate ~ pb(year) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                        family = GB2(), data = Nitrate,
                        mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

# all param
ni_all_GB2 <- gamlss(Nitrate ~ pb(year) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = GB2(), data = Nitrate,
                     mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate), nu.start = 2, tau.start = 2,
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))



                       