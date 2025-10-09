## READ ME ==============================================================================

# this script is for RQ2 on how moments are changing through time 
# starting with Ammonium, then Silicate, Phosphate, Nitrite, DIN, Nitrate

# 1. Run four-parameter family distribution that is on a positive real line distribution (BCT, BCTo, BCPE, BCPEo and GB2)
# 2. Run models with JSU or SST family distribution 
# 3. check which has the lowest AIC and use that family distribution 
# 4. Make sure that family can run all time-varying parameter models
# 5. IF NOT: try the second lowest AIC fam, so on and so forth 
# 5b. IF YES: check AIC for time varying parameters. 
# 6. Then check how the model predict the mean and plot them to see if it make sense 
# 7. IF NOT: try the second lowest AIC fam, so on and so forth
# 7b. IF YES: then that's your model baby!
# 8. EXTRA: try adding poly() to each parameter to the model and see how that does


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

# subset Silicate 
Silicate <- df %>% filter(!is.na(Silicate)) %>% 
  select(Date, Nitrite, year, month)

# subset Phosphate 
Phosphate <- df %>% filter(!is.na(Phosphate)) %>% 
  select(Date, Phosphate, year, month)

# subset Nitrite 
Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite, year, month)

# subset Nitrite 
DIN <- df %>% filter(!is.na(DIN)) %>% 
  select(Date, DIN, year, month)

# subset Nitrate 
Nitrate <- df %>% filter(!is.na(Nitrate)) %>% 
  select(Date, Nitrite, year, month)


## Ammonium ==============================================================================

# 1. Run four-parameter family distribution that is on a positive real line distribution --------------------------------------
#(BCT, BCTo, BCPE, BCPEo and GB2). run all param (time) 

#== BCT ==#
am_all_BCT <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCT(), data = Ammonium,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: right, BCT can't really run if values are close to 0, because the model turns mu negative during iteration, and 
## since BCT is on a positive real line distribution, and error pops up. 

#== BCTo ==#
am_all_BCTo <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCTo(), data = Ammonium,
                     mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got an error saying while (abs(G.dev.in - i.G.dev) > i.c.crit && i.iter < i.n.cyc) { : 
# missing value where TRUE/FALSE needed

#== BCPE ==#
am_all_BCPE <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCPE(), data = Ammonium,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== BCPEo ==#
am_all_BCPEo <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPEo(), data = Ammonium,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== GB2 ==#
am_all_GB2 <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = GB2(), data = Ammonium,
                       mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                       method = mixed(5,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
# all good here! 


# 2. Run models with JSU or SST family distribution ----------------------------------------------------------------------------

#== SST ==#
am_all_SST <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SST(), data = Ammonium,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: sigma must be positive error. always get this error when running SST on some nutrient! not sure why


#== JSU ==#
am_all_JSU <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = JSU(), data = Ammonium,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here!


# 3. check which has the lowest AIC and use that family distribution -----------------------------------------------------------

AIC(am_all_BCPE)
AIC(am_all_BCPEo)
AIC(am_all_GB2)
AIC(am_all_JSU)


# the rankings are:
## 1. BCPE
## 2. GB2
## 3. JSU
## 4. BCPEo


# 4. Make sure that family can run all time-varying parameter models ------------------------------------------------------------

## MEAN ONLY MODEL ##
am_meanonly_BCPE <- gamlss(Ammonium ~ year + month,
                      family = BCPE(), data = Ammonium,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
am_mean_and_sigma_BCPE <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month,
                           family = BCPE(), data = Ammonium,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
am_contskew_BCPE <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                                 family = BCPE(), data = Ammonium,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
am_contkurt_BCPE <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = BCPE(), data = Ammonium,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# all works. HOORAY! 


# 5b. IF YES: check AIC for time varying parameters -------------------------------------------------------------------------------

AIC(am_meanonly_BCPE)
AIC(am_mean_and_sigma_BCPE)
AIC(am_contskew_BCPE)
AIC(am_contkurt_BCPE)
AIC(am_all_BCPE)

# rankings are:
## 1. all param (time): 29124.04
## 2. constant kurtosis: 29156.9
## 3. constant skewness: 29369.81
## 4. mean and sigma: 29397.9
## 5. mean only: 29746.79


# 6. Then check how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Ammonium$mu_hat_all_BCPE <- predict(am_all_BCPE, what = "mu", type = "response")

ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPE), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()


# 8. EXTRA: try adding poly() to each parameter to the model and see how that does

##== poly to the mean ==## 
am_all_BCPE_meanpoly <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = BCPE(), data = Ammonium,
                       method = mixed(5,100),
                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: the AIC is lower then when I don't add poly...interesting: 29131.7


##== poly to the sigma ==## 
am_all_BCPE_sigmapoly <- gamlss(Ammonium ~ year + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCPE(), data = Ammonium,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: does very well: 29111.99


##== poly to the skewness ==## 
am_all_BCPE_skewpoly <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ year + month,
                                family = BCPE(), data = Ammonium,
                                method = mixed(5,100),
                                control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: the AIC is lower then when I don't add poly...interesting: 29134.41


##== poly to the kurtosis ==## 
am_all_BCPE_kurtpoly <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ poly(year,2) + month,
                               family = BCPE(), data = Ammonium,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: Does not do better: 29124.84

#######################################################
## CONCLUSION: adding polynomial to sigma is better! ##
#######################################################


# 6. RECHECK how the model predicts the mean and plot them to see if it make sense ----------------------------------------------


Ammonium$mu_hat_all_BCPE_sigmapoly <- predict(am_all_BCPE_sigmapoly, what = "mu", type = "response")

ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPE_sigmapoly), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()

## Looks the same COOL!


# 7b. IF YES: then that's your model baby! ---------------------------------------------------------------------------------------

# now let's check the summary of the model
summary(am_all_BCPE_sigmapoly)




