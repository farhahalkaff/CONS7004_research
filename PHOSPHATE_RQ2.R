## READ ME ==============================================================================

# this script is for RQ2 on how moments are changing through time 

## STEP 1 ##
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
# 9. Check PIT histogram and worm plot of the model

## STEP 2 ##
# 1. subset few years together to get period
# 2. get the estimated parameters for each period from the model
# 3. get non-parametric estimate of pdf 
# 4. Get PDF of those three year combination
# 5. plot that sucka to compare model prediction and non-parametric estimate 



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
library(ggridges)
library(stringr)

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
  select(Date, Silicate, year, month)

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
  select(Date, Nitrate, year, month)


## STEP 1 ==============================================================================================================================

# 1. Run four-parameter family distribution that is on a positive real line distribution --------------------------------------
#(BCT, BCTo, BCPE, BCPEo and GB2). run all param (time) 

#== BCT ==#
ph_all_BCT <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCT(), data = Phosphate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: right, BCT can't really run if values are close to 0, because the model turns mu negative during iteration, and 
## since BCT is on a positive real line distribution, and error pops up. 

#== BCTo ==#
ph_all_BCTo <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCTo(), data = Phosphate,
                      mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## all good here!

#== BCPE ==#
ph_all_BCPE <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPE(), data = Phosphate,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== BCPEo ==#
ph_all_BCPEo <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = BCPEo(), data = Phosphate,
                       method = mixed(5,100),
                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== GB2 ==#
ph_all_GB2 <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = GB2(), data = Phosphate,
                     mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
# all good here! 


# 2. Run models with JSU or SST family distribution ----------------------------------------------------------------------------

#== SST ==#
ph_all_SST <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SST(), data = Phosphate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: sigma must be positive error. always get this error when running SST on some nutrient! not sure why


#== JSU ==#
ph_all_JSU <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = JSU(), data = Phosphate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here!


# 3. check which has the lowest AIC and use that family distribution -----------------------------------------------------------

AIC(ph_all_BCTo)
AIC(ph_all_BCPE)
AIC(ph_all_BCPEo)
AIC(ph_all_GB2)
AIC(ph_all_JSU)


# the rankings are:
## 1. GB2: -507.8101
## 2. BCTo: 916.4988
## 3. BCPE: 991.1507
## 4. BCPEo: 1001.417
## 5. JSU: 1048.439


# 4. Make sure that family can run all time-varying parameter models ------------------------------------------------------------

## MEAN ONLY MODEL ##
ph_meanonly_GB2 <- gamlss(Phosphate ~ year + month,
                           family = GB2(), data = Phosphate,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
ph_mean_and_sigma_GB2 <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month,
                                 family = GB2(), data = Phosphate,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
ph_contskew_GB2 <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                           family = GB2(), data = Phosphate,
                          mu.start = mean(Phosphate$Phosphate), sigma.start  = sd(Phosphate$Phosphate),
                           method = mixed(5,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
ph_contkurt_GB2 <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = GB2(), data = Phosphate,
                          mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# all works. HOORAY! 


# 5b. IF YES: check AIC for time varying parameters -------------------------------------------------------------------------------

AIC(ph_meanonly_GB2)
AIC(ph_mean_and_sigma_GB2)
AIC(ph_contskew_GB2)
AIC(ph_contkurt_GB2)
AIC(ph_all_GB2)

# rankings are:
## 1. all param (time): -507.8101
## 2. constant skewness: 1043.067
## 3. constant kurtosis: 1059.469
## 4. mean and sigma: 1226.771
## 5. mean only: 2286.279


# 6. Then check how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Phosphate$mu_hat_all_GB2_ph <- predict(ph_all_GB2, what = "mu", type = "response")

ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_GB2_ph), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()


## NOTE: absolutely the fuck not! 


# 7. Let's try the second best AIC and see -----------------------------------------------------------------------------------------------------

### ==BCTo ==###

## MEAN ONLY MODEL ##
ph_meanonly_BCTo <- gamlss(Phosphate ~ year + month,
                          family = BCTo(), data = Phosphate,
                          method = mixed(5,100),
                          control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
ph_mean_and_sigma_BCTo <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month,
                                family = BCTo(), data = Phosphate,
                                method = mixed(5,100),
                                control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
ph_contskew_BCTo <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                          family = BCTo(), data = Phosphate,
                          mu.start = mean(Phosphate$Phosphate), sigma.start  = sd(Phosphate$Phosphate),
                          method = mixed(5,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
ph_contkurt_BCTo <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                          family = BCTo(), data = Phosphate,
                          mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                          method = mixed(5,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))


# now let's check AIC

AIC(ph_meanonly_BCTo)
AIC(ph_mean_and_sigma_BCTo)
AIC(ph_contskew_BCTo)
AIC(ph_contkurt_BCTo)
AIC(ph_all_BCTo)

# the rankings are 
## 1. All param (time): 916.4988
## 2. Constant kurtosis: 1003.649
## 3. Constant skewness: 1017.101
## 4. Mean and sigma: 1129.461
## 5. Mean only: 2208.809


# 6. Then check how AGAIN the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Phosphate$mu_hat_all_BCTo_ph <- predict(ph_all_BCTo, what = "mu", type = "response")

ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCTo_ph), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()



# 8. EXTRA: try adding poly() to each parameter to the model and see how that does -------------------------------------------------------------

##== poly to the mean ==## 
ph_all_BCTo_meanpoly <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCTo(), data = Phosphate,
                               mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: got the error: Error in while (abs(G.dev.in - i.G.dev) > i.c.crit && i.iter < i.n.cyc) { : 
# missing value where TRUE/FALSE needed

## NOTE: the mean dos look like it goes up and down, so might want to change family distirbution again


# 9. Do the third best family distribution  -------------------------------------------------------------

# first check whether we can add polynomial to the mean first
ph_all_BCPE_meanpoly <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCPE(), data = Phosphate,
                               mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: Ok works awesome! 

## Now let's make sure that all param (time) has the lowest AIC 

### ==BCPE ==###

## MEAN ONLY MODEL ##
ph_meanonly_BCPE <- gamlss(Phosphate ~ year + month,
                           family = BCPE(), data = Phosphate,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
ph_mean_and_sigma_BCPE <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month,
                                 family = BCPE(), data = Phosphate,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
ph_contskew_BCPE <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                           family = BCPE(), data = Phosphate,
                           mu.start = mean(Phosphate$Phosphate), sigma.start  = sd(Phosphate$Phosphate),
                           method = mixed(5,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
ph_contkurt_BCPE <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = BCPE(), data = Phosphate,
                           mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                           method = mixed(5,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))


# now let's check AIC

AIC(ph_meanonly_BCPE)
AIC(ph_mean_and_sigma_BCPE)
AIC(ph_contskew_BCPE)
AIC(ph_contkurt_BCPE)
AIC(ph_all_BCPE)

# the rankings are 
## 1. All param (time): 991.1507
## 2. Constant kurtosis: 1002.119
## 3. Constant skewness: 1110.707
## 4. Mean and sigma: 1159.482
## 5. Mean only: 2248.947


# 6. RECHECK how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Phosphate$mu_hat_all_BCPE_ph <- predict(ph_all_BCPE, what = "mu", type = "response")

ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPE_ph), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()



# 8. EXTRA: try adding poly() to each parameter to the model and see how that does -------------------------------------------------------------

##== poly to the mean ==## 
ph_all_BCPE_meanpoly <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCPE(), data = Phosphate,
                               mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: PRETTY GOOD: 168.2585


##== poly to the mean and sigma ==## 
ph_all_BCPE_sigmapoly <- gamlss(Phosphate ~ year + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCPE(), data = Phosphate,
                               mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: yeah does pretty well! : 835.042

##== poly to the skewness ==## 
ph_all_BCPE_skewpoly <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ year + month,
                                family = BCPE(), data = Phosphate,
                                mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                method = mixed(5,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: yeah did alright: 959.0701

##== poly to the kurtosis ==## 
ph_all_BCPE_kurtpoly <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ poly(year,2) + month,
                                family = BCPE(), data = Phosphate,
                                mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                method = mixed(5,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: nah not well: 992.9104


### == combine mean sigma and skew == ###
ph_all_BCPE_allpoly <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ year + month,
                               family = BCPE(), data = Phosphate,
                               mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: very low AIC: 4.46095, but gives me negative global deviance, so nah let's not use this 


# 6. RECHECK how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Phosphate$mu_hat_all_meanpoly_BCPE_ph <- predict(ph_all_BCPE_meanpoly, what = "mu", type = "response")

ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_meanpoly_BCPE_ph), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()

## Looks COOL!



# 9. Check PIT histogram, worm plot and moment bucket of the model -----------------------------------------------------------------------------

# fitted param for (ph_all_BCPE_meanpoly) -----------
ph_mu_hat_all_BCPE_meanpoly    <- predict(ph_all_BCPE_meanpoly, "mu", type = "response")
ph_sigma_hat_all_BCPE_meanpoly <- predict(ph_all_BCPE_meanpoly, "sigma", type = "response")
ph_nu_hat_all_BCPE_meanpoly    <- predict(ph_all_BCPE_meanpoly, "nu", type = "response")
ph_tau_hat_all_BCPE_meanpoly   <- predict(ph_all_BCPE_meanpoly, "tau", type = "response")

ph_pit <- pBCPE(Phosphate$Phosphate, mu = ph_mu_hat_all_BCPE_meanpoly, sigma = ph_sigma_hat_all_BCPE_meanpoly, 
                nu = ph_nu_hat_all_BCPE_meanpoly, tau = ph_tau_hat_all_BCPE_meanpoly)

# Plot PIT histogram
hist(ph_pit, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(ph_pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity
## NOTE: yeah looks pretty good

resid_wp(ph_all_BCPE_meanpoly)
## NOTE: hmm tails a bit off

moment_bucket(ph_all_BCPE_meanpoly)
## NOTE: a bit of heavy tailed residuals


# 7b. IF YES: then that's your model baby! ---------------------------------------------------------------------------------------

# now let's check the summary of the model
summary(am_all_BCPE_sigmapoly)

## ==================================================================================================================================


## STEP 2 ==============================================================================================================================

# 1. subset few years together to get period ---------------------------------------------------------------------------------------
# Add column to combine few years together 
Phosphate <- Phosphate %>% 
  mutate(period = case_when(
    between(year, 1962, 1966) ~ "1",
    between(year, 1967, 1970) ~ "2",
    between(year, 1971, 1974) ~ "3",
    between(year, 1975, 1978) ~ "4",
    between(year, 1979, 1982) ~ "5",
    between(year, 1983, 1986) ~ "6",
    between(year, 1987, 1990) ~ "7", 
    between(year, 1991, 1994) ~ "8"
  ))


# 2. get the estimated parameters for each period ---------------------------------------------------------------------------------

# create function to get predicted param from model ======
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- predict(model, what = "mu", type = "response")
  data$sigma_hat <- predict(model, what = "sigma", type = "response")
  data$nu_hat <- predict(model, what = "nu", type = "response")
  data$tau_hat <- predict(model, what = "tau", type = "response")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(period) %>%
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

# get the parameters for the year periods ======
ph_all_BCPE_meanpoly_param <- best_model_param(ph_all_BCPE_meanpoly, Phosphate)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
Phosphate_A <- Phosphate[Phosphate$period %in% c("1"), ]
Phosphate_A_pdf <- density(Phosphate_A$Phosphate)

# Period B
Phosphate_B <- Phosphate[Phosphate$period %in% c("2"), ]
Phosphate_B_pdf <- density(Phosphate_B$Phosphate)

# Period C
Phosphate_C <- Phosphate[Phosphate$period %in% c("3"), ]
Phosphate_C_pdf <- density(Phosphate_C$Phosphate)

# Period D
Phosphate_D <- Phosphate[Phosphate$period %in% c("4"), ]
Phosphate_D_pdf <- density(Phosphate_D$Phosphate)

# Period E
Phosphate_E <- Phosphate[Phosphate$period %in% c("5"), ]
Phosphate_E_pdf <- density(Phosphate_E$Phosphate)

# Period F
Phosphate_F <- Phosphate[Phosphate$period %in% c("6"), ]
Phosphate_F_pdf <- density(Phosphate_F$Phosphate)

# Period G
Phosphate_G <- Phosphate[Phosphate$period %in% c("7"), ]
Phosphate_G_pdf <- density(Phosphate_G$Phosphate)

# Period H
Phosphate_H <- Phosphate[Phosphate$period %in% c("8"), ]
Phosphate_H_pdf <- density(Phosphate_H$Phosphate)


# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(ph_all_BCPE_meanpoly$y, na.rm = TRUE),
         max(ph_all_BCPE_meanpoly$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
ph_all_BCPE_meanpoly_A_pdf  <- dBCPE(x, mu = 0.6049097, sigma = 0.4088816, nu = 0.4300382, tau = 1.844037)
ph_all_BCPE_meanpoly_B_pdf  <- dBCPE(x, mu = 0.6988106, sigma = 0.4243612, nu = 0.4307906, tau = 1.788255)
ph_all_BCPE_meanpoly_C_pdf  <- dBCPE(x, mu = 0.7750707, sigma = 0.4233480, nu = 0.4434835, tau = 1.715349)
ph_all_BCPE_meanpoly_D_pdf  <- dBCPE(x, mu = 0.8004955, sigma = 0.4299985, nu = 0.4533161, tau = 1.649334)
ph_all_BCPE_meanpoly_E_pdf  <- dBCPE(x, mu = 0.8021530, sigma = 0.4338234, nu = 0.4756160, tau = 1.596055)
ph_all_BCPE_meanpoly_F_pdf  <- dBCPE(x, mu = 0.7675777, sigma = 0.4394177, nu = 0.4880571, tau = 1.541047)
ph_all_BCPE_meanpoly_G_pdf  <- dBCPE(x, mu = 0.7047156, sigma = 0.4418867, nu = 0.4989334, tau = 1.482365)
ph_all_BCPE_meanpoly_H_pdf  <- dBCPE(x, mu = 0.5971119, sigma = 0.4509452, nu = 0.5067932, tau = 1.429937)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(Phosphate_A_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_meanpoly_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(Phosphate_B_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_meanpoly_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(Phosphate_C_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.13))
lines(ph_all_BCPE_meanpoly_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(Phosphate_D_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_meanpoly_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(Phosphate_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_meanpoly_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(Phosphate_F_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_meanpoly_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(Phosphate_G_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_meanpoly_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(Phosphate_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.3))
lines(ph_all_BCPE_meanpoly_H_pdf, type = "l", lwd = 2, col = "goldenrod") 


## NOTE: It seems that whenever I add polynomial, the pdf does pretty shitly, lemme check how it does
## with not polynomial 


# get the parameters for the year periods ======
ph_all_BCPE_param <- best_model_param(ph_all_BCPE, Phosphate)


# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(ph_all_BCPE$y, na.rm = TRUE),
         max(ph_all_BCPE$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
ph_all_BCPE_A_pdf  <- dBCPE(x, mu = 0.7700470, sigma = 0.4560479, nu = 0.4365239, tau = 2.243866)
ph_all_BCPE_B_pdf  <- dBCPE(x, mu = 0.7411489, sigma = 0.4601070, nu = 0.4460963, tau = 2.149845)
ph_all_BCPE_C_pdf  <- dBCPE(x, mu = 0.7432048, sigma = 0.4482091, nu = 0.4639943, tau = 2.052299)
ph_all_BCPE_D_pdf  <- dBCPE(x, mu = 0.7298875, sigma = 0.4469063, nu = 0.4777252, tau = 1.954415)
ph_all_BCPE_E_pdf  <- dBCPE(x, mu = 0.7244989, sigma = 0.4406815, nu = 0.5006972, tau = 1.874777)
ph_all_BCPE_F_pdf  <- dBCPE(x, mu = 0.7154570, sigma = 0.4365310, nu = 0.5174929, tau = 1.795690)
ph_all_BCPE_G_pdf  <- dBCPE(x, mu = 0.7116077, sigma = 0.4292786, nu = 0.5321648, tau = 1.716344)
ph_all_BCPE_H_pdf  <- dBCPE(x, mu = 0.6949121, sigma = 0.4280224, nu = 0.5455526, tau = 1.639262)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(Phosphate_A_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(Phosphate_B_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(Phosphate_C_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.13))
lines(ph_all_BCPE_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(Phosphate_D_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(Phosphate_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(Phosphate_F_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(Phosphate_G_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ph_all_BCPE_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(Phosphate_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.3))
lines(ph_all_BCPE_H_pdf, type = "l", lwd = 2, col = "goldenrod") 


## NOTE: does the same unfortunately 

ph_moments_summary <- Phosphate %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Phosphate, na.rm = TRUE),
    var = var(Phosphate, na.rm = TRUE),
    skew = skewness(Phosphate, na.rm = TRUE),
    kurt = kurtosis(Phosphate, na.rm = TRUE)
  )
print(ph_moments_summary, n=33)


################

# METHOD 1 ===========
logSurv(Phosphate$Phosphate, prob=0.90, tail="right")

# METHOD 2 ==========
par(mfrow = c(1, 3))
ph_m1 <- loglogSurv1(Phosphate$Phosphate, prob=0.90, title="(a) TYPE I")
ph_m2 <- loglogSurv2(Phosphate$Phosphate, prob=0.90, title="(b) TYPE II")
ph_m3 <- loglogSurv3(Phosphate$Phosphate, prob=0.90, title="(c) TYPE III")
par(mfrow = c(1, 1))

# METHOD 3 ===========

# use some truncated family distribution 
ph_m4 <- fitTail(Phosphate$Phosphate, family=WEI, percentage=10)
ph_m5 <- fitTail(Phosphate$Phosphate, family=LOGNO, percentage=10)
ph_m6 <- fitTail(Phosphate$Phosphate, family=BCPE, percentage=10)
ph_m7 <- fitTail(Phosphate$Phosphate, family=BCTo, percentage=10) # got warning 
ph_m8 <- fitTail(Phosphate$Phosphate, family=SEP4, percentage=10) # got warning 


# check AIC 
AIC(ph_m4, ph_m5, ph_m6, ph_m7, ph_m8)
## NOTE: BCPE does the best

wp(ph_m4, ylim.all = 1)
am_m6_2 <- fitTailAll(Phosphate$Phosphate, family=BCPE) # got like 50 warnings
plot(am_m6_2)


# METHOD 4 ===========
ph_f1 <- fitDist(Phosphate$Phosphate, ncpus=4, parallel="snow" )
ph_f1$fits[1:25]
wp(ph_f1, ylim.all=1.5)

# now log them
am_z<-log(Ammonium$Ammonium)
am_f2 <- fitDist(am_z, ncpus=4, parallel="snow")
am_f2$fits[1:5]

wp(am_f2, ylim.all=1)
am_f3 <- histDist(am_z, family="BCPE", nbins=30)

gen.Family("BCPE", "log")
am_f4 <- histDist(Ammonium$Ammonium, family="logBCPE", nbins=30)
AIC(am_f1, am_f4)





# 1. Total number of observations
total_obs <- length(Phosphate$Phosphate)
print(paste("Total observations:", total_obs))

# 2. Number of observations in the top 10%
num_top_10_percent <- ceiling(total_obs * 0.10)
print(paste("Number of observations in the top 10%:", num_top_10_percent))



# empirical estimate of moments
ph_moments_summary <- Phosphate %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Phosphate, na.rm = TRUE),
    var = var(Phosphate, na.rm = TRUE),
    skew = skewness(Phosphate, na.rm = TRUE),
    kurt = kurtosis(Phosphate, na.rm = TRUE)
  )
print(ph_moments_summary, n=33)








ph_all_BCPE_years_spline <- gamlss(Phosphate ~ cs(year) + month, sigma.fo = ~ cs(year) + month, nu.fo = ~ cs(year) + month, tau.fo = ~ cs(year) + month,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

Phosphate$mu_hat_years <- predict(ph_all_BCPE_years_spline, what ="tau", type = "response")

##### use bins in the model 
Phosphate$period <- as.numeric(Phosphate$period)

ph_all_BCPE_bins_spline <- gamlss(Phosphate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                   family = BCPE(), data = Phosphate,
                                   mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bins_spline_noseason <- gamlss(Phosphate ~ cs(period), sigma.fo = ~ cs(period), nu.fo = ~ cs(period), tau.fo = ~ cs(period),
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


ph_all_BCPE_bins_spline_param <- best_model_param(ph_all_BCPE_bins_spline, Phosphate)
ph_all_BCPE_bins_noseason_spline_param <- best_model_param(ph_all_BCPE_bins_spline_noseason, Phosphate)

## with seasonaility 
set.seed(123)
ph_all_BCPE_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5941064, sigma = 0.5179176, nu = 0.3315671, tau = 1.560181),
  "2" = rBCPE(532, mu = 0.6737867, sigma = 0.4540171, nu = 0.3044036, tau = 1.794850),
  "3" = rBCPE(562, mu = 0.7783245, sigma = 0.3705839, nu = 0.2936626, tau = 2.147731),
  "4" = rBCPE(862, mu = 0.8465924, sigma = 0.3437207, nu = 0.2624391, tau = 2.313994),
  "5" = rBCPE(959, mu = 0.8564020, sigma = 0.3603858, nu = 0.3093610, tau = 2.066227),
  "6" = rBCPE(980, mu = 0.7809474, sigma = 0.3857275, nu = 0.4149186, tau = 1.779417),
  "7" = rBCPE(973, mu = 0.6617136, sigma = 0.4187571, nu = 0.4527497, tau = 1.480159),
  "8" = rBCPE(968, mu = 0.5977861, sigma = 0.4953745, nu = 0.4091078, tau = 1.362330)
)


ph_all_BCPE_bins_spline_predict_long <- stack(ph_all_BCPE_bins_spline_predict)
names(ph_all_BCPE_bins_spline_predict_long) <- c("pred", "period")
ph_all_BCPE_bins_spline_predict_long$period <- as.numeric(ph_all_BCPE_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Phosphate, aes(x = Phosphate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = ph_all_BCPE_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  


## without seasonaility 
set.seed(123)
ph_all_BCPE_bins_noseason_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5223994, sigma = 0.5911594, nu = 0.6047738, tau = 2.499300),
  "2" = rBCPE(532, mu = 0.6402502, sigma = 0.5844982, nu = 0.8870768, tau = 2.310119),
  "3" = rBCPE(562, mu = 0.7352998, sigma = 0.5809116, nu = 1.1691256, tau = 2.651172),
  "4" = rBCPE(862, mu = 0.8010986, sigma = 0.6139708, nu = 1.3181606, tau = 3.007827),
  "5" = rBCPE(959, mu = 0.8050691, sigma = 0.6584863, nu = 1.3327172, tau = 2.852733),
  "6" = rBCPE(980, mu = 0.7365117, sigma = 0.6865062, nu = 1.2590187, tau = 2.631527),
  "7" = rBCPE(973, mu = 0.6169694, sigma = 0.6915215, nu = 1.1317877, tau = 2.314725),
  "8" = rBCPE(968, mu = 0.5352415, sigma = 0.7356121, nu = 0.9911742, tau = 2.034754)
)

ph_all_BCPE_bins_noseson_spline_predict_long <- stack(ph_all_BCPE_bins_noseason_spline_predict)
names(ph_all_BCPE_bins_noseson_spline_predict_long) <- c("pred", "period")
ph_all_BCPE_bins_noseson_spline_predict_long$period <- as.numeric(ph_all_BCPE_bins_noseson_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Phosphate, aes(x = Phosphate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = ph_all_BCPE_bins_noseson_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  



Phosphate$mu_hat_season <-  predict(ph_all_BCPE_bins_spline, what = "tau", type = "response")
Phosphate$mu_hat_noseason <-  predict(ph_all_BCPE_bins_spline_noseason, what = "tau", type = "response")


# method 1: use moments package (kurtosis() returns Pearson kurtosis; subtract 3 to get excess)
excess_kurtosis <- kurtosis(y_sim) - 3




### bins, changing paramter time-varying 

ph_all_BCPE_bins_spline <- gamlss(Phosphate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


ph_nokurt_BCPE_bins_spline <- gamlss(Phosphate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ 1,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_noskew_BCPE_bins_spline <- gamlss(Phosphate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ cs(period) + month,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_mean_sigma_BCPE_bins_spline <- gamlss(Phosphate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ 1,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


ph_meanonly_BCPE_bins_spline <- gamlss(Phosphate ~ cs(period) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


AIC(ph_all_BCPE_bins_spline,ph_nokurt_BCPE_bins_spline,
    ph_noskew_BCPE_bins_spline, ph_mean_sigma_BCPE_bins_spline,
    ph_meanonly_BCPE_bins_spline)


### RESULTS 
# ph_all_BCPE_bins_spline         -320.7370
# ph_nokurt_BCPE_bins_spline      -266.0001
# ph_noskew_BCPE_bins_spline      -212.3927
# ph_mean_sigma_BCPE_bins_spline  -118.4800
# ph_meanonly_BCPE_bins_spline     1366.2981



library(gamlss.tr)
gen.trun(family = "SEP4", type = "left")

ph_intercept_SEP4 <- gamlss(Phosphate ~  1,
                                       family = SEP4(), data = Phosphate,
                                       mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                       method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# truncate?
ph_intercept_SEP4 <- gamlss(Phosphate ~  1,
                            family = SEP4(), data = Phosphate,
                            mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                            method = mixed(5,100),
                            control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))




moment_bucket(ph_intercept_SEP4) +
resid_wp(ph_intercept_SEP4) + theme(plot.title = element_blank())









### model for each bin 

ph_all_BCPE_bin_1 <- gamlss(Phosphate ~ month, sigma.fo = ~  month, nu.fo = ~ month, tau.fo = ~ month,
                                  family = BCPE(), data = Phosphate_A,
                                  mu.start = mean(Phosphate_A$Phosphate), sigma.start = sd(Phosphate_A$Phosphate),
                                  #method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


ph_all_BCPE_bin_2 <- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~ month,
                                      family = BCPE(), data = Phosphate_B,
                                      mu.start = mean(Phosphate_B$Phosphate), sigma.start = sd(Phosphate_B$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bin_3 <- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~  month,
                                      family = BCPE(), data = Phosphate_C,
                                      mu.start = mean(Phosphate_C$Phosphate), sigma.start = sd(Phosphate_C$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bin_4 <- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~ month,
                                      family = BCPE(), data = Phosphate_D,
                                      mu.start = mean(Phosphate_D$Phosphate), sigma.start = sd(Phosphate_D$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bin_5 <- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~ month,
                                      family = BCPE(), data = Phosphate_E,
                                      mu.start = mean(Phosphate_E$Phosphate), sigma.start = sd(Phosphate_E$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bin_6 <- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~ month,
                                      family = BCPE(), data = Phosphate_F,
                                      mu.start = mean(Phosphate_F$Phosphate), sigma.start = sd(Phosphate_F$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bin_7<- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~  month,
                                      family = BCPE(), data = Phosphate_G,
                                      mu.start = mean(Phosphate_G$Phosphate), sigma.start = sd(Phosphate_G$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ph_all_BCPE_bin_8 <- gamlss(Phosphate ~ month, sigma.fo = ~ month, nu.fo = ~ month, tau.fo = ~  month,
                                      family = BCPE(), data = Phosphate_H,
                                      mu.start = mean(Phosphate_H$Phosphate), sigma.start = sd(Phosphate_H$Phosphate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


bins_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- predict(model, what = "mu", type = "response")
  data$sigma_hat <- predict(model, what = "sigma", type = "response")
  data$nu_hat <- predict(model, what = "nu", type = "response")
  data$tau_hat <- predict(model, what = "tau", type = "response")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(period) %>%
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
  #return(average_by_year)
  return(average_by_month)
  
}


ph_all_BCPE_bin_1_param <- best_model_param(ph_all_BCPE_bin_1, Phosphate_A)
ph_all_BCPE_bin_2_param <- best_model_param(ph_all_BCPE_bin_2, Phosphate_B)
ph_all_BCPE_bin_3_param <- best_model_param(ph_all_BCPE_bin_3, Phosphate_C)
ph_all_BCPE_bin_4_param <- best_model_param(ph_all_BCPE_bin_4, Phosphate_D)
ph_all_BCPE_bin_5_param <- best_model_param(ph_all_BCPE_bin_5, Phosphate_E)
ph_all_BCPE_bin_6_param <- best_model_param(ph_all_BCPE_bin_6, Phosphate_F)
ph_all_BCPE_bin_7_param <- best_model_param(ph_all_BCPE_bin_7, Phosphate_G)
ph_all_BCPE_bin_8_param <- best_model_param(ph_all_BCPE_bin_8, Phosphate_H)


separate_bin_models_param <- rbind(ph_all_BCPE_bin_1_param, ph_all_BCPE_bin_2_param, 
                                   ph_all_BCPE_bin_3_param, ph_all_BCPE_bin_4_param, 
                                   ph_all_BCPE_bin_5_param, ph_all_BCPE_bin_6_param, 
                                   ph_all_BCPE_bin_7_param, ph_all_BCPE_bin_8_param)


set.seed(123)
ph_all_BCPE_separate_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5282223, sigma = 0.4185144, nu = 0.23108453, tau = 2.170455),
  "2" = rBCPE(532, mu = 0.6315479, sigma = 0.4187895, nu = 0.27968359, tau = 1.999609),
  "3" = rBCPE(562, mu = 0.7701517, sigma = 0.3275054, nu = 0.67620624, tau = 90995.156949),
  "4" = rBCPE(862, mu = 0.8642670, sigma = 0.2925768, nu = 0.21534561, tau = 2.672778),
  "5" = rBCPE(959, mu = 0.8756509, sigma = 0.3518320, nu = 0.38041203, tau = 2.416395),
  "6" = rBCPE(980, mu = 0.8572213, sigma = 0.3376212, nu = 0.83626801, tau = 3.160718),
  "7" = rBCPE(973, mu = 0.6251784, sigma = 0.4027832, nu = 0.60655177, tau = 2.315820),
  "8" = rBCPE(968, mu = 0.5591880, sigma = 0.9247212, nu = 0.05446912, tau = 1.760222)
)

ph_all_BCPE_separate_bins_spline_predict_long <- stack(ph_all_BCPE_separate_bins_spline_predict)
names(ph_all_BCPE_separate_bins_spline_predict_long) <- c("pred", "period")
ph_all_BCPE_separate_bins_spline_predict_long$period <- as.numeric(ph_all_BCPE_separate_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Phosphate, aes(x = Phosphate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = ph_all_BCPE_separate_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  






### doing linear and quadratic 

# linear
ph_all_BCPE_bins_linear <- gamlss(Phosphate ~ period + month, sigma.fo = ~ period + month, nu.fo = ~ period + month, tau.fo = ~ period + month,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# quadratic
ph_all_BCPE_bins_quad <- gamlss(Phosphate ~ poly(period,2) + month, sigma.fo = ~ poly(period,2) + month, nu.fo = ~ poly(period,2) + month, tau.fo = ~ poly(period,2) + month,
                                  family = BCPE(), data = Phosphate,
                                  mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


ph_all_BCPE_bins_linear_param <- best_model_param(ph_all_BCPE_bins_linear, Phosphate)
ph_all_BCPE_bins_quad_param <- best_model_param(ph_all_BCPE_bins_quad, Phosphate)












# simulation for seasonality 

ph_all_BCPE_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5941064, sigma = 0.5179176, nu = 0.3315671, tau = 1.560181),
  "2" = rBCPE(532, mu = 0.6737867, sigma = 0.4540171, nu = 0.3044036, tau = 1.794850),
  "3" = rBCPE(562, mu = 0.7783245, sigma = 0.3705839, nu = 0.2936626, tau = 2.147731),
  "4" = rBCPE(862, mu = 0.8465924, sigma = 0.3437207, nu = 0.2624391, tau = 2.313994),
  "5" = rBCPE(959, mu = 0.8564020, sigma = 0.3603858, nu = 0.3093610, tau = 2.066227),
  "6" = rBCPE(980, mu = 0.7809474, sigma = 0.3857275, nu = 0.4149186, tau = 1.779417),
  "7" = rBCPE(973, mu = 0.6617136, sigma = 0.4187571, nu = 0.4527497, tau = 1.480159),
  "8" = rBCPE(968, mu = 0.5977861, sigma = 0.4953745, nu = 0.4091078, tau = 1.362330)
)


library(gamlss.dist)
library(moments)

############################## bin 1
mu_ph_bin1   <- 0.5941064    
sigma_ph_bin1 <- 0.5179176   
nu_ph_bin1  <- 0.3315671   
tau_ph_bin1  <- 1.560181

set.seed(1)
N <- 5e5                    
ph_bin1_sim <- rBCPE(N, mu = mu_ph_bin1, sigma = sigma_ph_bin1, nu = nu_ph_bin1, tau = tau_ph_bin1)

# get skewness and kurtosis 
kurtosis(ph_bin1_sim) - 3 # 3.716351
skewness(ph_bin1_sim) # 1.350677
var(ph_bin1_sim) # 0.1092073
mean(ph_bin1_sim) # 0.647022

############################# bin 2
mu_ph_bin2   <- 0.6737867    
sigma_ph_bin2 <- 0.4540171   
nu_ph_bin2  <- 0.3044036   
tau_ph_bin2  <- 1.794850

set.seed(2)
ph_bin2_sim <- rBCPE(N, mu = mu_ph_bin2, sigma = sigma_ph_bin2, nu = nu_ph_bin2, tau = tau_ph_bin2)

# get skewness and kurtosis 
kurtosis(ph_bin2_sim) - 3 # 2.716351
skewness(ph_bin2_sim) # 1.051625
var(ph_bin2_sim) # 0.1049501
mean(ph_bin2_sim) # 0.7228153

############################# bin 3
mu_ph_bin3   <- 0.7783245    
sigma_ph_bin3 <- 0.3705839  
nu_ph_bin3  <- 0.2936626   
tau_ph_bin3  <- 2.147731

set.seed(3)
ph_bin3_sim <- rBCPE(N, mu = mu_ph_bin3, sigma = sigma_ph_bin3, nu = nu_ph_bin3, tau = tau_ph_bin3)

# get skewness and kurtosis 
kurtosis(ph_bin3_sim) - 3 # 0.6636381
skewness(ph_bin3_sim) # 0.7217592
var(ph_bin3_sim) # 0.08891965
mean(ph_bin3_sim) # 0.8159112


############################# bin 4
mu_ph_bin4   <- 0.8465924    
sigma_ph_bin4 <- 0.3437207  
nu_ph_bin4  <- 0.2624391   
tau_ph_bin4  <- 2.313994

set.seed(4)
ph_bin4_sim <- rBCPE(N, mu = mu_ph_bin4, sigma = sigma_ph_bin4, nu = nu_ph_bin4, tau = tau_ph_bin4)

# get skewness and kurtosis 
kurtosis(ph_bin4_sim) - 3 # 0.4070751
skewness(ph_bin4_sim) # 0.6572835
var(ph_bin4_sim) # 0.09034172
mean(ph_bin4_sim) # 0.8833055


############################# bin 5
mu_ph_bin5   <- 0.8564020    
sigma_ph_bin5 <- 0.3603858  
nu_ph_bin5  <- 0.3093610   
tau_ph_bin5  <- 2.066227

set.seed(5)
ph_bin5_sim <- rBCPE(N, mu = mu_ph_bin5, sigma = sigma_ph_bin5, nu = nu_ph_bin5, tau = tau_ph_bin5)

# get skewness and kurtosis 
kurtosis(ph_bin5_sim) - 3 # 0.7509475
skewness(ph_bin5_sim) # 0.7230283
var(ph_bin5_sim) # 0.1019335
mean(ph_bin5_sim) # 0.8946723


############################# bin 6
mu_ph_bin6   <- 0.7809474    
sigma_ph_bin6 <- 0.3857275  
nu_ph_bin6  <- 0.4149186   
tau_ph_bin6  <- 1.779417

set.seed(6)
ph_bin6_sim <- rBCPE(N, mu = mu_ph_bin6, sigma = sigma_ph_bin6, nu = nu_ph_bin6, tau = tau_ph_bin6)

# get skewness and kurtosis 
kurtosis(ph_bin6_sim) - 3 # 1.201999
skewness(ph_bin6_sim) # 0.7695995
var(ph_bin6_sim) # 0.09512644
mean(ph_bin6_sim) # 0.8157169


############################# bin 7
mu_ph_bin7   <- 0.6617136    
sigma_ph_bin7 <- 0.4187571  
nu_ph_bin7  <- 0.4527497   
tau_ph_bin7  <- 1.480159

set.seed(7)
ph_bin7_sim <- rBCPE(N, mu = mu_ph_bin7, sigma = sigma_ph_bin7, nu = nu_ph_bin7, tau = tau_ph_bin7)

# get skewness and kurtosis 
kurtosis(ph_bin7_sim) - 3 # 2.115725
skewness(ph_bin7_sim) # 0.932626
var(ph_bin7_sim) # 0.08039081
mean(ph_bin7_sim) # 0.6936117


############################# bin 8
mu_ph_bin8   <- 0.5977861    
sigma_ph_bin8 <- 0.4953745  
nu_ph_bin8  <- 0.4091078   
tau_ph_bin8  <- 1.362330

set.seed(8)
ph_bin8_sim <- rBCPE(N, mu = mu_ph_bin8, sigma = sigma_ph_bin8, nu = nu_ph_bin8, tau = tau_ph_bin8)

# get skewness and kurtosis 
kurtosis(ph_bin8_sim) - 3 # 3.937416
skewness(ph_bin8_sim) # 1.314376
var(ph_bin8_sim) # 0.09672909
mean(ph_bin8_sim) # 0.6408134


#=====================================================================================================

# simulation without seasonaility 

ph_all_BCPE_bins_noseason_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5223994, sigma = 0.5911594, nu = 0.6047738, tau = 2.499300),
  "2" = rBCPE(532, mu = 0.6402502, sigma = 0.5844982, nu = 0.8870768, tau = 2.310119),
  "3" = rBCPE(562, mu = 0.7352998, sigma = 0.5809116, nu = 1.1691256, tau = 2.651172),
  "4" = rBCPE(862, mu = 0.8010986, sigma = 0.6139708, nu = 1.3181606, tau = 3.007827),
  "5" = rBCPE(959, mu = 0.8050691, sigma = 0.6584863, nu = 1.3327172, tau = 2.852733),
  "6" = rBCPE(980, mu = 0.7365117, sigma = 0.6865062, nu = 1.2590187, tau = 2.631527),
  "7" = rBCPE(973, mu = 0.6169694, sigma = 0.6915215, nu = 1.1317877, tau = 2.314725),
  "8" = rBCPE(968, mu = 0.5352415, sigma = 0.7356121, nu = 0.9911742, tau = 2.034754)
)


############################## bin 1
mu_ph_bin1_noseason   <-    0.5223994 
sigma_ph_bin1_noseason  <-    0.5911594
nu_ph_bin1_noseason   <-    0.6047738
tau_ph_bin1_noseason   <- 2.499300

set.seed(11)
N <- 5e5                    
ph_bin1_noseason_sim <- rBCPE(N, mu = mu_ph_bin1_noseason, sigma = sigma_ph_bin1_noseason, 
                               nu = nu_ph_bin1_noseason , tau = tau_ph_bin1_noseason )

# get skewness and kurtosis 
kurtosis(ph_bin1_noseason_sim) - 3 # -0.01946382
skewness(ph_bin1_noseason_sim) # 0.5787896
var(ph_bin1_noseason_sim) # 0.0949779
mean(ph_bin1_noseason_sim) # 0.5603139

############################# bin 2
mu_ph_bin2_noseason   <-     0.6402502
sigma_ph_bin2_noseason <-    0.5911594
nu_ph_bin2_noseason  <-    0.8870768
tau_ph_bin2_noseason  <-  2.310119

set.seed(22)
ph_bin2_noseason_sim <- rBCPE(N, mu = mu_ph_bin2_noseason, sigma = sigma_ph_bin2_noseason, 
                              nu = nu_ph_bin2_noseason, tau = tau_ph_bin2_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin2_noseason_sim) - 3 # -0.3427192
skewness(ph_bin2_noseason_sim) # 0.3441104
var(ph_bin2_noseason_sim) # 0.1266226
mean(ph_bin2_noseason_sim) # 0.6742878


############################# bin 3
mu_ph_bin3_noseason   <-     0.7352998
sigma_ph_bin3_noseason <-   0.5809116
nu_ph_bin3_noseason  <-    1.1691256
tau_ph_bin3_noseason  <- 2.651172

set.seed(33)
ph_bin3_noseason_sim <- rBCPE(N, mu = mu_ph_bin3_noseason, sigma = sigma_ph_bin3_noseason,
                              nu = nu_ph_bin3_noseason, tau = tau_ph_bin3_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin3_noseason_sim) - 3 # -0.6733169
skewness(ph_bin3_noseason_sim) # 0.09022965
var(ph_bin3_noseason_sim) # 0.1437281
mean(ph_bin3_noseason_sim) # 0.7784258


############################# bin 4
mu_ph_bin4_noseason   <-     0.8010986
sigma_ph_bin4_noseason <-   0.6139708
nu_ph_bin4_noseason  <-    1.3181606
tau_ph_bin4_noseason  <- 3.007827

set.seed(44)
ph_bin4_noseason_sim <- rBCPE(N, mu = mu_ph_bin4_noseason, sigma = sigma_ph_bin4_noseason,
                              nu = nu_ph_bin4_noseason, tau = tau_ph_bin4_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin4_noseason_sim) - 3 # -0.7806745
skewness(ph_bin4_noseason_sim) # -0.01996173
var(ph_bin4_noseason_sim) # 0.1706082
mean(ph_bin4_noseason_sim) # 0.8733238


############################# bin 5
mu_ph_bin5_noseason   <-     0.8050691
sigma_ph_bin5_noseason <-   0.6584863
nu_ph_bin5_noseason  <-    1.3327172
tau_ph_bin5_noseason  <-  2.852733

set.seed(55)
ph_bin5_noseason_sim <- rBCPE(N, mu = mu_ph_bin5_noseason, sigma = sigma_ph_bin5_noseason,
                              nu = nu_ph_bin5_noseason, tau = tau_ph_bin5_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin5_noseason_sim) - 3 # -0.7559807
skewness(ph_bin5_noseason_sim) # 0.008365335
var(ph_bin5_noseason_sim) # 0.1841075
mean(ph_bin5_noseason_sim) # 0.8990165


############################# bin 6
mu_ph_bin6_noseason   <-     0.7365117
sigma_ph_bin6_noseason <-   0.6865062
nu_ph_bin6_noseason  <-    1.2590187
tau_ph_bin6_noseason  <- 2.631527

set.seed(66)
ph_bin6_noseason_sim <- rBCPE(N, mu = mu_ph_bin6_noseason, sigma = sigma_ph_bin6_noseason,
                              nu = nu_ph_bin6_noseason, tau = tau_ph_bin6_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin6_noseason_sim) - 3 # -0.694727
skewness(ph_bin6_noseason_sim) # 0.09545519
var(ph_bin6_noseason_sim) # 0.1690114
mean(ph_bin6_noseason_sim) # 0.8294085


############################# bin 7
mu_ph_bin7_noseason   <-     0.6169694
sigma_ph_bin7_noseason <-   0.6915215
nu_ph_bin7_noseason  <-    1.1317877
tau_ph_bin7_noseason  <- 2.314725

set.seed(77)
ph_bin7_noseason_sim <- rBCPE(N, mu = mu_ph_bin7_noseason, sigma = sigma_ph_bin7_noseason,
                              nu = nu_ph_bin7_noseason, tau = tau_ph_bin7_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin7_noseason_sim) - 3 # -0.5221081
skewness(ph_bin7_noseason_sim) # 0.2393307
var(ph_bin7_noseason_sim) # 0.1286196
mean(ph_bin7_noseason_sim) # 0.688109

 
############################# bin 8
mu_ph_bin8_noseason   <-     0.5352415
sigma_ph_bin8_noseason <-   0.7356121
nu_ph_bin8_noseason  <-    0.9911742
tau_ph_bin8_noseason  <- 2.034754

set.seed(88)
ph_bin8_noseason_sim <- rBCPE(N, mu = mu_ph_bin8_noseason, sigma = sigma_ph_bin8_noseason,
                              nu = nu_ph_bin8_noseason, tau = tau_ph_bin8_noseason)

# get skewness and kurtosis 
kurtosis(ph_bin8_noseason_sim) - 3 # -0.1899719
skewness(ph_bin8_noseason_sim) # 0.4386638
var(ph_bin8_noseason_sim) # 0.1147229
mean(ph_bin8_noseason_sim) # 0.6038945


#=====================================================================================================

# linear predictor 

ph_all_BCPE_bins_linear_predict <- list(
  "1" = rBCPE(656, mu = 0.7792665, sigma = 0.4596888, nu = 0.4367311, tau = 2.280653),
  "2" = rBCPE(532, mu = 0.7488316, sigma = 0.4632582, nu = 0.4456891, tau = 2.183464),
  "3" = rBCPE(562, mu = 0.7483788, sigma = 0.4501536, nu = 0.4644804, tau = 2.072435),
  "4" = rBCPE(862, mu = 0.7327450, sigma = 0.4477479, nu = 0.4793662, tau = 1.962916),
  "5" = rBCPE(959, mu = 0.7248292, sigma = 0.4404432, nu = 0.5036419, tau = 1.871124),
  "6" = rBCPE(980, mu = 0.7132266, sigma = 0.4351375, nu = 0.5217891, tau = 1.780893),
  "7" = rBCPE(973, mu = 0.7069516, sigma = 0.4268789, nu = 0.5375946, tau = 1.692632),
  "8" = rBCPE(968, mu = 0.6877554, sigma = 0.4245272, nu = 0.5522392, tau = 1.606726)
)


############################## bin 1
mu_ph_bin1_linear   <-     0.7792665
sigma_ph_bin1_linear  <-    0.4596888
nu_ph_bin1_linear   <-    0.4367311
tau_ph_bin1_linear   <- 2.280653

set.seed(111)
N <- 5e5                    
ph_bin1_linear_sim <- rBCPE(N, mu = mu_ph_bin1_linear, sigma = sigma_ph_bin1_linear, 
                              nu = nu_ph_bin1_linear , tau = tau_ph_bin1_linear)

# get skewness and kurtosis 
kurtosis(ph_bin1_linear_sim) - 3 # 0.3711447
skewness(ph_bin1_linear_sim) # 0.6693149
var(ph_bin1_linear_sim) # 0.1342779
mean(ph_bin1_linear_sim) # 0.8257801

############################# bin 2
mu_ph_bin2_linear  <-     0.7488316
sigma_ph_bin2_linear <-    0.4632582
nu_ph_bin2_linear  <-    0.4456891
tau_ph_bin2_linear  <-  2.183464

set.seed(222)
ph_bin2_linear_sim <- rBCPE(N, mu = mu_ph_bin2_linear, sigma = sigma_ph_bin2_linear, 
                              nu = nu_ph_bin2_linear, tau = tau_ph_bin2_linear)

# get skewness and kurtosis 
kurtosis(ph_bin2_linear_sim) - 3 # 0.481956
skewness(ph_bin2_linear_sim) # 0.6879231
var(ph_bin2_linear_sim) # 0.1256405
mean(ph_bin2_linear_sim) # 0.7939405


############################# bin 3
mu_ph_bin3_linear   <-     0.7483788
sigma_ph_bin3_linear <-   0.4501536
nu_ph_bin3_linear  <-    0.4644804
tau_ph_bin3_linear  <- 2.072435

set.seed(333)
ph_bin3_linear_sim <- rBCPE(N, mu = mu_ph_bin3_linear, sigma = sigma_ph_bin3_linear,
                              nu = nu_ph_bin3_linear, tau = tau_ph_bin3_linear)

# get skewness and kurtosis 
kurtosis(ph_bin3_linear_sim) - 3 # 0.586546
skewness(ph_bin3_linear_sim) # 0.6835901
var(ph_bin3_linear_sim) # 0.1172834
mean(ph_bin3_linear_sim) # 0.7887759


############################# bin 4
mu_ph_bin4_linear   <-     0.7327450
sigma_ph_bin4_linear <-   0.4477479
nu_ph_bin4_linear  <-    0.4793662
tau_ph_bin4_linear  <- 1.962916

set.seed(444)
ph_bin4_linear_sim <- rBCPE(N, mu = mu_ph_bin4_linear, sigma = sigma_ph_bin4_linear,
                              nu = nu_ph_bin4_linear, tau = tau_ph_bin4_linear)

# get skewness and kurtosis 
kurtosis(ph_bin4_linear_sim) - 3 # 0.6842439
skewness(ph_bin4_linear_sim) # 0.6944822
var(ph_bin4_linear_sim) # 0.1111389
mean(ph_bin4_linear_sim) # 0.7706508


############################# bin 5
mu_ph_bin5_linear   <-     0.7248292
sigma_ph_bin5_linear <-   0.4404432
nu_ph_bin5_linear  <-    0.5036419
tau_ph_bin5_linear  <-  1.871124

set.seed(555)
ph_bin5_linear_sim <- rBCPE(N, mu = mu_ph_bin5_linear, sigma = sigma_ph_bin5_linear,
                              nu = nu_ph_bin5_linear, tau = tau_ph_bin5_linear)

# get skewness and kurtosis 
kurtosis(ph_bin5_linear_sim) - 3 # 0.7674997
skewness(ph_bin5_linear_sim) # 0.6855897
var(ph_bin5_linear_sim) # 0.1046883
mean(ph_bin5_linear_sim) # 0.760207


############################# bin 6
mu_ph_bin6_linear   <-     0.7132266
sigma_ph_bin6_linear <-   0.4351375
nu_ph_bin6_linear  <-    0.5217891
tau_ph_bin6_linear  <- 1.780893

set.seed(666)
ph_bin6_linear_sim <- rBCPE(N, mu = mu_ph_bin6_linear, sigma = sigma_ph_bin6_linear,
                              nu = nu_ph_bin6_linear, tau = tau_ph_bin6_linear)

# get skewness and kurtosis 
kurtosis(ph_bin6_linear_sim) - 3 # 0.8711347
skewness(ph_bin6_linear_sim) # 0.6806238
var(ph_bin6_linear_sim) # 0.09780262
mean(ph_bin6_linear_sim) # 0.7457044


############################# bin 7
mu_ph_bin7_linear   <-     0.7069516
sigma_ph_bin7_linear <-   0.4268789
nu_ph_bin7_linear  <-    0.5375946
tau_ph_bin7_linear  <- 1.692632

set.seed(777)
ph_bin7_linear_sim <- rBCPE(N, mu = mu_ph_bin7_linear, sigma = sigma_ph_bin7_linear,
                              nu = nu_ph_bin7_linear, tau = tau_ph_bin7_linear)

# get skewness and kurtosis 
kurtosis(ph_bin7_linear_sim) - 3 # 1.005653
skewness(ph_bin7_linear_sim) # 0.6887053
var(ph_bin7_linear_sim) # 0.09244795
mean(ph_bin7_linear_sim) # 0.7372133


############################# bin 8
mu_ph_bin8_linear   <-     0.6877554
sigma_ph_bin8_linear <-   0.4245272
nu_ph_bin8_linear  <-    0.5522392
tau_ph_bin8_linear  <- 1.606726

set.seed(888)
ph_bin8_linear_sim <- rBCPE(N, mu = mu_ph_bin8_linear, sigma = sigma_ph_bin8_linear,
                              nu = nu_ph_bin8_linear, tau = tau_ph_bin8_linear)

# get skewness and kurtosis 
kurtosis(ph_bin8_linear_sim) - 3 # 1.262481
skewness(ph_bin8_linear_sim) # 0.7239299
var(ph_bin8_linear_sim) # 0.08651441
mean(ph_bin8_linear_sim) # 0.7157071



#=====================================================================================================

# quadratic predictor 

ph_all_BCPE_bins_quad_predict <- list(
  "1" = rBCPE(656, mu = 0.5837121, sigma = 0.5280267, nu = 0.3540376, tau = 1.511927),
  "2" = rBCPE(532, mu = 0.6934904, sigma = 0.4427230, nu = 0.3115548, tau = 1.809922),
  "3" = rBCPE(562, mu = 0.7904432, sigma = 0.3807979, nu = 0.2872804, tau = 1.996296),
  "4" = rBCPE(862, mu = 0.8282551, sigma = 0.3571958, nu = 0.2875039, tau = 2.081026),
  "5" = rBCPE(959, mu = 0.8325877, sigma = 0.3552568, nu = 0.3100270, tau = 2.038800),
  "6" = rBCPE(980, mu = 0.7905541, sigma = 0.3800358, nu = 0.3485406, tau = 1.865902),
  "7" = rBCPE(973, mu = 0.7112296, sigma = 0.4310772, nu = 0.4017649, tau = 1.597868),
  "8" = rBCPE(968, mu = 0.5766264, sigma = 0.5314282, nu = 0.4849630, tau = 1.285604)
)

############################## bin 1
mu_ph_bin1_quad   <-     0.5837121
sigma_ph_bin1_quad   <-    0.5280267
nu_ph_bin1_quad    <-    0.3540376
tau_ph_bin1_quad    <- 1.511927

set.seed(1111)
N <- 5e5                    
ph_bin1_quad_sim <- rBCPE(N, mu = mu_ph_bin1_quad , sigma = sigma_ph_bin1_quad , 
                            nu = nu_ph_bin1_quad  , tau = tau_ph_bin1_quad )

# get skewness and kurtosis 
kurtosis(ph_bin1_quad_sim) - 3 # 4.057714
skewness(ph_bin1_quad_sim) # 1.390423
var(ph_bin1_quad_sim) # 0.1099008
mean(ph_bin1_quad_sim) # 0.6367964

############################# bin 2
mu_ph_bin2_quad  <-     0.6934904
sigma_ph_bin2_quad <-    0.4427230
nu_ph_bin2_quad  <-    0.3115548
tau_ph_bin2_quad  <-  1.809922

set.seed(2222)
ph_bin2_quad_sim <- rBCPE(N, mu = mu_ph_bin2_quad, sigma = sigma_ph_bin2_quad, 
                            nu = nu_ph_bin2_quad, tau = tau_ph_bin2_quad)

# get skewness and kurtosis 
kurtosis(ph_bin2_quad_sim) - 3 # 1.867631
skewness(ph_bin2_quad_sim) # 1.005016
var(ph_bin2_quad_sim) # 0.1045187
mean(ph_bin2_quad_sim) # 0.7407339


############################# bin 3
mu_ph_bin3_quad   <-     0.7904432
sigma_ph_bin3_quad <-   0.3807979
nu_ph_bin3_quad  <-    0.2872804
tau_ph_bin3_quad  <- 1.996296

set.seed(3333)
ph_bin3_quad_sim <- rBCPE(N, mu = mu_ph_bin3_quad, sigma = sigma_ph_bin3_quad,
                            nu = nu_ph_bin3_quad, tau = tau_ph_bin3_quad)

# get skewness and kurtosis 
kurtosis(ph_bin3_quad_sim) - 3 # 1.068061
skewness(ph_bin3_quad_sim) # 0.8150953
var(ph_bin3_quad_sim) # 0.09843538
mean(ph_bin3_quad_sim) # 0.8318957


############################# bin 4
mu_ph_bin4_quad   <-     0.8282551
sigma_ph_bin4_quad <-   0.3571958
nu_ph_bin4_quad  <-    0.2875039
tau_ph_bin4_quad  <- 2.081026

set.seed(4444)
ph_bin4_quad_sim <- rBCPE(N, mu = mu_ph_bin4_quad, sigma = sigma_ph_bin4_quad,
                            nu = nu_ph_bin4_quad, tau = tau_ph_bin4_quad)

# get skewness and kurtosis 
kurtosis(ph_bin4_quad_sim) - 3 # 0.7366888
skewness(ph_bin4_quad_sim) # 0.7233398
var(ph_bin4_quad_sim) # 0.0933613
mean(ph_bin4_quad_sim) # 0.8660879


############################# bin 5
mu_ph_bin5_quad   <-     0.8325877
sigma_ph_bin5_quad <-   0.3552568
nu_ph_bin5_quad  <-    0.3100270
tau_ph_bin5_quad  <-  2.038800

set.seed(5555)
ph_bin5_quad_sim <- rBCPE(N, mu = mu_ph_bin5_quad, sigma = sigma_ph_bin5_quad,
                            nu = nu_ph_bin5_quad, tau = tau_ph_bin5_quad)

# get skewness and kurtosis 
kurtosis(ph_bin5_quad_sim) - 3 # 0.8035028
skewness(ph_bin5_quad_sim) # 0.7201877
var(ph_bin5_quad_sim) # 0.09261315
mean(ph_bin5_quad_sim) # 0.86844


############################# bin 6
mu_ph_bin6_quad   <-     0.7905541
sigma_ph_bin6_quad <-   0.3800358
nu_ph_bin6_quad  <-    0.3485406
tau_ph_bin6_quad  <- 1.865902

set.seed(6666)
ph_bin6_quad_sim <- rBCPE(N, mu = mu_ph_bin6_quad, sigma = sigma_ph_bin6_quad,
                            nu = nu_ph_bin6_quad, tau = tau_ph_bin6_quad)

# get skewness and kurtosis 
kurtosis(ph_bin6_quad_sim) - 3 # 1.141508
skewness(ph_bin6_quad_sim) # 0.7891089
var(ph_bin6_quad_sim) # 0.09619799
mean(ph_bin6_quad_sim) # 0.8280794


############################# bin 7
mu_ph_bin7_quad   <-     0.7112296
sigma_ph_bin7_quad <-   0.4310772
nu_ph_bin7_quad  <-    0.4017649
tau_ph_bin7_quad  <- 1.597868

set.seed(7777)
ph_bin7_quad_sim <- rBCPE(N, mu = mu_ph_bin7_quad, sigma = sigma_ph_bin7_quad,
                            nu = nu_ph_bin7_quad, tau = tau_ph_bin7_quad)

# get skewness and kurtosis 
kurtosis(ph_bin7_quad_sim) - 3 # 2.066997
skewness(ph_bin7_quad_sim) # 0.9681598
var(ph_bin7_quad_sim) # 0.09993119
mean(ph_bin7_quad_sim) # 0.7493972


############################# bin 8
mu_ph_bin8_quad   <-     0.5766264
sigma_ph_bin8_quad <-   0.5314282
nu_ph_bin8_quad  <-    0.4849630
tau_ph_bin8_quad  <- 1.285604

set.seed(8888)
ph_bin8_quad_sim <- rBCPE(N, mu = mu_ph_bin8_quad, sigma = sigma_ph_bin8_quad,
                            nu = nu_ph_bin8_quad, tau = tau_ph_bin8_quad)

# get skewness and kurtosis 
kurtosis(ph_bin8_quad_sim) - 3 # 4.451172
skewness(ph_bin8_quad_sim) # 1.358754
var(ph_bin8_quad_sim) # 0.1000485
mean(ph_bin8_quad_sim) # 0.6184842






##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 

set.seed(1)
bootstrap_moments_sim_only <- function(y, B = 2000) {
  n <- length(y)
  sk <- numeric(B); exk <- numeric(B); mu <- numeric(B); s2 <- numeric(B)
  for(b in 1:B) {
    samp <- sample(y, size = n, replace = TRUE)
    mu[b] <- mean(samp)
    s2[b] <- var(samp)              # sample variance (n-1)
    sk[b] <- mean((samp - mean(samp))^3)/sd(samp)^3
    exk[b] <- mean((samp - mean(samp))^4)/sd(samp)^4 - 3
  }
  list(
    se = c(mean = sd(mu), var = sd(s2), skew = sd(sk), exkurt = sd(exk)),
    ci = rbind(
      mean = quantile(mu, c(0.025,0.975)),
      var  = quantile(s2, c(0.025,0.975)),
      skew = quantile(sk, c(0.025,0.975)),
      exk  = quantile(exk, c(0.025,0.975))
    ),
    boot_reps = list(mean = mu, var = s2, skew = sk, exkurt = exk)
  )
}

# get them errors gurll (for seasonaility)
ph_bin1_sim_se <- bootstrap_moments_sim_only(ph_bin1_sim) 
ph_bin2_sim_se <- bootstrap_moments_sim_only(ph_bin2_sim) 
ph_bin3_sim_se <- bootstrap_moments_sim_only(ph_bin3_sim) 
ph_bin4_sim_se <- bootstrap_moments_sim_only(ph_bin4_sim) 
ph_bin5_sim_se <- bootstrap_moments_sim_only(ph_bin5_sim) 
ph_bin6_sim_se <- bootstrap_moments_sim_only(ph_bin6_sim) 
ph_bin7_sim_se <- bootstrap_moments_sim_only(ph_bin7_sim) 
ph_bin8_sim_se <- bootstrap_moments_sim_only(ph_bin8_sim) 

ph_bin1_sim_se$se

ph_bins_df <- tibble::tibble(
  year_mid = c(1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  mean = c(0.647,0.723,0.816,0.883,0.894,0.815,0.694,0.640),
  variance = c(0.109,0.105,0.089,0.090,0.102,0.095,0.080,0.097),
  skewness = c(1.351,1.052,0.722,0.657,0.723,0.770,0.933,1.314),
  kurtosis = c(3.716,2.716,0.664,0.407,0.750,1.202,2.115,3.937)
)

ph_se_df <- tibble::tibble(
  year_mid = ph_bins_df$year_mid,
  se_mean = c(0.0004739638,0.0004553520,0.0004190247,0.0004307061,
              0.0004425233,0.0004387041,0.000402048,0.0004405365),
  se_variance = c(0.0003669594,0.0002975767,0.0002083750,0.0002004666,
                  0.0002367871,0.0002410548,0.000227958,0.0003337661),
  se_skew = c(0.0108495722,0.0076980327,0.0043515345,0.0040797363,
              0.0046396849,0.0056142459,0.007746112,0.0111093414),
  se_kurt = c(0.0910671947,0.0556927519,0.0187235385,0.0169571335,
              0.0214152375,0.0287815614,0.048765570,0.0922350562)
)

# join and pivot to long form
ph_df_long <- ph_bins_df %>%
  left_join(ph_se_df, by = "year_mid") %>%
  pivot_longer(cols = c(mean, variance, skewness, kurtosis),
               names_to = "stat", values_to = "value") %>%
  # attach the corresponding se column based on stat name
  rowwise() %>%
  mutate(se = case_when(
    stat == "mean"     ~ se_mean,
    stat == "variance" ~ se_variance,
    stat == "skewness" ~ se_skew,
    stat == "kurtosis" ~ se_kurt,
    TRUE ~ NA_real_
  )) %>%
  ungroup() %>%
  # compute 95% normal-approx CIs; if you already have lo/hi use them
  mutate(lo = value - 1.96 * se,
         hi = value + 1.96 * se)


###########

### standard error for no seasonality
ph_bin1_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin1_noseason_sim) 
ph_bin2_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin2_noseason_sim) 
ph_bin3_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin3_noseason_sim) 
ph_bin4_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin4_noseason_sim) 
ph_bin5_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin5_noseason_sim) 
ph_bin6_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin6_noseason_sim) 
ph_bin7_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin7_noseason_sim) 
ph_bin8_noseason_sim_se <- bootstrap_moments_sim_only(ph_bin8_noseason_sim) 

ph_bin8_noseason_sim_se$se

### standard error for linear
ph_bin1_linear_sim_se <- bootstrap_moments_sim_only(ph_bin1_linear_sim) 
ph_bin2_linear_sim_se <- bootstrap_moments_sim_only(ph_bin2_linear_sim)
ph_bin3_linear_sim_se <- bootstrap_moments_sim_only(ph_bin3_linear_sim)
ph_bin4_linear_sim_se <- bootstrap_moments_sim_only(ph_bin4_linear_sim)
ph_bin5_linear_sim_se <- bootstrap_moments_sim_only(ph_bin5_linear_sim)
ph_bin6_linear_sim_se <- bootstrap_moments_sim_only(ph_bin6_linear_sim)
ph_bin7_linear_sim_se <- bootstrap_moments_sim_only(ph_bin7_linear_sim)
ph_bin8_linear_sim_se <- bootstrap_moments_sim_only(ph_bin8_linear_sim)

ph_bin8_linear_sim_se$se


### standard error for quadratic
ph_bin1_quad_sim_se <- bootstrap_moments_sim_only(ph_bin1_quad_sim) 
ph_bin2_quad_sim_se <- bootstrap_moments_sim_only(ph_bin2_quad_sim) 
ph_bin3_quad_sim_se <- bootstrap_moments_sim_only(ph_bin3_quad_sim) 
ph_bin4_quad_sim_se <- bootstrap_moments_sim_only(ph_bin4_quad_sim) 
ph_bin5_quad_sim_se <- bootstrap_moments_sim_only(ph_bin5_quad_sim) 
ph_bin6_quad_sim_se <- bootstrap_moments_sim_only(ph_bin6_quad_sim) 
ph_bin7_quad_sim_se <- bootstrap_moments_sim_only(ph_bin7_quad_sim) 
ph_bin8_quad_sim_se <- bootstrap_moments_sim_only(ph_bin8_quad_sim) 

ph_bin8_quad_sim_se$se


kurt_phosphate <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Seasonaility_k = c(3.716,2.716,0.664,0.407,0.750,1.202,2.115,3.937),
  No_seasonality_k = c(-0.019,-0.343,-0.673,-0.781,-0.756,-0.695,-0.522,-0.190),
  Linear_k = c(0.371,0.482,0.587,0.684,0.767,0.871,1.006,1.262),
  Quadratic_k = c(4.058,1.867,1.068,0.737,0.804,1.142,2.067,4.451),
)

kurt_phosphate_se <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Seasonaility_k_se = c(0.0910671947,0.0556927519,0.0187235385,0.0169571335,
               0.0214152375,0.0287815614,0.048765570,0.0922350562),
  No_seasonality_k_se = c(0.0099786871,0.0065190524,0.0034419759,0.0028188508,
                 0.0029500065,0.0033937414,0.0047913039,0.0078501632),
  Linear_k_se = c(0.014496574,0.0161810675,0.0179418010,0.0188053116,
                  0.0202831486,0.0238863533,0.0245856798,0.0297939926),
  Quadratic_k_se = c(0.1104776192,0.0402800909,0.0272042724,0.0205674841,
                0.0210399117,0.0261173161,0.0494359997,0.1369690735),
)


kurt_phosphate_df <- merge(kurt_phosphate, kurt_phosphate_se, by = "year_mid")

kurt_phosphate_long <- kurt_phosphate_df %>%
  pivot_longer(
    cols = -year_mid,
    names_to = c("predictor", "metric"),
    names_pattern = "^(.*)_(k(?:_se)?)$",
    values_drop_na = FALSE
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(Kurtosis = k, se = k_se) %>%
  mutate(predictor = str_remove(predictor, "_$")) %>%
  arrange(predictor, year_mid)



# plot that sucka
ggplot(kurt_phosphate_long, aes(x = year_mid, y = Kurtosis, color = predictor, group = predictor)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = Kurtosis - 1.96 * se, ymax = Kurtosis + 1.96 * se),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  scale_color_manual(values = c("Seasonaility" = "darkolivegreen", "No_seasonality" = "darkseagreen",
                                "Linear" = "chocolate", "Quadratic" = "chocolate4")) +
  labs(x = "Year bins", y = "Excess kurtosis", color = "Predictors") +
  scale_x_continuous(breaks = unique(kurt_phosphate_long$year_mid)) +
  #scale_x_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              #"1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "right") 



