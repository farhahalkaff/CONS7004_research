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
    between(year, 1962, 1966) ~ "A",
    between(year, 1967, 1970) ~ "B",
    between(year, 1971, 1974) ~ "C",
    between(year, 1975, 1978) ~ "D",
    between(year, 1979, 1982) ~ "E",
    between(year, 1983, 1986) ~ "F",
    between(year, 1987, 1990) ~ "G", 
    between(year, 1991, 1994) ~ "H"
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
Phosphate_A <- Phosphate[Phosphate$period %in% c("A"), ]
Phosphate_A_pdf <- density(Phosphate_A$Phosphate)

# Period B
Phosphate_B <- Phosphate[Phosphate$period %in% c("B"), ]
Phosphate_B_pdf <- density(Phosphate_B$Phosphate)

# Period C
Phosphate_C <- Phosphate[Phosphate$period %in% c("C"), ]
Phosphate_C_pdf <- density(Phosphate_C$Phosphate)

# Period D
Phosphate_D <- Phosphate[Phosphate$period %in% c("D"), ]
Phosphate_D_pdf <- density(Phosphate_D$Phosphate)

# Period E
Phosphate_E <- Phosphate[Phosphate$period %in% c("E"), ]
Phosphate_E_pdf <- density(Phosphate_E$Phosphate)

# Period F
Phosphate_F <- Phosphate[Phosphate$period %in% c("F"), ]
Phosphate_F_pdf <- density(Phosphate_F$Phosphate)

# Period G
Phosphate_G <- Phosphate[Phosphate$period %in% c("G"), ]
Phosphate_G_pdf <- density(Phosphate_G$Phosphate)

# Period H
Phosphate_H <- Phosphate[Phosphate$period %in% c("H"), ]
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







