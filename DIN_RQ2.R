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
DIN_all_BCT <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCT(), data = DIN,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: all good here

#== BCTo ==#
DIN_all_BCTo <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCTo(), data = DIN,
                      mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got an error saying while (abs(G.dev.in - i.G.dev) > i.c.crit && i.iter < i.n.cyc) { : 
# missing value where TRUE/FALSE needed

#== BCPE ==#
DIN_all_BCPE <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPE(), data = DIN,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== BCPEo ==#
DIN_all_BCPEo <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = BCPEo(), data = DIN,
                       method = mixed(5,100),
                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== GB2 ==#
DIN_all_GB2 <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = GB2(), data = DIN,
                     mu.start = mean(DIN), sigma.start = sd(DIN$DIN),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: got an error saying while (abs(G.dev.in - i.G.dev) > i.c.crit && i.iter < i.n.cyc) { : 
# missing value where TRUE/FALSE needed


# 2. Run models with JSU or SST family distribution ----------------------------------------------------------------------------

#== SST ==#
DIN_all_SST <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SST(), data = DIN,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: sigma must be positive error. always get this error when running SST on some nutrient! not sure why


#== JSU ==#
DIN_all_JSU <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = JSU(), data = DIN,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here!


# 3. check which has the lowest AIC and use that family distribution -----------------------------------------------------------

AIC(DIN_all_BCT)
AIC(DIN_all_BCPE)
AIC(DIN_all_BCPEo)
AIC(DIN_all_JSU)


# the rankings are:
## 1. BCPEo
## 2. BCT
## 3. BCPE
## 4. JSU


# 4. Make sure that family can run all time-varying parameter models ------------------------------------------------------------

## MEAN ONLY MODEL ##
DIN_meanonly_BCPEo <- gamlss(DIN ~ year + month,
                           family = BCPEo(), data = DIN,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
DIN_mean_and_sigma_BCPEo <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month,
                                 family = BCPEo(), data = DIN,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
DIN_contskew_BCPEo <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                           family = BCPEo(), data = DIN,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
DIN_contkurt_BCPEo <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = BCPEo(), data = DIN,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# all works. HOORAY! 


# 5b. IF YES: check AIC for time varying parameters -------------------------------------------------------------------------------

AIC(DIN_meanonly_BCPEo)
AIC(DIN_mean_and_sigma_BCPEo)
AIC(DIN_contskew_BCPEo)
AIC(DIN_contkurt_BCPEo)
AIC(DIN_all_BCPEo)

# rankings are:
## 1. all param (time): 45990.01
## 2. constant kurtosis: 46005.01
## 3. constant skewness: 46041.35
## 4. mean and sigma: 46060.61
## 5. mean only: 46419.5


# 6. Then check how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

DIN$mu_hat_all_BCPEo_DIN <- predict(DIN_all_BCPEo, what = "mu", type = "response")

ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPEo_DIN), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "DIN (µmol/l)") +
  theme_minimal()


# 8. EXTRA: try adding poly() to each parameter to the model and see how that does

##== poly to the mean ==## 
DIN_all_BCPEo_meanpoly <- gamlss(DIN ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCPEo(), data = DIN,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: lower AIC: 45876.68


##== poly to the sigma ==## 
DIN_all_BCPEo_sigmapoly <- gamlss(DIN ~ year + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                family = BCPEo(), data = DIN,
                                method = mixed(5,100),
                                control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: even lower AIC: 45646.67


##== poly to the skewness ==## 
DIN_all_BCPEo_skewpoly <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ year + month,
                               family = BCPEo(), data = DIN,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: a little lower AIC: 45978.39


##== poly to the kurtosis ==## 
DIN_all_BCPE_kurtpoly <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ poly(year,2) + month,
                               family = BCPEo(), data = DIN,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: a little lower AIC: 45905.67


## let's go with sigma poly ##

# 6. RECHECK how the model predicts the mean and plot them to see if it make sense ----------------------------------------------


DIN$mu_hat_all_BCPEo_sigmapoly <- predict(DIN_all_BCPEo_sigmapoly, what = "mu", type = "response")

ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPEo_sigmapoly), color = "goldenrod", linewidth = 1) +
  #geom_line(aes(y = mu_hat_all_BCT), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "DIN (µmol/l)") +
  theme_minimal()

## Looks the same COOL!



# 9. Check PIT histogram, worm plot and moment bucket of the model -----------------------------------------------------------------------------

# fitted param for (am_all_SSTtr) -----------
DIN_mu_hat_all_BCPEo_sigmapoly    <- predict(DIN_all_BCPEo_sigmapoly, "mu", type = "response")
DIN_sigma_hat_all_BCPEo_sigmapoly <- predict(DIN_all_BCPEo_sigmapoly, "sigma", type = "response")
DIN_nu_hat_all_BCPEo_sigmapoly    <- predict(DIN_all_BCPEo_sigmapoly, "nu", type = "response")
DIN_tau_hat_all_BCPEo_sigmapoly   <- predict(DIN_all_BCPEo_sigmapoly, "tau", type = "response")

DIN_pit <- pBCPEo(DIN$DIN, mu = DIN_mu_hat_all_BCPEo_sigmapoly, sigma = DIN_sigma_hat_all_BCPEo_sigmapoly, 
                nu = DIN_nu_hat_all_BCPEo_sigmapoly, tau = DIN_tau_hat_all_BCPEo_sigmapoly)

# Plot PIT histogram
hist(DIN_pit, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(DIN_pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity
## NOTE: hmm not too great

resid_wp(DIN_all_BCPEo_sigmapoly)
## NOTE: very much skewed holy SHIT

moment_bucket(DIN_all_BCPEo_sigmapoly)
## NOTE: very much skewed holy SHIT


# let's check if we don't add polynomial
resid_wp(DIN_all_BCPEo)
## NOTE: looks a bit better, not skewed, but tails not so much 

moment_bucket(DIN_all_BCPEo)
## NOTE: looks better, a bit of platy


##### LET'S STICK WITH NO POLY TO SIGMA ######



# 7b. IF YES: then that's your model baby! ---------------------------------------------------------------------------------------

# now let's check the summary of the model
summary(am_all_BCPE_sigmapoly)

## ==================================================================================================================================


## STEP 2 ==============================================================================================================================

# 1. subset few years together to get period ---------------------------------------------------------------------------------------
# Add column to combine few years together 
DIN <- DIN %>% 
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
DIN_all_BCPEo_param <- best_model_param(DIN_all_BCPEo, DIN)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
DIN_A <- DIN[DIN$period %in% c("A"), ]
DIN_A_pdf <- density(DIN_A$DIN)

# Period B
DIN_B <- DIN[DIN$period %in% c("B"), ]
DIN_B_pdf <- density(DIN_B$DIN)

# Period C
DIN_C <- DIN[DIN$period %in% c("C"), ]
DIN_C_pdf <- density(DIN_C$DIN)

# Period D
DIN_D <- DIN[DIN$period %in% c("D"), ]
DIN_D_pdf <- density(DIN_D$DIN)

# Period E
DIN_E <- DIN[DIN$period %in% c("E"), ]
DIN_E_pdf <- density(DIN_E$DIN)

# Period F
DIN_F <- DIN[DIN$period %in% c("F"), ]
DIN_F_pdf <- density(DIN_F$DIN)

# Period G
DIN_G <- DIN[DIN$period %in% c("G"), ]
DIN_G_pdf <- density(DIN_G$DIN)

# Period H
DIN_H <- DIN[DIN$period %in% c("H"), ]
DIN_H_pdf <- density(DIN_H$DIN)


# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(DIN_all_BCPEo$y, na.rm = TRUE),
         max(DIN_all_BCPEo$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
DIN_all_BCPEo_A_pdf  <- dBCPEo(x, mu = 14.77722, sigma = 0.3939043, nu = 0.2915655, tau = 1.850251)
DIN_all_BCPEo_B_pdf  <- dBCPEo(x, mu = 16.58895, sigma = 0.4080650, nu = 0.2546235, tau = 1.885676)
DIN_all_BCPEo_C_pdf  <- dBCPEo(x, mu = 18.19902, sigma = 0.4250647, nu = 0.2180932, tau = 1.899121)
DIN_all_BCPEo_D_pdf  <- dBCPEo(x, mu = 19.31376, sigma = 0.4448139, nu = 0.2115163, tau = 1.909659)
DIN_all_BCPEo_E_pdf  <- dBCPEo(x, mu = 20.99940, sigma = 0.4608142, nu = 0.1875452, tau = 1.895523)
DIN_all_BCPEo_F_pdf  <- dBCPEo(x, mu = 22.69424, sigma = 0.4802670, nu = 0.1671654, tau = 1.905469)
DIN_all_BCPEo_G_pdf  <- dBCPEo(x, mu = 24.44354, sigma = 0.5004322, nu = 0.1495924, tau = 1.911277)
DIN_all_BCPEo_H_pdf  <- dBCPEo(x, mu = 26.56958, sigma = 0.5215259, nu = 0.1227661, tau = 1.919613)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(DIN_A_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(DIN_all_BCPEo_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(DIN_B_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.07))
lines(DIN_all_BCPEo_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(DIN_C_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(DIN_all_BCPEo_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(DIN_D_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(DIN_all_BCPEo_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(DIN_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(DIN_all_BCPEo_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(DIN_F_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(DIN_all_BCPEo_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(DIN_G_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.04))
lines(DIN_all_BCPEo_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(DIN_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.032))
lines(DIN_all_BCPEo_H_pdf, type = "l", lwd = 2, col = "goldenrod") 

