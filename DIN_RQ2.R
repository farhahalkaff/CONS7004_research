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



################

# METHOD 1 =====
logSurv(DIN$DIN, prob=0.90, tail="right")


# METHOD 2 =====
par(mfrow = c(1, 3))
DIN_m1 <- loglogSurv1(DIN$DIN, prob=0.90, title="(a) TYPE I")
DIN_m2 <- loglogSurv2(DIN$DIN, prob=0.90, title="(b) TYPE II")
DIN_m3 <- loglogSurv3(DIN$DIN, prob=0.90, title="(c) TYPE III")
par(mfrow = c(1, 1))


# METHOD 3 ======

# use some truncated family distribution 
DIN_m4 <- fitTail(DIN$DIN, family=WEI, percentage=10)
DIN_m5 <- fitTail(DIN$DIN, family=LOGNO, percentage=10)
DIN_m6 <- fitTail(DIN$DIN, family=BCPE, percentage=10)
DIN_m7 <- fitTail(DIN$DIN, family=BCPEo, percentage=10)
DIN_m8 <- fitTail(DIN$DIN, family=BCT, percentage=10) # warning 
DIN_m9 <- fitTail(DIN$DIN, family=BCTo, percentage=10) # warning 
DIN_m10 <- fitTail(DIN$DIN, family=JSU, percentage=10) # warning 


# check AIC 
AIC(DIN_m4, DIN_m5, DIN_m6, DIN_m7, DIN_m8, DIN_m9, DIN_m10)
## NOTE: BCPE does the best

wp(DIN_m4, ylim.all = 1)
DIN_m4_2 <- fitTailAll(DIN$DIN, family=WEI) # got like 50 warnings
plot(DIN_m4_2)




# 1. Total number of observations
total_obs <- length(DIN$DIN)
print(paste("Total observations:", total_obs))

# 2. Number of observations in the top 10%
num_top_10_percent <- ceiling(total_obs * 0.10)
print(paste("Number of observations in the top 10%:", num_top_10_percent))




# empirical estimate of moments
DIN_moments_summary <- DIN %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(DIN, na.rm = TRUE),
    var = var(DIN, na.rm = TRUE),
    skew = skewness(DIN, na.rm = TRUE),
    kurt = kurtosis(DIN, na.rm = TRUE)
  )
print(DIN_moments_summary, n=33)



















## adding smoothing spline


DIN_all_BCPE_spline <- gamlss(DIN ~ cs(year) + month, sigma.fo = ~ cs(year) + month, nu.fo = ~ cs(year) + month, tau.fo = ~ cs(year) + month,
                                family = BCPE(), data = DIN,
                                mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                #method = mixed(5,100),
                                control = gamlss.control(n.cyc = 1000, c.crit = 0.01, trace = TRUE))


DIN_all_BCPEo_spline_param <- best_model_param(DIN_all_BCPE_spline, DIN)

set.seed(123)
DIN_spline_predict <- list(
  A = rBCPEo(656, mu = 0.5531981, sigma = 0.9942300, nu = 1.0523044, tau = 5.746193),
  B = rBCPEo(532, mu = 0.6605832, sigma = 0.6494657, nu = 0.3551003, tau = 2.101785),
  C = rBCPEo(562, mu = 0.7260273, sigma = 1.0089550, nu = 1.1233813, tau = 3.132722),
  D = rBCPEo(862, mu = 0.6154539, sigma = 1.6903863, nu = 1.9191222, tau = 1.468543),
  E = rBCPEo(959, mu = 0.8598428, sigma = 0.6254443, nu = 0.6322192, tau = 2.734192),
  F = rBCPEo(980, mu = 0.7792452, sigma = 0.5956562, nu = 0.3224346, tau = 3.898431),
  G = rBCPEo(973, mu = 0.6595062, sigma = 0.6746916, nu = 0.5138789, tau = 5.569839),
  H = rBCPEo(968, mu = 0.6165068, sigma = 0.6516583, nu = 0.6742924, tau = 8.234633)
)


DIN_spline_predict_long <- stack(DIN_spline_predict)
names(DIN_spline_predict_long) <- c("pred", "period")


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = DIN, aes(x = DIN, y = period), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = DIN_spline_predict_long, aes(x = pred, y = period), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  







##### use bins in the model 
DIN$period <- as.numeric(DIN$period)

DIN_all_BCPE_bins_spline <- gamlss(DIN  ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                    family = BCPE(), data = DIN,
                                    #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                    method = mixed(5,100),
                                    control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


DIN_all_BCPE_bins_spline_noseason <- gamlss(DIN  ~ cs(period), sigma.fo = ~ cs(period), nu.fo = ~ cs(period), tau.fo = ~ cs(period),
                                   family = BCPE(), data = DIN,
                                   #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))



DIN_all_BCPE_bins_spline_param <- best_model_param(DIN_all_BCPE_bins_spline, DIN)
DIN_all_BCPE_bins_noseason_spline_param <- best_model_param(DIN_all_BCPE_bins_spline_noseason, DIN)

## with seasonality
set.seed(123)
DIN_all_BCPE_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 17.16373, sigma = 0.4238357, nu = 0.5371271, tau = 2.787200),
  "2" = rBCPE(532, mu = 16.99580, sigma = 0.4475782, nu = 0.4465462, tau = 2.235752),
  "3" = rBCPE(562, mu = 17.09090, sigma = 0.4385957, nu = 0.3565606, tau = 1.993635),
  "4" = rBCPE(862, mu = 18.25920, sigma = 0.3932305, nu = 0.1813605, tau = 1.895521),
  "5" = rBCPE(959, mu = 20.18476, sigma = 0.3633727, nu = 0.1823012, tau = 1.792768),
  "6" = rBCPE(980, mu = 23.09552, sigma = 0.3759787, nu = 0.3308157, tau = 1.951908),
  "7" = rBCPE(973, mu = 26.43533, sigma = 0.4575213, nu = 0.3221296, tau = 2.523017),
  "8" = rBCPE(968, mu = 28.21057, sigma = 0.6696904, nu = 0.1331123, tau = 4.439182)
)


DIN_all_BCPE_bins_spline_predict_long <- stack(DIN_all_BCPE_bins_spline_predict)
names(DIN_all_BCPE_bins_spline_predict_long) <- c("pred", "period")
DIN_all_BCPE_bins_spline_predict_long$period <- as.numeric(DIN_all_BCPE_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = DIN, aes(x = DIN, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = DIN_all_BCPE_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "DIN (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  


## without seasonality
set.seed(123)
DIN_all_BCPE_bins_noseason_spline_predict <- list(
  "1" = rBCPE(656, mu = 15.35985, sigma = 0.4722517, nu = -0.248743410, tau = 1.576428),
  "2" = rBCPE(532, mu = 15.13255, sigma = 0.5348510, nu = 0.002369268, tau = 1.712444),
  "3" = rBCPE(562, mu = 15.24094, sigma = 0.5300953, nu = 0.225169428, tau = 1.636829),
  "4" = rBCPE(862, mu = 17.53925, sigma = 0.5012095, nu = 0.225809414, tau = 1.948899),
  "5" = rBCPE(959, mu = 20.25049, sigma = 0.4871626, nu = 0.299247661, tau = 2.430468),
  "6" = rBCPE(980, mu = 23.24447, sigma = 0.5015433, nu = 0.501801058, tau = 2.295098),
  "7" = rBCPE(973, mu = 25.98744, sigma = 0.5796489, nu = 0.491649628, tau = 2.046770),
  "8" = rBCPE(968, mu = 26.18295, sigma = 0.7839106, nu = 0.275310304, tau = 2.414724)
)


DIN_all_BCPE_bins_noseason_spline_predict_long <- stack(DIN_all_BCPE_bins_noseason_spline_predict)
names(DIN_all_BCPE_bins_noseason_spline_predict_long) <- c("pred", "period")
DIN_all_BCPE_bins_noseason_spline_predict_long$period <- as.numeric(DIN_all_BCPE_bins_noseason_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = DIN, aes(x = DIN, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = DIN_all_BCPE_bins_noseason_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "DIN (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  +
  ggplot() +
  geom_density_ridges(data = DIN, aes(x = DIN, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = DIN_all_BCPE_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "DIN (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  







### bins, changing paramter time-varying 

DIN_meanonly_BCPE_bins_spline <- gamlss(DIN  ~ cs(period) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                                   family = BCPE(), data = DIN,
                                   #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


DIN_mean_sigma_BCPE_bins_spline <- gamlss(DIN  ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ 1,
                                        family = BCPE(), data = DIN,
                                        #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                        method = mixed(5,100),
                                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


DIN_nokurt_BCPE_bins_spline <- gamlss(DIN  ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ 1,
                                          family = BCPE(), data = DIN,
                                          #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                          method = mixed(5,100),
                                          control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


DIN_noskew_BCPE_bins_spline <- gamlss(DIN  ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ cs(period) + month,
                                      family = BCPE(), data = DIN,
                                      #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                      method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

DIN_all_BCPE_bins_spline <- gamlss(DIN  ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                   family = BCPE(), data = DIN,
                                   #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


AIC(DIN_meanonly_BCPE_bins_spline, DIN_mean_sigma_BCPE_bins_spline,
    DIN_nokurt_BCPE_bins_spline, DIN_noskew_BCPE_bins_spline, DIN_all_BCPE_bins_spline)

##### RESULTS
# DIN_all_BCPE_bins_spline         45428.31
# DIN_nokurt_BCPE_bins_spline      45457.24
# DIN_noskew_BCPE_bins_spline      45558.09
# DIN_mean_sigma_BCPE_bins_spline  45601.53
# DIN_meanonly_BCPE_bins_spline    46347.40




DIN_intercept_JSU <- gamlss(DIN  ~  1,
                                   family = JSU(), data = DIN,
                                   #mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


moment_bucket(DIN_intercept_JSU) +
resid_wp(DIN_intercept_JSU) + theme(plot.title = element_blank())








## monte carlo simulation

DIN_all_BCPE_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 17.16373, sigma = 0.4238357, nu = 0.5371271, tau = 2.787200),
  "2" = rBCPE(532, mu = 16.99580, sigma = 0.4475782, nu = 0.4465462, tau = 2.235752),
  "3" = rBCPE(562, mu = 17.09090, sigma = 0.4385957, nu = 0.3565606, tau = 1.993635),
  "4" = rBCPE(862, mu = 18.25920, sigma = 0.3932305, nu = 0.1813605, tau = 1.895521),
  "5" = rBCPE(959, mu = 20.18476, sigma = 0.3633727, nu = 0.1823012, tau = 1.792768),
  "6" = rBCPE(980, mu = 23.09552, sigma = 0.3759787, nu = 0.3308157, tau = 1.951908),
  "7" = rBCPE(973, mu = 26.43533, sigma = 0.4575213, nu = 0.3221296, tau = 2.523017),
  "8" = rBCPE(968, mu = 28.21057, sigma = 0.6696904, nu = 0.1331123, tau = 4.439182)
)

library(gamlss.dist)
library(moments)

############################## bin 1
mu_DIN_bin1   <- 17.16373    
sigma_DIN_bin1 <- 0.4238357   
nu_DIN_bin1  <- 0.5371271   
tau_DIN_bin1  <- 2.787200

set.seed(1)
N <- 5e5                    
DIN_bin1_sim <- rBCPE(N, mu = mu_DIN_bin1, sigma = sigma_DIN_bin1, nu = nu_DIN_bin1, tau = tau_DIN_bin1)

# get skewness and kurtosis 
kurtosis(DIN_bin1_sim) - 3 # -0.283515
skewness(DIN_bin1_sim) # 0.4358236
var(DIN_bin1_sim) # 53.37748
mean(DIN_bin1_sim) # 17.87089


############################# bin 2
mu_DIN_bin2   <-    16.99580 
sigma_DIN_bin2 <-    0.4475782
nu_DIN_bin2  <-    0.4465462
tau_DIN_bin2  <- 2.235752

set.seed(2)
DIN_bin2_sim <- rBCPE(N, mu = mu_DIN_bin2, sigma = sigma_DIN_bin2, nu = nu_DIN_bin2, tau = tau_DIN_bin2)

# get skewness and kurtosis 
kurtosis(DIN_bin2_sim) - 3 # 0.3757609
skewness(DIN_bin2_sim) # 0.6503901
var(DIN_bin2_sim) # 60.25346
mean(DIN_bin2_sim) # 17.95159


############################# bin 3
mu_DIN_bin3   <-     17.09090
sigma_DIN_bin3 <-    0.4385957
nu_DIN_bin3  <-    0.3565606
tau_DIN_bin3  <- 1.993635

set.seed(3)
DIN_bin3_sim <- rBCPE(N, mu = mu_DIN_bin3, sigma = sigma_DIN_bin3, nu = nu_DIN_bin3, tau = tau_DIN_bin3)

# get skewness and kurtosis 
kurtosis(DIN_bin3_sim) - 3 # 1.038866
skewness(DIN_bin3_sim) # 0.82783
var(DIN_bin3_sim) # 60.31834
mean(DIN_bin3_sim) # 18.14148


############################# bin 4
mu_DIN_bin4   <-     18.25920
sigma_DIN_bin4 <-    0.3932305
nu_DIN_bin4  <-    0.1813605
tau_DIN_bin4  <- 1.895521

set.seed(4)
DIN_bin4_sim <- rBCPE(N, mu = mu_DIN_bin4, sigma = sigma_DIN_bin4, nu = nu_DIN_bin4, tau = tau_DIN_bin4)

# get skewness and kurtosis 
kurtosis(DIN_bin4_sim) - 3 # 2.00124
skewness(DIN_bin4_sim) # 1.039221
var(DIN_bin4_sim) # 58.97537
mean(DIN_bin4_sim) # 19.42133


############################# bin 5
mu_DIN_bin5   <- 20.18476    
sigma_DIN_bin5 <-    0.3633727
nu_DIN_bin5  <-    0.1823012
tau_DIN_bin5  <- 1.792768

set.seed(5)
DIN_bin5_sim <- rBCPE(N, mu = mu_DIN_bin5, sigma = sigma_DIN_bin5, nu = nu_DIN_bin5, tau = tau_DIN_bin5)

# get skewness and kurtosis 
kurtosis(DIN_bin5_sim) - 3 # 2.122309
skewness(DIN_bin5_sim) # 1.030426
var(DIN_bin5_sim) # 61.10735
mean(DIN_bin5_sim) # 21.28347


############################# bin 6
mu_DIN_bin6   <-     23.09552
sigma_DIN_bin6 <-    0.3759787
nu_DIN_bin6  <-    0.3308157
tau_DIN_bin6  <- 1.951908

set.seed(6)
DIN_bin6_sim <- rBCPE(N, mu = mu_DIN_bin6, sigma = sigma_DIN_bin6, nu = nu_DIN_bin6, tau = tau_DIN_bin6)

# get skewness and kurtosis 
kurtosis(DIN_bin6_sim) - 3 # 1.043993
skewness(DIN_bin6_sim) # 0.7814984
var(DIN_bin6_sim) # 80.59317
mean(DIN_bin6_sim) # 24.21176


############################# bin 7
mu_DIN_bin7   <-     26.43533
sigma_DIN_bin7 <-    0.4575213
nu_DIN_bin7  <-    0.3221296
tau_DIN_bin7  <- 2.523017

set.seed(7)
DIN_bin7_sim <- rBCPE(N, mu = mu_DIN_bin7, sigma = sigma_DIN_bin7, nu = nu_DIN_bin7, tau = tau_DIN_bin7)

# get skewness and kurtosis 
kurtosis(DIN_bin7_sim) - 3 # 0.3952519
skewness(DIN_bin7_sim) # 0.735436
var(DIN_bin7_sim) # 158.5022
mean(DIN_bin7_sim) # 28.32225


############################# bin 8
mu_DIN_bin8   <-     28.21057
sigma_DIN_bin8 <-    0.6696904
nu_DIN_bin8  <-    0.1331123
tau_DIN_bin8  <- 4.439182

set.seed(8)
DIN_bin8_sim <- rBCPE(N, mu = mu_DIN_bin8, sigma = sigma_DIN_bin8, nu = nu_DIN_bin8, tau = tau_DIN_bin8)

# get skewness and kurtosis 
kurtosis(DIN_bin8_sim) - 3 # 0.4191868
skewness(DIN_bin8_sim) # 0.9668979
var(DIN_bin8_sim) # 474.1645
mean(DIN_bin8_sim) # 33.88857



## without seasonality ======================================================================

DIN_all_BCPE_bins_noseason_spline_predict <- list(
  "1" = rBCPE(656, mu = 15.35985, sigma = 0.4722517, nu = -0.248743410, tau = 1.576428),
  "2" = rBCPE(532, mu = 15.13255, sigma = 0.5348510, nu = 0.002369268, tau = 1.712444),
  "3" = rBCPE(562, mu = 15.24094, sigma = 0.5300953, nu = 0.225169428, tau = 1.636829),
  "4" = rBCPE(862, mu = 17.53925, sigma = 0.5012095, nu = 0.225809414, tau = 1.948899),
  "5" = rBCPE(959, mu = 20.25049, sigma = 0.4871626, nu = 0.299247661, tau = 2.430468),
  "6" = rBCPE(980, mu = 23.24447, sigma = 0.5015433, nu = 0.501801058, tau = 2.295098),
  "7" = rBCPE(973, mu = 25.98744, sigma = 0.5796489, nu = 0.491649628, tau = 2.046770),
  "8" = rBCPE(968, mu = 26.18295, sigma = 0.7839106, nu = 0.275310304, tau = 2.414724)
)

############################## bin 1
mu_DIN_bin1_noseason   <- 15.35985    
sigma_DIN_bin1_noseason <- 0.4722517   
nu_DIN_bin1_noseason  <- -0.248743410   
tau_DIN_bin1_noseason  <- 1.576428

set.seed(11)
DIN_bin1_noseason_sim <- rBCPE(N, mu = mu_DIN_bin1_noseason, sigma = sigma_DIN_bin1_noseason,
                               nu = nu_DIN_bin1_noseason, tau = tau_DIN_bin1_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin1_noseason_sim) - 3 # 4133.439
skewness(DIN_bin1_noseason_sim) # 26.87133
var(DIN_bin1_noseason_sim) # 161.9267
mean(DIN_bin1_noseason_sim) # 18.02925


############################# bin 2
mu_DIN_bin2_noseason   <- 15.13255    
sigma_DIN_bin2_noseason <- 0.5348510   
nu_DIN_bin2_noseason  <- 0.002369268   
tau_DIN_bin2_noseason  <- 1.712444

set.seed(22)
DIN_bin2_noseason_sim <- rBCPE(N, mu = mu_DIN_bin2_noseason, sigma = sigma_DIN_bin2_noseason,
                               nu = nu_DIN_bin2_noseason, tau = tau_DIN_bin2_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin2_noseason_sim) - 3 # 13.56433
skewness(DIN_bin2_noseason_sim) # 2.473163
var(DIN_bin2_noseason_sim) # 108.1929
mean(DIN_bin2_noseason_sim) # 17.46857


############################# bin 3
mu_DIN_bin3_noseason   <- 15.24094    
sigma_DIN_bin3_noseason <- 0.5300953   
nu_DIN_bin3_noseason  <- 0.225169428   
tau_DIN_bin3_noseason  <- 1.636829

set.seed(33)
DIN_bin3_noseason_sim <- rBCPE(N, mu = mu_DIN_bin3_noseason, sigma = sigma_DIN_bin3_noseason,
                               nu = nu_DIN_bin3_noseason, tau = tau_DIN_bin3_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin3_noseason_sim) - 3 # 5.167582
skewness(DIN_bin3_noseason_sim) # 1.594645
var(DIN_bin3_noseason_sim) # 82.03065
mean(DIN_bin3_noseason_sim) # 16.91762


############################# bin 4
mu_DIN_bin4_noseason   <- 17.53925    
sigma_DIN_bin4_noseason <- 0.5012095   
nu_DIN_bin4_noseason  <- 0.225809414   
tau_DIN_bin4_noseason  <- 1.948899

set.seed(44)
DIN_bin4_noseason_sim <- rBCPE(N, mu = mu_DIN_bin4_noseason, sigma = sigma_DIN_bin4_noseason,
                               nu = nu_DIN_bin4_noseason, tau = tau_DIN_bin4_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin4_noseason_sim) - 3 # 2.431368
skewness(DIN_bin4_noseason_sim) # 1.198697
var(DIN_bin4_noseason_sim) # 92.08855
mean(DIN_bin4_noseason_sim) # 19.26622


############################# bin 5
mu_DIN_bin5_noseason   <- 20.25049    
sigma_DIN_bin5_noseason <- 0.4871626   
nu_DIN_bin5_noseason  <- 0.299247661   
tau_DIN_bin5_noseason  <- 2.430468

set.seed(55)
DIN_bin5_noseason_sim <- rBCPE(N, mu = mu_DIN_bin5_noseason, sigma = sigma_DIN_bin5_noseason,
                               nu = nu_DIN_bin5_noseason, tau = tau_DIN_bin5_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin5_noseason_sim) - 3 # 0.7242271
skewness(DIN_bin5_noseason_sim) # 0.8370996
var(DIN_bin5_noseason_sim) # 107.6741
mean(DIN_bin5_noseason_sim) # 21.92153


############################# bin 6
mu_DIN_bin6_noseason   <- 23.24447    
sigma_DIN_bin6_noseason <- 0.5015433   
nu_DIN_bin6_noseason  <- 0.501801058   
tau_DIN_bin6_noseason  <- 2.295098

set.seed(66)
DIN_bin6_noseason_sim <- rBCPE(N, mu = mu_DIN_bin6_noseason, sigma = sigma_DIN_bin6_noseason,
                               nu = nu_DIN_bin6_noseason, tau = tau_DIN_bin6_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin6_noseason_sim) - 3 # 0.2841474
skewness(DIN_bin6_noseason_sim) # 0.6428815
var(DIN_bin6_noseason_sim) # 139.7718
mean(DIN_bin6_noseason_sim) # 24.71611


############################# bin 7
mu_DIN_bin7_noseason   <- 25.98744    
sigma_DIN_bin7_noseason <- 0.5796489   
nu_DIN_bin7_noseason  <- 0.491649628   
tau_DIN_bin7_noseason  <- 2.046770

set.seed(77)
DIN_bin7_noseason_sim <- rBCPE(N, mu = mu_DIN_bin7_noseason, sigma = sigma_DIN_bin7_noseason,
                               nu = nu_DIN_bin7_noseason, tau = tau_DIN_bin7_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin7_noseason_sim) - 3 # 0.9180623
skewness(DIN_bin7_noseason_sim) # 0.8405822
var(DIN_bin7_noseason_sim) # 237.5082
mean(DIN_bin7_noseason_sim) # 28.22634


############################# bin 8
mu_DIN_bin8_noseason   <- 26.18295    
sigma_DIN_bin8_noseason <- 0.7839106   
nu_DIN_bin8_noseason  <- 0.275310304   
tau_DIN_bin8_noseason  <- 2.414724

set.seed(88)
DIN_bin8_noseason_sim <- rBCPE(N, mu = mu_DIN_bin8_noseason, sigma = sigma_DIN_bin8_noseason,
                               nu = nu_DIN_bin8_noseason, tau = tau_DIN_bin8_noseason)

# get skewness and kurtosis 
kurtosis(DIN_bin8_noseason_sim) - 3 # 2.601064
skewness(DIN_bin8_noseason_sim) # 1.381829
var(DIN_bin8_noseason_sim) # 561.0685
mean(DIN_bin8_noseason_sim) # 32.11708













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

# get them errors gurll (seasonaility)
DIN_bin1_sim_se <- bootstrap_moments_sim_only(DIN_bin1_sim) 
DIN_bin2_sim_se <- bootstrap_moments_sim_only(DIN_bin2_sim) 
DIN_bin3_sim_se <- bootstrap_moments_sim_only(DIN_bin3_sim) 
DIN_bin4_sim_se <- bootstrap_moments_sim_only(DIN_bin4_sim) 
DIN_bin5_sim_se <- bootstrap_moments_sim_only(DIN_bin5_sim) 
DIN_bin6_sim_se <- bootstrap_moments_sim_only(DIN_bin6_sim) 
DIN_bin7_sim_se <- bootstrap_moments_sim_only(DIN_bin7_sim) 
DIN_bin8_sim_se <- bootstrap_moments_sim_only(DIN_bin8_sim) 

DIN_bin_sim_se$se

# get them errors gurll (no seasonaility)
DIN_bin1_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin1_noseason_sim) 
DIN_bin2_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin2_noseason_sim) 
DIN_bin3_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin3_noseason_sim) 
DIN_bin4_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin4_noseason_sim) 
DIN_bin5_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin5_noseason_sim) 
DIN_bin6_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin6_noseason_sim) 
DIN_bin7_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin7_noseason_sim) 
DIN_bin8_noseason_sim_se <- bootstrap_moments_sim_only(DIN_bin8_noseason_sim) 

DIN_bin8_noseason_sim_se$se

DIN_bins_df <- tibble::tibble(
  year_mid = c(1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  mean = c(17.871,17.951,18.141,19.421,21.283,24.211,28.322,33.889),
  variance = c(53.377,60.253,60.318,58.975,61.107,80.593,158.502,474.165),
  skewness = c(0.436,0.650,0.828,1.039,1.030,0.781,0.735,0.967),
  kurtosis = c(-0.283,0.376,1.038,2.001,2.122,1.044,0.395,0.419)
)

DIN_se_df <- tibble::tibble(
  year_mid = df$year_mid,
  se_mean = c(0.010306205,0.01093010,0.010975780,0.010889395,
              0.011115839,0.01278470,0.018202419,0.03081060),
  se_variance = c(0.098832581,0.13149763,0.149360467,0.167194544,
                  0.177430045,0.19541194,0.349876593,1.03674250),
  se_skew = c(0.002754898,0.00394651,0.005058384,0.008070872,
              0.008002874,0.00535313,0.003689457,0.00320749),
  se_kurt = c(0.006601600,0.01536872,0.025791753,0.069681021,
              0.056682017,0.02646837,0.014115509,0.01174154)
)

# join and pivot to long form
DIN_df_long <- DIN_bins_df %>%
  left_join(DIN_se_df, by = "year_mid") %>%
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




DIN_bins_noseason_df <- tibble::tibble(
  year_mid = c(1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  kurtosis = c(4133.2,13.5,5.2,2.4,0.7,0.3,0.9,2.6),
  se_kurt = c(2.779079e+03,0.75530708,0.18146585,0.059107893,
              0.02008191,0.013860975,0.024428169,0.047358488)
)


ggplot(data = DIN_bins_noseason_df, aes(x = year_mid, y = kurtosis)) +
  geom_line(aes(y = kurtosis), color = "azure4") +
  geom_point(aes(y = kurtosis), color = "azure4") +
  geom_errorbar(aes(ymin = kurtosis - 1.96 * se_kurt, ymax = kurtosis + 1.96 * se_kurt),
                  width = 0.7, color = "azure4") +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  theme_bw()
  
  
  
  
  
  

