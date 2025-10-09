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



# 9. Check PIT histogram, worm plot and moment bucket of the model -----------------------------------------------------------------------------

# fitted param for (am_all_SSTtr) -----------
am_mu_hat_all_BCPE_sigmapoly    <- predict(am_all_BCPE_sigmapoly, "mu", type = "response")
am_sigma_hat_all_BCPE_sigmapoly <- predict(am_all_BCPE_sigmapoly, "sigma", type = "response")
am_nu_hat_all_BCPE_sigmapoly    <- predict(am_all_BCPE_sigmapoly, "nu", type = "response")
am_tau_hat_all_BCPE_sigmapoly   <- predict(am_all_BCPE_sigmapoly, "tau", type = "response")

am_pit <- pBCPE(Ammonium$Ammonium, mu = am_mu_hat_all_BCPE_sigmapoly, sigma = am_sigma_hat_all_BCPE_sigmapoly, 
                nu = am_nu_hat_all_BCPE_sigmapoly, tau = am_tau_hat_all_BCPE_sigmapoly)

# Plot PIT histogram
hist(am_pit, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(am_pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity
## NOTE: yeah looks pretty good

resid_wp(am_all_BCPE_sigmapoly)
## NOTE: looks pretty dang good 

moment_bucket(am_all_BCPE_sigmapoly)
## NOTE: looks pretty dang good 


# 7b. IF YES: then that's your model baby! ---------------------------------------------------------------------------------------

# now let's check the summary of the model
summary(am_all_BCPE_sigmapoly)

## ==================================================================================================================================


## STEP 2 ==============================================================================================================================

# 1. subset few years together to get period ---------------------------------------------------------------------------------------
# Add column to combine few years together 
Ammonium <- Ammonium %>% 
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
am_all_BCPE_sigmapoly_param <- best_model_param(am_all_BCPE_sigmapoly, Ammonium)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
Ammonium_A <- Ammonium[Ammonium$period %in% c("A"), ]
Ammonium_A_pdf <- density(Ammonium_A$Ammonium)

# Period B
Ammonium_B <- Ammonium[Ammonium$period %in% c("B"), ]
Ammonium_B_pdf <- density(Ammonium_B$Ammonium)

# Period C
Ammonium_C <- Ammonium[Ammonium$period %in% c("C"), ]
Ammonium_C_pdf <- density(Ammonium_C$Ammonium)

# Period D
Ammonium_D <- Ammonium[Ammonium$period %in% c("D"), ]
Ammonium_D_pdf <- density(Ammonium_D$Ammonium)

# Period E
Ammonium_E <- Ammonium[Ammonium$period %in% c("E"), ]
Ammonium_E_pdf <- density(Ammonium_E$Ammonium)

# Period F
Ammonium_F <- Ammonium[Ammonium$period %in% c("F"), ]
Ammonium_F_pdf <- density(Ammonium_F$Ammonium)

# Period G
Ammonium_G <- Ammonium[Ammonium$period %in% c("G"), ]
Ammonium_G_pdf <- density(Ammonium_G$Ammonium)

# Period H
Ammonium_H <- Ammonium[Ammonium$period %in% c("H"), ]
Ammonium_H_pdf <- density(Ammonium_H$Ammonium)


# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(am_all_BCPE_sigmapoly$y, na.rm = TRUE),
         max(am_all_BCPE_sigmapoly$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
am_all_BCPE_sigmapoly_A_pdf  <- dBCPE(x, mu = 9.054708, sigma = 0.4029700, nu = 0.9245075, tau = 1.950330)
am_all_BCPE_sigmapoly_B_pdf  <- dBCPE(x, mu = 8.020892, sigma = 0.4488957, nu = 0.8023085, tau = 1.884630)
am_all_BCPE_sigmapoly_C_pdf  <- dBCPE(x, mu = 7.151599, sigma = 0.4798892, nu = 0.7293900, tau = 1.805774)
am_all_BCPE_sigmapoly_D_pdf  <- dBCPE(x, mu = 6.264545, sigma = 0.5123985, nu = 0.6269449, tau = 1.739688)
am_all_BCPE_sigmapoly_E_pdf  <- dBCPE(x, mu = 5.395271, sigma = 0.5328389, nu = 0.5652574, tau = 1.669930)
am_all_BCPE_sigmapoly_F_pdf  <- dBCPE(x, mu = 4.527159, sigma = 0.5450763, nu = 0.4770014, tau = 1.606992)
am_all_BCPE_sigmapoly_G_pdf  <- dBCPE(x, mu = 3.670040, sigma = 0.5459937, nu = 0.3921758, tau = 1.541246)
am_all_BCPE_sigmapoly_H_pdf  <- dBCPE(x, mu = 2.759726, sigma = 0.5463837, nu = 0.3031668, tau = 1.488819)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(Ammonium_A_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(am_all_BCPE_sigmapoly_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(Ammonium_B_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(am_all_BCPE_sigmapoly_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(Ammonium_C_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.13))
lines(am_all_BCPE_sigmapoly_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(Ammonium_D_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(am_all_BCPE_sigmapoly_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(Ammonium_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(am_all_BCPE_sigmapoly_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(Ammonium_F_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(am_all_BCPE_sigmapoly_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(Ammonium_G_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(am_all_BCPE_sigmapoly_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(Ammonium_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.3))
lines(am_all_BCPE_sigmapoly_H_pdf, type = "l", lwd = 2, col = "goldenrod") 



