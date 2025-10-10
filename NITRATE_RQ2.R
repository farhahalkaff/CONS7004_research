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
  select(Date, Nitrate, year, month)


## STEP 1 ==============================================================================================================================

# 1. Run four-parameter family distribution that is on a positive real line distribution --------------------------------------
#(BCT, BCTo, BCPE, BCPEo and GB2). run all param (time) 

#== BCT ==#
ni_all_BCT <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCT(), data = Nitrate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: right, BCT can't really run if values are close to 0, because the model turns mu negative during iteration, and 
## since BCT is on a positive real line distribution, and error pops up. 

#== BCTo ==#
ni_all_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCTo(), data = Nitrate,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## all good here

#== BCPE ==#
ni_all_BCPE <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPE(), data = Nitrate,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: mu must be positive 

#== BCPEo ==#
ni_all_BCPEo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = BCPEo(), data = Nitrate,
                       method = mixed(5,100),
                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got an error saying while (abs(G.dev.in - i.G.dev) > i.c.crit && i.iter < i.n.cyc) { : 
# missing value where TRUE/FALSE needed

#== GB2 ==#
ni_all_GB2 <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = GB2(), data = Nitrate,
                     mu.start = mean(Nitrate$AmmoNitratenium), sigma.start = sd(Nitrate$Nitrate),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
# NOTE: mu must be positive 


# 2. Run models with JSU or SST family distribution ----------------------------------------------------------------------------

#== SST ==#
ni_all_SST <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SST(), data = Nitrate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: sigma must be positive error. always get this error when running SST on some nutrient! not sure why


#== JSU ==#
ni_all_JSU <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = JSU(), data = Nitrate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here!


# 3. check which has the lowest AIC and use that family distribution -----------------------------------------------------------

AIC(ni_all_JSU)
AIC(ni_all_BCTo)


# the rankings are:
## 1. BCTo
## 2. JSU



# 4. Make sure that family can run all time-varying parameter models ------------------------------------------------------------

## MEAN ONLY MODEL ##
ni_meanonly_BCTo <- gamlss(Nitrate ~ year + month,
                           family = BCTo(), data = Nitrate,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
ni_mean_and_sigma_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month,
                                 family = BCTo(), data = Nitrate,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
ni_contskew_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                           family = BCTo(), data = Nitrate,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
ni_contkurt_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = BCTo(), data = Nitrate,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# all works. HOORAY! 


# 5b. IF YES: check AIC for time varying parameters -------------------------------------------------------------------------------

AIC(ni_meanonly_BCTo)
AIC(ni_mean_and_sigma_BCTo)
AIC(ni_contskew_BCTo)
AIC(ni_contkurt_BCTo)
AIC(ni_all_BCTo)

# rankings are:
## 1. all param (time): 43665.37
## 2. constant skewness: 43675.66
## 3. constant kurtosis: 43710.15
## 4. mean and sigma: 43727.12
## 5. mean only: 44396.56


# 6. Then check how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Nitrate$mu_hat_all_BCTo_ni <- predict(ni_all_BCTo, what = "mu", type = "response")

ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCTo_ni), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()



# 9. Check PIT histogram, worm plot and moment bucket of the model -----------------------------------------------------------------------------

# fitted param for (am_all_SSTtr) -----------
ni_mu_hat_all_BCTo    <- predict(ni_all_BCTo, "mu", type = "response")
ni_sigma_hat_all_BCTo <- predict(ni_all_BCTo, "sigma", type = "response")
ni_nu_hat_all_BCTo    <- predict(ni_all_BCTo, "nu", type = "response")
ni_tau_hat_all_BCTo   <- predict(ni_all_BCTo, "tau", type = "response")

ni_pit <- pBCTo(Nitrate$Nitrate, mu = ni_mu_hat_all_BCTo, sigma = ni_sigma_hat_all_BCTo, 
                nu = ni_nu_hat_all_BCTo, tau = ni_tau_hat_all_BCTo)

# Plot PIT histogram
hist(ni_pit, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(ni_pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity
## NOTE: yeah looks pretty good

resid_wp(ni_all_BCTo)
## NOTE: yeah absolutely not

moment_bucket(ni_all_BCTo)
## NOTE: platy!!!!!


# 7b. IF YES: then that's your model baby! ---------------------------------------------------------------------------------------

# now let's check the summary of the model
summary(am_all_BCPE_sigmapoly)

## ==================================================================================================================================


## STEP 2 ==============================================================================================================================

# 1. subset few years together to get period ---------------------------------------------------------------------------------------
# Add column to combine few years together 
Nitrate <- Nitrate %>% 
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
ni_all_BCTo_param <- best_model_param(ni_all_BCTo, Nitrate)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
Nitrate_A <- Nitrate[Nitrate$period %in% c("A"), ]
Nitrate_A_pdf <- density(Nitrate_A$Nitrate)

# Period B
Nitrate_B <- Nitrate[Nitrate$period %in% c("B"), ]
Nitrate_B_pdf <- density(Nitrate_B$Nitrate)

# Period C
Nitrate_C <- Nitrate[Nitrate$period %in% c("C"), ]
Nitrate_C_pdf <- density(Nitrate_C$Nitrate)

# Period D
Nitrate_D <- Nitrate[Nitrate$period %in% c("D"), ]
Nitrate_D_pdf <- density(Nitrate_D$Nitrate)

# Period E
Nitrate_E <- Nitrate[Nitrate$period %in% c("E"), ]
Nitrate_E_pdf <- density(Nitrate_E$Nitrate)

# Period F
Nitrate_F <- Nitrate[Nitrate$period %in% c("F"), ]
Nitrate_F_pdf <- density(Nitrate_F$Nitrate)

# Period G
Nitrate_G <- Nitrate[Nitrate$period %in% c("G"), ]
Nitrate_G_pdf <- density(Nitrate_G$Nitrate)

# Period H
Nitrate_H <- Nitrate[Nitrate$period %in% c("H"), ]
Nitrate_H_pdf <- density(Nitrate_H$Nitrate)


# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for ni_all_BCTo -----
x <- seq(min(ni_all_BCTo$y, na.rm = TRUE),
         max(ni_all_BCTo$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
ni_all_BCTo_A_pdf  <- dBCTo(x, mu = 5.409477, sigma = 0.8046698, nu = 0.2206034, tau = 1.891183e+07)
ni_all_BCTo_B_pdf  <- dBCTo(x, mu = 6.931207, sigma = 0.7584897, nu = 0.2131525, tau = 1.978359e+06)
ni_all_BCTo_C_pdf  <- dBCTo(x, mu = 8.756472, sigma = 0.7125375, nu = 0.2131189, tau = 2.671558e+05)
ni_all_BCTo_D_pdf  <- dBCTo(x, mu = 10.509681, sigma = 0.6874634, nu = 0.2157017, tau = 3.893915e+04)
ni_all_BCTo_E_pdf  <- dBCTo(x, mu = 13.194903, sigma = 0.6472785, nu = 0.2128342, tau = 5.296185e+03)
ni_all_BCTo_F_pdf  <- dBCTo(x, mu = 16.266895, sigma = 0.6154905, nu = 0.2139006, tau = 7.842388e+02)
ni_all_BCTo_G_pdf  <- dBCTo(x, mu = 19.944658, sigma = 0.5856530, nu = 0.2170114, tau = 1.092946e+02)
ni_all_BCTo_H_pdf  <- dBCTo(x, mu = 24.854995, sigma = 0.5560830, nu = 0.2118573, tau = 1.569397e+01)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(Nitrate_A_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.12))
lines(ni_all_BCTo_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(Nitrate_B_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.1))
lines(ni_all_BCTo_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(Nitrate_C_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ni_all_BCTo_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(Nitrate_D_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.065))
lines(ni_all_BCTo_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(Nitrate_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(ni_all_BCTo_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(Nitrate_F_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.045))
lines(ni_all_BCTo_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(Nitrate_G_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.04))
lines(ni_all_BCTo_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(Nitrate_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.032))
lines(ni_all_BCTo_H_pdf, type = "l", lwd = 2, col = "goldenrod") 
