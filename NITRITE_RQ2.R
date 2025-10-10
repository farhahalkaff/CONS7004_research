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
nii_all_BCT <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCT(), data = Nitrite,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: right, BCT can't really run if values are close to 0, because the model turns mu negative during iteration, and 
## since BCT is on a positive real line distribution, and error pops up. 

#== BCTo ==#
nii_all_BCTo <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCTo(), data = Nitrite,
                      mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got an error saying while (abs(G.dev.in - i.G.dev) > i.c.crit && i.iter < i.n.cyc) { : 
# missing value where TRUE/FALSE needed

#== BCPE ==#
nii_all_BCPE <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPE(), data = Nitrite,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== BCPEo ==#
nii_all_BCPEo <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = BCPEo(), data = Nitrite,
                       method = mixed(5,100),
                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here! 

#== GB2 ==#
nii_all_GB2 <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = GB2(), data = Nitrite,
                     mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
# convergence issue! 


# 2. Run models with JSU or SST family distribution ----------------------------------------------------------------------------

#== SST ==#
nii_all_SST <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SST(), data = Nitrite,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: sigma must be positive error. always get this error when running SST on some nutrient! not sure why


#== JSU ==#
nii_all_JSU <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = JSU(), data = Nitrite,
                     mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
# all good here!


# 3. check which has the lowest AIC and use that family distribution -----------------------------------------------------------

AIC(nii_all_BCPE)
AIC(nii_all_BCPEo)
AIC(nii_all_JSU)


# the rankings are:
## 1. BCPEo: 5679.88
## 2. BCPE: 5712.723
## 3. JSU: 5946.483


# 4. Make sure that family can run all time-varying parameter models ------------------------------------------------------------

## MEAN ONLY MODEL ##
nii_meanonly_BCPEo <- gamlss(Nitrite ~ year + month,
                           family = BCPEo(), data = Nitrite,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
nii_mean_and_sigma_BCPEo <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month,
                                 family = BCPEo(), data = Nitrite,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
nii_contskew_BCPEo <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                           family = BCPEo(), data = Nitrite,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
nii_contkurt_BCPEo <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = BCPEo(), data = Nitrite,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# all works. HOORAY! 


# 5b. IF YES: check AIC for time varying parameters -------------------------------------------------------------------------------

AIC(nii_meanonly_BCPEo)
AIC(nii_mean_and_sigma_BCPEo)
AIC(nii_contskew_BCPEo)
AIC(nii_contkurt_BCPEo)
AIC(nii_all_BCPEo)

# rankings are:
## 1. all param (time): 5679.88
## 2. constant skewness: 5711.564
## 3. constant kurtosis: 5742.003
## 4. mean and sigma: 5745.846
## 5. mean only: 6032.581


## NOTE: constant skewness has lower AIC than constant kurtosis, means skewness not really changing
## through time much? 


# 6. Then check how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Nitrite$mu_hat_all_BCPEo_nii <- predict(nii_all_BCPEo, what = "mu", type = "response")

ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPEo_nii), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()


# 8. EXTRA: try adding poly() to each parameter to the model and see how that does

##== poly to the mean ==## 
nii_all_BCPEo_meanpoly <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = BCPEo(), data = Nitrite,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: AIC is much better!: 5153.492


##== poly to the sigma ==## 
nii_all_BCPEo_sigmapoly <- gamlss(Nitrite ~ year + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                family = BCPEo(), data = Nitrite,
                                method = mixed(5,100),
                                control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: AIC is a little lower! : 5671.306


##== poly to the skewness ==## 
nii_all_BCPEo_skewpoly <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ year + month,
                               family = BCPEo(), data = Nitrite,
                               mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: doesn't do too well either: 5670.382


##== poly to the kurtosis ==## 
nii_all_BCPEo_kurtpoly <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ poly(year,2) + month,
                               family = BCPEo(), data = Nitrite,
                               mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: AIC is a little lower: 5686.603


##== add poly to mean and sigma ==#
nii_all_BCPEo_meansigmapoly <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                 family = BCPEo(), data = Nitrite,
                                 mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: SO GOOD: 5082.55

##== add poly to mean and sigma and skew and kurt ==#
nii_all_BCPEo_allpoly <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ poly(year,2) + month,
                                      family = BCPEo(), data = Nitrite,
                                      mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                      method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: THE BEST: 5072.28


# 6. RECHECK how the model predicts the mean and plot them to see if it make sense ----------------------------------------------


Nitrite$mu_hat_all_BCPEo_allpoly <- predict(nii_all_BCPEo_allpoly, what = "mu", type = "response")

ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCPEo_allpoly), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()

## Looks noice! COOL!



# 9. Check PIT histogram, worm plot and moment bucket of the model -----------------------------------------------------------------------------

# fitted param for (am_all_SSTtr) -----------
nii_mu_hat_all_BCPEo_allpoly    <- predict(nii_all_BCPEo_allpoly, "mu", type = "response")
nii_sigma_hat_all_BCPEo_allpoly <- predict(nii_all_BCPEo_allpoly, "sigma", type = "response")
nii_nu_hat_all_BCPEo_allpoly    <- predict(nii_all_BCPEo_allpoly, "nu", type = "response")
nii_tau_hat_all_BCPEo_allpoly   <- predict(nii_all_BCPEo_allpoly, "tau", type = "response")

nii_pit <- pBCPEo(Nitrite$Nitrite, mu = nii_mu_hat_all_BCPEo_allpoly, sigma = nii_sigma_hat_all_BCPEo_allpoly, 
                nu = nii_nu_hat_all_BCPEo_allpoly, tau = nii_tau_hat_all_BCPEo_allpoly)

# Plot PIT histogram
hist(nii_pit, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(nii_pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity
## NOTE: yeah looks pretty good

resid_wp(nii_all_BCPEo_allpoly)
## NOTE: looks alright, but taily

moment_bucket(nii_all_BCPEo_allpoly)
## NOTE: looks alright, a little platy


# 7b. IF YES: then that's your model baby! ---------------------------------------------------------------------------------------

# now let's check the summary of the model
summary(nii_all_BCPEo_allpoly)

## ==================================================================================================================================


## STEP 2 ==============================================================================================================================

# 1. subset few years together to get period ---------------------------------------------------------------------------------------
# Add column to combine few years together 
Nitrite <- Nitrite %>% 
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
nii_all_BCPEo_allpoly_param <- best_model_param(nii_all_BCPEo_allpoly, Nitrite)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
Nitrite_A <- Nitrite[Nitrite$period %in% c("A"), ]
Nitrite_A_pdf <- density(Nitrite_A$Nitrite)

# Period B
Nitrite_B <- Nitrite[Nitrite$period %in% c("B"), ]
Nitrite_B_pdf <- density(Nitrite_B$Nitrite)

# Period C
Nitrite_C <- Nitrite[Nitrite$period %in% c("C"), ]
Nitrite_C_pdf <- density(Nitrite_C$Nitrite)

# Period D
Nitrite_D <- Nitrite[Nitrite$period %in% c("D"), ]
Nitrite_D_pdf <- density(Nitrite_D$Nitrite)

# Period E
Nitrite_E <- Nitrite[Nitrite$period %in% c("E"), ]
Nitrite_E_pdf <- density(Nitrite_E$Nitrite)

# Period F
Nitrite_F <- Nitrite[Nitrite$period %in% c("F"), ]
Nitrite_F_pdf <- density(Nitrite_F$Nitrite)

# Period G
Nitrite_G <- Nitrite[Nitrite$period %in% c("G"), ]
Nitrite_G_pdf <- density(Nitrite_G$Nitrite)

# Period H
Nitrite_H <- Nitrite[Nitrite$period %in% c("H"), ]
Nitrite_H_pdf <- density(Nitrite_H$Nitrite)


# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(nii_all_BCPEo_allpoly$y, na.rm = TRUE),
         max(nii_all_BCPEo_allpoly$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
nii_all_BCPEo_allpoly_A_pdf  <- dBCPEo(x, mu = 0.5812946, sigma = 0.7301296, nu = 0.6030427, tau = 1.828145)
nii_all_BCPEo_allpoly_B_pdf  <- dBCPEo(x, mu = 0.7264105, sigma = 0.6104033, nu = 0.4488233, tau = 1.936335)
nii_all_BCPEo_allpoly_C_pdf  <- dBCPEo(x, mu = 0.8283477, sigma = 0.5481248, nu = 0.3555217, tau = 2.013209)
nii_all_BCPEo_allpoly_D_pdf  <- dBCPEo(x, mu = 0.8722824, sigma = 0.5228625, nu = 0.3107930, tau = 2.063076)
nii_all_BCPEo_allpoly_E_pdf  <- dBCPEo(x, mu = 0.8640284, sigma = 0.5134612, nu = 0.2980582, tau = 2.127777)
nii_all_BCPEo_allpoly_F_pdf  <- dBCPEo(x, mu = 0.7971493, sigma = 0.5312123, nu = 0.3240426, tau = 2.170300)
nii_all_BCPEo_allpoly_G_pdf  <- dBCPEo(x, mu = 0.6878732, sigma = 0.5758329, nu = 0.3938743, tau = 2.190660)
nii_all_BCPEo_allpoly_H_pdf  <- dBCPEo(x, mu = 0.5469675, sigma = 0.6501158, nu = 0.5088351, tau = 2.164449)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(Nitrite_A_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(Nitrite_B_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(Nitrite_C_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(Nitrite_D_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(Nitrite_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(Nitrite_F_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(Nitrite_G_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(Nitrite_H_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(nii_all_BCPEo_allpoly_H_pdf, type = "l", lwd = 2, col = "goldenrod") 



plot(nii_all_BCPEo_allpoly_A_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_B_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_C_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_D_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_E_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_F_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_G_pdf, type = "l")
plot(nii_all_BCPEo_allpoly_H_pdf, type = "l")




