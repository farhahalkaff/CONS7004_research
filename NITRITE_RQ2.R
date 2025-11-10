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
nii_all_BCPEo_allpoly_param <- best_model_param(nii_all_BCPEo_allpoly, Nitrite)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
Nitrite_A <- Nitrite[Nitrite$period %in% c("1"), ]
Nitrite_A_pdf <- density(Nitrite_A$Nitrite)

# Period B
Nitrite_B <- Nitrite[Nitrite$period %in% c("2"), ]
Nitrite_B_pdf <- density(Nitrite_B$Nitrite)

# Period C
Nitrite_C <- Nitrite[Nitrite$period %in% c("3"), ]
Nitrite_C_pdf <- density(Nitrite_C$Nitrite)

# Period D
Nitrite_D <- Nitrite[Nitrite$period %in% c("4"), ]
Nitrite_D_pdf <- density(Nitrite_D$Nitrite)

# Period E
Nitrite_E <- Nitrite[Nitrite$period %in% c("5"), ]
Nitrite_E_pdf <- density(Nitrite_E$Nitrite)

# Period F
Nitrite_F <- Nitrite[Nitrite$period %in% c("6"), ]
Nitrite_F_pdf <- density(Nitrite_F$Nitrite)

# Period G
Nitrite_G <- Nitrite[Nitrite$period %in% c("7"), ]
Nitrite_G_pdf <- density(Nitrite_G$Nitrite)

# Period H
Nitrite_H <- Nitrite[Nitrite$period %in% c("8"), ]
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





################

# METHOD 1 
logSurv(Nitrite$Nitrite, prob=0.90, tail="right")

# METHOD 2
par(mfrow = c(1, 3))
nii_m1 <- loglogSurv1(Nitrite$Nitrite, prob=0.90, title="(a) TYPE I")
nii_m2 <- loglogSurv2(Nitrite$Nitrite, prob=0.90, title="(b) TYPE II")
nii_m3 <- loglogSurv3(Nitrite$Nitrite, prob=0.90, title="(c) TYPE III")
par(mfrow = c(1, 1))



# 1. Total number of observations
total_obs <- length(Nitrite$Nitrite)
print(paste("Total observations:", total_obs))

# 2. Number of observations in the top 10%
num_top_10_percent <- ceiling(total_obs * 0.10)
print(paste("Number of observations in the top 10%:", num_top_10_percent))




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













## adding smoothing spline 

nii_all_BCPE_spline <- gamlss(Nitrite ~ cs(year) + month, sigma.fo = ~ cs(year) + month, nu.fo = ~ cs(year) + month, tau.fo = ~ cs(year) + month,
                                family = BCPE(), data = Nitrite,
                                mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite))
                                #method = mixed(5,100),
                                #control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


nii_all_BCPE_spline_param <- best_model_param(nii_all_BCPE_spline, Nitrite)
Nitrite$mu_hat_year <- predict(nii_all_BCPE_spline, what = "tau", type = "response")


mean_dplyr_filtered <- Nitrite %>%
  filter(year %in% c("1966")) %>%
  summarise(mean_value1 = mean(mu_hat_year))
print(mean_dplyr_filtered)

set.seed(123)
nii_spline_predict <- list(
  A = rBCPE(656, mu = 0.5531981, sigma = 0.9942300, nu = 1.0523044, tau = 5.746193),
  B = rBCPE(532, mu = 0.6605832, sigma = 0.6494657, nu = 0.3551003, tau = 2.101785),
  C = rBCPE(562, mu = 0.7260273, sigma = 1.0089550, nu = 1.1233813, tau = 3.132722),
  D = rBCPE(862, mu = 0.6154539, sigma = 1.6903863, nu = 1.9191222, tau = 1.468543),
  E = rBCPE(959, mu = 0.8598428, sigma = 0.6254443, nu = 0.6322192, tau = 2.734192),
  F = rBCPE(980, mu = 0.7792452, sigma = 0.5956562, nu = 0.3224346, tau = 3.898431),
  G = rBCPE(973, mu = 0.6595062, sigma = 0.6746916, nu = 0.5138789, tau = 5.569839),
  H = rBCPE(968, mu = 0.6165068, sigma = 0.6516583, nu = 0.6742924, tau = 8.234633)
)


nii_spline_predict_long <- stack(nii_spline_predict)
names(nii_spline_predict_long) <- c("pred", "period")


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Nitrite, aes(x = Nitrite, y = period), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = nii_spline_predict_long, aes(x = pred, y = period), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  






citation("stats")


##### use bins in the model 
Nitrite$period <- as.numeric(Nitrite$period)

nii_all_BCPE_bins_spline <- gamlss(Nitrite ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                              family = BCPE(), data = Nitrite,
                              mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                              method = mixed(5,100),
                              control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


nii_all_BCPE_bins_spline_noseason <- gamlss(Nitrite ~ cs(period), sigma.fo = ~ cs(period), nu.fo = ~ cs(period), tau.fo = ~ cs(period),
                                   family = BCPE(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


nii_all_BCPE_bins_spline_param <- best_model_param(nii_all_BCPE_bins_spline, Nitrite)
nii_all_BCPE_bins_noseason_spline_param <- best_model_param(nii_all_BCPE_bins_spline_noseason, Nitrite)


Nitrite$mu_hat_bins <- predict(nii_all_BCPE_bins_spline, what ="tau", type = "response")


## with seasonaility 
set.seed(123)
nii_all_BCPE_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5968993, sigma = 0.8015084, nu = 0.7175134, tau = 1.653869),
  "2" = rBCPE(532, mu = 0.6823589, sigma = 0.6679866, nu = 0.3935962, tau = 2.243619),
  "3" = rBCPE(562, mu = 0.7388949, sigma = 0.5896781, nu = 0.2804600, tau = 2.545870),
  "4" = rBCPE(862, mu = 0.7681278, sigma = 0.5483098, nu = 0.3590524, tau = 2.258531),
  "5" = rBCPE(959, mu = 0.8142381, sigma = 0.5509836, nu = 0.4133060, tau = 1.755169),
  "6" = rBCPE(980, mu = 0.7736620, sigma = 0.6087074, nu = 0.4236331, tau = 1.827649),
  "7" = rBCPE(973, mu = 0.6910492, sigma = 0.6751632, nu = 0.3825676, tau = 2.494553),
  "8" = rBCPE(968, mu = 0.5908073, sigma = 0.6275279, nu = 0.4280078, tau = 2.981172)
)

nii_all_BCPE_bins_spline_predict_long <- stack(nii_all_BCPE_bins_spline_predict)
names(nii_all_BCPE_bins_spline_predict_long) <- c("pred", "period")
nii_all_BCPE_bins_spline_predict_long$period <- as.numeric(nii_all_BCPE_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Nitrite, aes(x = Nitrite, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = nii_all_BCPE_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Nitrite (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  


## without seasonality 
set.seed(123)
nii_all_BCPE_bins_noseason_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5458773, sigma = 0.8063905, nu = 0.31925398, tau = 2.740700),
  "2" = rBCPE(532, mu = 0.6407777, sigma = 0.7736957, nu = 0.29019910, tau = 2.778094),
  "3" = rBCPE(562, mu = 0.7392532, sigma = 0.7671997, nu = 0.22006368, tau = 2.883818),
  "4" = rBCPE(862, mu = 0.7918942, sigma = 0.7403506, nu = 0.11735362, tau = 2.376556),
  "5" = rBCPE(959, mu = 0.7990549, sigma = 0.6684704, nu = 0.04978452, tau = 1.761647),
  "6" = rBCPE(980, mu = 0.6989036, sigma = 0.6784926, nu = 0.12317597, tau = 1.458834),
  "7" = rBCPE(973, mu = 0.5907477, sigma = 0.7762875, nu = 0.20933141, tau = 1.610333),
  "8" = rBCPE(968, mu = 0.4968001, sigma = 0.7911600, nu = 0.23693456, tau = 2.250268)
)


nii_all_BCPE_bins_noseason_spline_predict_long <- stack(nii_all_BCPE_bins_noseason_spline_predict)
names(nii_all_BCPE_bins_noseason_spline_predict_long) <- c("pred", "period")
nii_all_BCPE_bins_noseason_spline_predict_long $period <- as.numeric(nii_all_BCPE_bins_noseason_spline_predict_long $period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Nitrite, aes(x = Nitrite, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = nii_all_BCPE_bins_noseason_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Nitrite (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic() 









### bins, changing paramter time-varying 

nii_all_BCPE_bins_spline <- gamlss(Nitrite ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                   family = BCPE(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


nii_nokurt_BCPE_bins_spline <- gamlss(Nitrite ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ 1,
                                   family = BCPE(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

nii_noskew_BCPE_bins_spline <- gamlss(Nitrite ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ cs(period) + month,
                                   family = BCPE(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

nii_mean_sigma_BCPE_bins_spline <- gamlss(Nitrite ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ 1,
                                   family = BCPE(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


nii_meanonly_BCPE_bins_spline <- gamlss(Nitrite ~ cs(period) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                                   family = BCPE(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


AIC(nii_all_BCPE_bins_spline, nii_nokurt_BCPE_bins_spline,
    nii_noskew_BCPE_bins_spline, nii_mean_sigma_BCPE_bins_spline,
    nii_meanonly_BCPE_bins_spline)


##### RESULTS
# nii_all_BCPE_bins_spline         5045.469
# nii_noskew_BCPE_bins_spline      5117.305
# nii_nokurt_BCPE_bins_spline      5149.380
# nii_mean_sigma_BCPE_bins_spline  5183.075
# nii_meanonly_BCPE_bins_spline    5634.140



library(gamlss.tr)
gen.trun(family = "JSU", type = "left")

nii_intercept_JSU <- gamlss(Nitrite ~ 1,
                                        family = JSU(), data = Nitrite,
                                        mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                        method = mixed(5,100),
                                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

# truncate it and see?
nii_intercept_JSUtr <- gamlss(Nitrite ~ 1,
                            family = JSUtr(), data = Nitrite,
                            mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                            method = mixed(5,100),
                            control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))




moment_bucket(nii_intercept_JSU) +
resid_wp(nii_intercept_JSU) + theme(plot.title = element_blank())









### model for each bin 

nii_all_BCPE_bin_1 <- gamlss(Nitrite ~ month, sigma.fo = ~  month, nu.fo = ~ month, tau.fo = ~ month,
                            family = BCPE(), data = Nitrite_A,
                            mu.start = mean(Nitrite_A$Nitrite), sigma.start = sd(Nitrite_A$Nitrite),
                            #method = mixed(5,100),
                            control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


nii_all_BCPE_bin_2 <- gamlss(Nitrite ~ month, sigma.fo = ~  month, nu.fo = ~ month, tau.fo = ~ month,
                             family = BCPE(), data = Nitrite_B,
                             mu.start = mean(Nitrite_B$Nitrite), sigma.start = sd(Nitrite_B$Nitrite),
                             #method = mixed(5,100),
                             control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

nii_all_BCPE_bin_3 <- gamlss(Nitrite ~ month, sigma.fo = ~  month, nu.fo = ~ month, tau.fo = ~ month,
                             family = BCPE(), data = Nitrite_C,
                             mu.start = mean(Nitrite_C$Nitrite), sigma.start = sd(Nitrite_C$Nitrite),
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













# monte carlo simulation

nii_all_BCPE_bins_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5968993, sigma = 0.8015084, nu = 0.7175134, tau = 1.653869),
  "2" = rBCPE(532, mu = 0.6823589, sigma = 0.6679866, nu = 0.3935962, tau = 2.243619),
  "3" = rBCPE(562, mu = 0.7388949, sigma = 0.5896781, nu = 0.2804600, tau = 2.545870),
  "4" = rBCPE(862, mu = 0.7681278, sigma = 0.5483098, nu = 0.3590524, tau = 2.258531),
  "5" = rBCPE(959, mu = 0.8142381, sigma = 0.5509836, nu = 0.4133060, tau = 1.755169),
  "6" = rBCPE(980, mu = 0.7736620, sigma = 0.6087074, nu = 0.4236331, tau = 1.827649),
  "7" = rBCPE(973, mu = 0.6910492, sigma = 0.6751632, nu = 0.3825676, tau = 2.494553),
  "8" = rBCPE(968, mu = 0.5908073, sigma = 0.6275279, nu = 0.4280078, tau = 2.981172)
)

library(gamlss.dist)
library(moments)

############################## bin 1
mu_nii_bin1   <-     0.5968993
sigma_nii_bin1 <-    0.8015084
nu_nii_bin1  <-    0.7175134
tau_nii_bin1  <- 1.653869

set.seed(1)
nii_bin1_sim <- rBCPE(N, mu = mu_nii_bin1, sigma = sigma_nii_bin1, nu = nu_nii_bin1, tau = tau_nii_bin1)

# get skewness and kurtosis 
kurtosis(nii_bin1_sim) - 3 # 1.371239
skewness(nii_bin1_sim) # 0.9752737
var(nii_bin1_sim) # 0.2010821
mean(nii_bin1_sim) # 0.6871321


############################# bin 2
mu_nii_bin2   <-     0.6823589
sigma_nii_bin2 <-    0.6679866
nu_nii_bin2  <-    0.3935962
tau_nii_bin2  <-  2.243619
  
set.seed(2)
nii_bin2_sim <- rBCPE(N, mu = mu_nii_bin2, sigma = sigma_nii_bin2, nu = nu_nii_bin2, tau = tau_nii_bin2)

# get skewness and kurtosis 
kurtosis(nii_bin2_sim) - 3 # 1.342921
skewness(nii_bin2_sim) # 1.036183
var(nii_bin2_sim) # 0.2349244
mean(nii_bin2_sim) # 0.7753099


############################# bin 3
mu_nii_bin3   <-   0.7388949  
sigma_nii_bin3 <-    0.5896781
nu_nii_bin3  <-    0.2804600
tau_nii_bin3  <- 2.545870
  
set.seed(3)
nii_bin3_sim <- rBCPE(N, mu = mu_nii_bin3, sigma = sigma_nii_bin3, nu = nu_nii_bin3, tau = tau_nii_bin3)

# get skewness and kurtosis 
kurtosis(nii_bin3_sim) - 3 # 1.075109
skewness(nii_bin3_sim) # 0.9949351
var(nii_bin3_sim) # 0.2220873
mean(nii_bin3_sim) # 0.8314637



############################# bin 4
mu_nii_bin4   <-   0.7681278  
sigma_nii_bin4 <-    0.5483098
nu_nii_bin4  <-    0.3590524
tau_nii_bin4  <- 2.258531
  
set.seed(4)
nii_bin4_sim <- rBCPE(N, mu = mu_nii_bin4, sigma = sigma_nii_bin4, nu = nu_nii_bin4, tau = tau_nii_bin4)

# get skewness and kurtosis 
kurtosis(nii_bin4_sim) - 3 # 0.9965382
skewness(nii_bin4_sim) # 0.908888
var(nii_bin4_sim) # 0.1961895
mean(nii_bin4_sim) # 0.8416406


############################# bin 5
mu_nii_bin5   <-     0.8142381
sigma_nii_bin5 <-    0.5509836
nu_nii_bin5  <-    0.4133060
tau_nii_bin5  <- 1.755169
  
set.seed(5)
nii_bin5_sim <- rBCPE(N, mu = mu_nii_bin5, sigma = sigma_nii_bin5, nu = nu_nii_bin5, tau = tau_nii_bin5)

# get skewness and kurtosis 
kurtosis(nii_bin5_sim) - 3 # 2.100385
skewness(nii_bin5_sim) # 1.083375
var(nii_bin5_sim) # 0.2216022
mean(nii_bin5_sim) # 0.8864044


############################# bin 6
mu_nii_bin6   <-     0.7736620
sigma_nii_bin6 <-    0.6087074
nu_nii_bin6  <-    0.4236331
tau_nii_bin6  <- 1.827649
  
set.seed(6)
nii_bin6_sim <- rBCPE(N, mu = mu_nii_bin6, sigma = sigma_nii_bin6, nu = nu_nii_bin6, tau = tau_nii_bin6)

# get skewness and kurtosis 
kurtosis(nii_bin6_sim) - 3 # 2.161513
skewness(nii_bin6_sim) # 1.130199
var(nii_bin6_sim) # 0.2456636
mean(nii_bin6_sim) # 0.8573267


############################# bin 7
mu_nii_bin7   <-    0.6910492 
sigma_nii_bin7 <-    0.6751632
nu_nii_bin7  <-    0.3825676
tau_nii_bin7  <- 2.494553
  
set.seed(7)
nii_bin7_sim <- rBCPE(N, mu = mu_nii_bin7, sigma = sigma_nii_bin7, nu = nu_nii_bin7, tau = tau_nii_bin7)

# get skewness and kurtosis 
kurtosis(nii_bin7_sim) - 3 # 0.9394921
skewness(nii_bin7_sim) # 0.9698492
var(nii_bin7_sim) # 0.2458738
mean(nii_bin7_sim) # 0.7883491


############################# bin 8
mu_nii_bin8   <-     0.5908073
sigma_nii_bin8 <-    0.6275279
nu_nii_bin8  <-    0.4280078
tau_nii_bin8  <- 2.981172
  
set.seed(8)
nii_bin8_sim <- rBCPE(N, mu = mu_nii_bin8, sigma = sigma_nii_bin8, nu = nu_nii_bin8, tau = tau_nii_bin8)

# get skewness and kurtosis 
kurtosis(nii_bin8_sim) - 3 # 0.1050014
skewness(nii_bin8_sim) # 0.7317325
var(nii_bin8_sim) # 0.1476112
mean(nii_bin8_sim) # 0.65703



#### without seasonaility =======================================================================

nii_all_BCPE_bins_noseason_spline_predict <- list(
  "1" = rBCPE(656, mu = 0.5458773, sigma = 0.8063905, nu = 0.31925398, tau = 2.740700),
  "2" = rBCPE(532, mu = 0.6407777, sigma = 0.7736957, nu = 0.29019910, tau = 2.778094),
  "3" = rBCPE(562, mu = 0.7392532, sigma = 0.7671997, nu = 0.22006368, tau = 2.883818),
  "4" = rBCPE(862, mu = 0.7918942, sigma = 0.7403506, nu = 0.11735362, tau = 2.376556),
  "5" = rBCPE(959, mu = 0.7990549, sigma = 0.6684704, nu = 0.04978452, tau = 1.761647),
  "6" = rBCPE(980, mu = 0.6989036, sigma = 0.6784926, nu = 0.12317597, tau = 1.458834),
  "7" = rBCPE(973, mu = 0.5907477, sigma = 0.7762875, nu = 0.20933141, tau = 1.610333),
  "8" = rBCPE(968, mu = 0.4968001, sigma = 0.7911600, nu = 0.23693456, tau = 2.250268)
)


############################## bin 1
mu_nii_bin1_noseason   <-     0.5458773
sigma_nii_bin1_noseason <-    0.8063905
nu_nii_bin1_noseason  <-    0.31925398
tau_nii_bin1_noseason  <- 2.740700

set.seed(11)
nii_bin1_noseason_sim <- rBCPE(N, mu = mu_nii_bin1_noseason, sigma = sigma_nii_bin1_noseason,
                               nu = nu_nii_bin1_noseason, tau = tau_nii_bin1_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin1_noseason_sim) - 3 # 1.552272
skewness(nii_bin1_noseason_sim) # 1.188707
var(nii_bin1_noseason_sim) # 0.2445721
mean(nii_bin1_noseason_sim) # 0.6680511


############################# bin 2
mu_nii_bin2_noseason   <-     0.6407777
sigma_nii_bin2_noseason <-    0.7736957
nu_nii_bin2_noseason  <-    0.29019910
tau_nii_bin2_noseason  <- 2.778094

set.seed(22)
nii_bin2_noseason_sim <- rBCPE(N, mu = mu_nii_bin2_noseason, sigma = sigma_nii_bin2_noseason,
                               nu = nu_nii_bin2_noseason, tau = tau_nii_bin2_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin2_noseason_sim) - 3 # 1.57843
skewness(nii_bin2_noseason_sim) # 1.192064
var(nii_bin2_noseason_sim) # 0.3133341
mean(nii_bin2_noseason_sim) # 0.7772897



############################# bin 3
mu_nii_bin3_noseason   <-     0.7392532
sigma_nii_bin3_noseason <-    0.7671997
nu_nii_bin3_noseason  <-    0.22006368
tau_nii_bin3_noseason  <- 2.883818

set.seed(33)
nii_bin3_noseason_sim <- rBCPE(N, mu = mu_nii_bin3_noseason, sigma = sigma_nii_bin3_noseason,
                               nu = nu_nii_bin3_noseason, tau = tau_nii_bin3_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin3_noseason_sim) - 3 # 1.953887
skewness(nii_bin3_noseason_sim) # 1.288705
var(nii_bin3_noseason_sim) # 0.4397112
mean(nii_bin3_noseason_sim) # 0.9121842



############################# bin 4
mu_nii_bin4_noseason   <-     0.7918942
sigma_nii_bin4_noseason <-    0.7403506
nu_nii_bin4_noseason  <-    0.11735362
tau_nii_bin4_noseason  <- 2.376556

set.seed(44)
nii_bin4_noseason_sim <- rBCPE(N, mu = mu_nii_bin4_noseason, sigma = sigma_nii_bin4_noseason,
                               nu = nu_nii_bin4_noseason, tau = tau_nii_bin4_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin4_noseason_sim) - 3 # 5.037747
skewness(nii_bin4_noseason_sim) # 1.785265
var(nii_bin4_noseason_sim) # 0.5561944
mean(nii_bin4_noseason_sim) # 0.9954447



############################# bin 5
mu_nii_bin5_noseason   <-     0.7990549
sigma_nii_bin5_noseason <-    0.6684704
nu_nii_bin5_noseason  <-    0.04978452
tau_nii_bin5_noseason  <- 1.761647

set.seed(55)
nii_bin5_noseason_sim <- rBCPE(N, mu = mu_nii_bin5_noseason, sigma = sigma_nii_bin5_noseason,
                               nu = nu_nii_bin5_noseason, tau = tau_nii_bin5_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin5_noseason_sim) - 3 # 18.50362
skewness(nii_bin5_noseason_sim) # 2.880271
var(nii_bin5_noseason_sim) # 0.5385658
mean(nii_bin5_noseason_sim) # 0.9839833



############################# bin 6
mu_nii_bin6_noseason   <-     0.6989036
sigma_nii_bin6_noseason <-    0.6784926
nu_nii_bin6_noseason  <-    0.12317597
tau_nii_bin6_noseason  <- 1.458834

set.seed(66)
nii_bin6_noseason_sim <- rBCPE(N, mu = mu_nii_bin6_noseason, sigma = sigma_nii_bin6_noseason,
                               nu = nu_nii_bin6_noseason, tau = tau_nii_bin6_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin6_noseason_sim) - 3 # 27.06774
skewness(nii_bin6_noseason_sim) # 3.362734
var(nii_bin6_noseason_sim) # 0.4174064
mean(nii_bin6_noseason_sim) # 0.8509449



############################# bin 7
mu_nii_bin7_noseason   <-     0.5907477
sigma_nii_bin7_noseason <-    0.7762875
nu_nii_bin7_noseason  <-    0.20933141
tau_nii_bin7_noseason  <- 1.610333

set.seed(77)
nii_bin7_noseason_sim <- rBCPE(N, mu = mu_nii_bin7_noseason, sigma = sigma_nii_bin7_noseason,
                               nu = nu_nii_bin7_noseason, tau = tau_nii_bin7_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin7_noseason_sim) - 3 # 14.54437
skewness(nii_bin7_noseason_sim) # 2.622423
var(nii_bin7_noseason_sim) # 0.3517711
mean(nii_bin7_noseason_sim) # 0.7376265



############################# bin 8
mu_nii_bin8_noseason   <-     0.4968001
sigma_nii_bin8_noseason <-    0.7911600
nu_nii_bin8_noseason  <-    0.23693456
tau_nii_bin8_noseason  <- 2.250268

set.seed(88)
nii_bin8_noseason_sim <- rBCPE(N, mu = mu_nii_bin8_noseason, sigma = sigma_nii_bin8_noseason,
                               nu = nu_nii_bin8_noseason, tau = tau_nii_bin8_noseason)

# get skewness and kurtosis 
kurtosis(nii_bin8_noseason_sim) - 3 # 3.91345
skewness(nii_bin8_noseason_sim) # 1.604267
var(nii_bin8_noseason_sim) # 0.2213472
mean(nii_bin8_noseason_sim) # 0.6188528















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

# get them errors gurll
nii_bin1_sim_se <- bootstrap_moments_sim_only(nii_bin1_sim) 
nii_bin2_sim_se <- bootstrap_moments_sim_only(nii_bin2_sim) 
nii_bin3_sim_se <- bootstrap_moments_sim_only(nii_bin3_sim) 
nii_bin4_sim_se <- bootstrap_moments_sim_only(nii_bin4_sim) 
nii_bin5_sim_se <- bootstrap_moments_sim_only(nii_bin5_sim) 
nii_bin6_sim_se <- bootstrap_moments_sim_only(nii_bin6_sim) 
nii_bin7_sim_se <- bootstrap_moments_sim_only(nii_bin7_sim) 
nii_bin8_sim_se <- bootstrap_moments_sim_only(nii_bin8_sim) 

nii_bin9_sim_se$se

nii_bins_df <- tibble::tibble(
  year_mid = c(1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  mean = c(0.687,0.775,0.831,0.842,0.886,0.857,0.788,0.657),
  variance = c(0.201,0.235,0.222,0.196,0.222,0.246,0.246,0.148),
  skewness = c(0.975,1.036,0.995,0.909,1.084,1.130,0.970,0.732),
  kurtosis = c(1.371,1.342,1.075,0.997,2.100,2.161,0.939,0.105)
)

nii_se_df <- tibble::tibble(
  year_mid = nii_bins_df$year_mid,
  se_mean = c(0.0006448352,0.0006854193,0.0006687207,0.0006316300,
              0.0006792097,0.0007030647,0.0007034426,0.0005531132),
  se_variance = c(0.0005252876,0.0006202958,0.0005477940,0.0004864865,
                  0.0006394150,0.0007024575,0.0006084682,0.0003040004),
  se_skew = c(0.0053590092,0.0052255988,0.0046475068,0.0049345745,
              0.0073008419,0.0070670580,0.0043976990,0.0031003094),
  se_kurt = c(0.0285856634,0.0301624652,0.0236317279,0.0273073230,
              0.0500960830,0.0477571198,0.0206259617,0.0103299206)
)

# join and pivot to long form
nii_df_long <- nii_bins_df %>%
  left_join(nii_se_df, by = "year_mid") %>%
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


ggplot(nii_df_long, aes(x = year_mid, y = value)) +
  geom_line(size = 0.6) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.6) +   # vertical error bars
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "Year bins", y = NULL, title = "Distribution moments through time") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, hjust = 0.02)
  )









