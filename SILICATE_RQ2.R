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
  select(Date, Nitrite, year, month)


## STEP 1 ==============================================================================================================================

# 1. Run four-parameter family distribution that is on a positive real line distribution --------------------------------------
#(BCT, BCTo, BCPE, BCPEo and GB2). run all param (time) 

#== BCT ==#
si_all_BCT <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = BCT(), data = Silicate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got the error of: response variable out of range

#== BCTo ==#
si_all_BCTo <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCTo(), data = Silicate,
                      mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got the error of: response variable out of range

#== BCPE ==#
si_all_BCPE <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPE(), data = Silicate,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got the error of: response variable out of range

#== BCPEo ==#
si_all_BCPEo <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                       family = BCPEo(), data = Silicate,
                       method = mixed(5,100),
                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: got the error of: response variable out of range
 

#== GB2 ==#
si_all_GB2 <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = GB2(), data = Silicate,
                     #mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: NA's in the working vector or weights for parameter sigma


# 2. Run models with JSU or SST family distribution ----------------------------------------------------------------------------

#== SST ==#
si_all_SST <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = SST(), data = Silicate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: sigma must be positive error. always get this error when running SST on some nutrient! not sure why


#== JSU ==#
si_all_JSU <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = JSU(), data = Silicate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# all good here!


# 4. Make sure that family can run all time-varying parameter models ------------------------------------------------------------

## MEAN ONLY MODEL ##
si_meanonly_JSU <- gamlss(Silicate ~ year + month,
                           family = JSU(), data = Silicate,
                           method = mixed(5,100),
                           control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## MEAN + SIGMA MODEL ##
si_mean_and_sigma_JSU <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month,
                                 family = JSU(), data = Silicate,
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

## CONSTANT SKEWNESS ##
si_contskew_JSU <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month,
                           family = JSU(), data = Silicate,
                          mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                           method = mixed(5,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

## CONSTANT KURTOSIS ##
si_contkurt_JSU <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month,
                           family = JSU(), data = Silicate,
                          mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                           method = mixed(5,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))

# all works. HOORAY! 


# 5b. IF YES: check AIC for time varying parameters -------------------------------------------------------------------------------

AIC(si_meanonly_JSU)
AIC(si_mean_and_sigma_JSU)
AIC(si_contskew_JSU)
AIC(si_contkurt_JSU)
AIC(si_all_JSU)

# rankings are:
## 1. all param (time): 30996.57
## 2. constant kurtosis: 31246.31
## 3. constant skewness: 31676.63
## 4. mean and sigma: 32083.27
## 5. mean only: 33201.34


# 6. Then check how the model predicts the mean and plot them to see if it make sense ----------------------------------------------

Silicate$mu_hat_all_JSU_si <- predict(si_all_JSU, what = "mu", type = "response")

ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_JSU_si), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()


# 8. EXTRA: try adding poly() to each parameter to the model and see how that does

##== poly to the mean ==## 
si_all_JSU_meanpoly <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                               family = JSU(), data = Silicate,
                               method = mixed(5,100),
                               control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: Tis is better!: 30895.63


##== poly to the sigma ==## 
si_all_JSU_sigmapoly <- gamlss(Silicate ~ year + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                family = JSU(), data = Silicate,
                               mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                method = mixed(5,200),
                                control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: AIC lowered ever so slightly: 30920.64


##== poly to the skewness ==## 
si_all_JSU_skewpoly <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ poly(year,2) + month, tau.fo = ~ year + month,
                               family = JSU(), data = Silicate,
                              mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: does terribly! : 31050.65


##== poly to the kurtosis ==## 
si_all_JSU_kurtpoly <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ poly(year,2) + month,
                               family = JSU(), data = Silicate,
                              mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                               method = mixed(5,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: does really terrible: 31008.95

## NOTE: since adding poly to mean and sigma does well, let's add both into one model 

##== poly to the mean and sigma ==## 
si_all_JSU_meansigmapoly <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ poly(year,2) + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                              family = JSU(), data = Silicate,
                              mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                              method = mixed(5,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: convergence issue

#######################################################
## CONCLUSION: adding polynomial to sigma is better! ##
#######################################################


# 6. RECHECK how the model predicts the mean and plot them to see if it make sense ----------------------------------------------


Silicate$mu_hat_all_JSU_meanpoly_si <- predict(si_all_JSU_meanpoly, what = "mu", type = "response")

ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_JSU_meanpoly_si), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()

## Looks the same COOL!



# 9. Check PIT histogram, worm plot and moment bucket of the model -----------------------------------------------------------------------------

# fitted param for (am_all_SSTtr) -----------
si_mu_hat_all_JSU_meanpoly    <- predict(si_all_JSU_meanpoly, "mu", type = "response")
si_sigma_hat_all_JSU_meanpoly <- predict(si_all_JSU_meanpoly, "sigma", type = "response")
si_nu_hat_all_JSU_meanpoly    <- predict(si_all_JSU_meanpoly, "nu", type = "response")
si_tau_hat_all_JSU_meanpoly   <- predict(si_all_JSU_meanpoly, "tau", type = "response")

si_pit <- pJSU(Silicate$Silicate, mu = si_mu_hat_all_JSU_meanpoly, sigma = si_sigma_hat_all_JSU_meanpoly, 
                nu = si_nu_hat_all_JSU_meanpoly, tau = si_tau_hat_all_JSU_meanpoly)

# Plot PIT histogram
hist(si_pit, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(si_pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity
## NOTE: not too great

resid_wp(si_all_JSU_meanpoly)
## NOTE: tails ain't looking good

moment_bucket(si_all_JSU_meanpoly)
## NOTE: platy not considered, which makes sense, because JSU does not allow platy to occur, only
## either normal tails or lepto


## 10. Lets use SEP2 instead, since it's the best model for static distribution ----------------------------------------------

si_all_SEP2 <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                   family = SEP2(), data = Silicate,
                                   #mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: got error: missing value where TRUE/FALSE needed

si_all_SEP3 <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = SEP3(), data = Silicate,
                      mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                      method = mixed(5,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: convergence issue. Let's instead try other flexible family distribution like SHASHo


si_all_SHASHo2 <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = SHASHo2(), data = Silicate,
                      mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                      method = mixed(5,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
# NOTE: The AIC is larger than JSU, but das ok let's see how it does: 31582.36

# let's add poly first 
si_all_SHASHo2_meanpoly <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                         family = SHASHo2(), data = Silicate,
                         mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                         method = mixed(5,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: does really well! : 30743.17, better than JSU


## 11. see how well it predict the mean but looking at the plot -----------------------------------------------------------------

Silicate$mu_hat_all_SHASHo2_meanpoly_si <- predict(si_all_SHASHo2_meanpoly, what = "mu", type = "response")

ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_SHASHo2_meanpoly_si), color = "goldenrod", linewidth = 1) +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()

## HMM, it seem SHASHo2 really likes to predict negative values. So let's NOT use this 


## 12. let's check other flexible family distribution  ---------------------------------------------------------------------------

# SHASH 
si_all_SHASH_meanpoly <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                  family = SHASH(), data = Silicate,
                                  mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                  method = mixed(5,300),
                                  control = gamlss.control(n.cyc = 300, c.crit = 0.01, trace = TRUE))
# NOTE: hmm convergence issue


##### NOTE: ok, so I found a zero inflated/ zero-adjusted package from the book, so let's use that 
library(gamlss.inf)
## BCT
gen.Family(family="BCT", type="logit")
gen.Family(family="BCPE", type="logit")
gen.Family(family="SST", type="logit")
gen.Family(family="JSU", type="logit")




#== BCT ==#
si_all_BCTadj <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                     family = logitBCT(), data = Silicate,
                     method = mixed(5,100),
                     control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: still the same error of: response variable out of range


#== BCPE ==#
si_all_BCPEadj <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = logitBCPE(), data = Silicate,
                        method = mixed(5,100),
                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: still the same error of: response variable out of range


## OK now let's try real line distribution --------------------------------------------------------------

#== SST ==#
si_all_SSTadj <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = logitSST(), data = Silicate,
                        method = mixed(5,100),
                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: still the same error of: response variable out of range

#== JSU ==#
si_all_JSUadj <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = logitJSU(), data = Silicate,
                        method = mixed(5,100),
                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
## NOTE: still the same error of: response variable out of range


## left-truncate --------------------------------------------------------------

# really don't want to left truncate it, but right now don't have much of a choice. so let's try with SHASHo2

library(gamlss.tr)
gen.trun(0,"SHASHo2",type="left")

si_all_SHASHo2tr <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = SHASHo2tr(), data = Silicate,
                        mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                        method = mixed(5,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))
## NOTE: ERROR!: missing value where TRUE/FALSE needed


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



