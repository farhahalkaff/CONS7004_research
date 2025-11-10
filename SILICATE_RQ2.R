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
Silicate$mu_hat_all_SHASHo2_si <- predict(si_all_SHASHo2, what = "mu", type = "response")

ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_SHASHo2_meanpoly_si), color = "goldenrod", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_SHASHo2_si), color = "steelblue", linewidth = 1) +
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
Silicate <- Silicate %>% 
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
si_all_JSU_meanpoly_param <- best_model_param(si_all_JSU_meanpoly, Silicate)


# 3. get non-parametric estimate of pdf ----------------------------------------------------------------------------------------

# Period A
Silicate_A <- Silicate[Silicate$period %in% c("1"), ]
Silicate_A_pdf <- density(Silicate_A$Silicate)

# Period B
Silicate_B <- Silicate[Silicate$period %in% c("2"), ]
Silicate_B_pdf <- density(Silicate_B$Silicate)

# Period C
Silicate_C <- Silicate[Silicate$period %in% c("3"), ]
Silicate_C_pdf <- density(Silicate_C$Silicate)

# Period D
Silicate_D <- Silicate[Silicate$period %in% c("4"), ]
Silicate_D_pdf <- density(Silicate_D$Silicate)

# Period E
Silicate_E <- Silicate[Silicate$period %in% c("5"), ]
Silicate_E_pdf <- density(Silicate_E$Silicate)

# Period F
Silicate_F <- Silicate[Silicate$period %in% c("6"), ]
Silicate_F_pdf <- density(Silicate_F$Silicate)

# Period G
Silicate_G <- Silicate[Silicate$period %in% c("7"), ]
Silicate_G_pdf <- density(Silicate_G$Silicate)

# Period H
Silicate_H <- Silicate[Silicate$period %in% c("8"), ]
Silicate_H_pdf <- density(Silicate_H$Silicate)


# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(si_all_JSU_meanpoly$y, na.rm = TRUE),
         max(si_all_JSU_meanpoly$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
si_all_JSU_meanpoly_A_pdf  <- dJSU(x, mu = 6.582839, sigma = 5.137672, nu = 6.111668, tau = 10.16478)
si_all_JSU_meanpoly_B_pdf  <- dJSU(x, mu = 6.575298, sigma = 5.227842, nu = 7.591288, tau = 11.21919)
si_all_JSU_meanpoly_C_pdf  <- dJSU(x, mu = 5.724009, sigma = 5.308134, nu = 7.114745, tau = 10.67988)
si_all_JSU_meanpoly_D_pdf  <- dJSU(x, mu = 5.416849, sigma = 5.282102, nu = 7.209079, tau = 10.03653)
si_all_JSU_meanpoly_E_pdf  <- dJSU(x, mu = 6.196346, sigma = 5.668848, nu = 7.828953, tau = 10.78032)
si_all_JSU_meanpoly_F_pdf  <- dJSU(x, mu = 7.020484, sigma = 5.825220, nu = 7.953717, tau = 11.48802)
si_all_JSU_meanpoly_G_pdf  <- dJSU(x, mu = 8.325362, sigma = 5.963357, nu = 8.130878, tau = 12.35281)
si_all_JSU_meanpoly_H_pdf  <- dJSU(x, mu = 10.056940, sigma = 6.194196, nu = 8.408597, tau = 11.78617)



# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

# plot density plots for period A ------------------------------
plot(Silicate_A_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.077))
lines(si_all_JSU_meanpoly_A_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period B ------------------------------
plot(Silicate_B_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_B_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period C ------------------------------
plot(Silicate_C_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_C_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period D ------------------------------
plot(Silicate_D_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_D_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period E ------------------------------
plot(Silicate_E_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_E_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period F ------------------------------
plot(Silicate_F_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_F_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period G ------------------------------
plot(Silicate_G_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_G_pdf, type = "l", lwd = 2, col = "goldenrod") 

# plot density plots for period H ------------------------------
plot(Silicate_H_pdf, lwd = 2, lty = 2, col = "black", type = "l")
lines(si_all_JSU_meanpoly_H_pdf, type = "l", lwd = 2, col = "goldenrod") 





################

# METHOD 1 ==========
logSurv(Silicate$Silicate, prob=0.90, tail="right")


# METHOD 2 ==========
par(mfrow = c(1, 3))
si_m1 <- loglogSurv1(Silicate$Silicate, prob=0.90, title="(a) TYPE I")
si_m2 <- loglogSurv2(Silicate$Silicate, prob=0.90, title="(b) TYPE II")
si_m3 <- loglogSurv3(Silicate$Silicate, prob=0.90, title="(c) TYPE III")
par(mfrow = c(1, 1))


# METHOD 3 ===========

# use some truncated family distribution 
si_m4 <- fitTail(Silicate$Silicate, family=WEI, percentage=10)
si_m5 <- fitTail(Silicate$Silicate, family=LOGNO, percentage=10)
si_m6 <- fitTail(Silicate$Silicate, family=BCPE, percentage=10)
si_m7 <- fitTail(Silicate$Silicate, family=BCT, percentage=10) # got warning message 
si_m8 <- fitTail(Silicate$Silicate, family=JSU, percentage=10)
si_m9 <- fitTail(Silicate$Silicate, family=SHASHo2, percentage=10) # got warning message
si_m10 <- fitTail(Silicate$Silicate, family=GIG, percentage=10)
si_m11 <- fitTail(Silicate$Silicate, family=SEP2, percentage=10)

# check AIC 
AIC(si_m4, si_m5, si_m6, si_m7, si_m8, si_m9, si_m10, si_m11)
## NOTE: BCPE does the best

wp(si_m4, ylim.all = 1)
si_m4_2 <- fitTailAll(Silicate$Silicate, family=WEI) # got like 50 warnings
plot(si_m4_2)



# METHOD 4 ===========
si_f1 <- fitDist(Silicate$Silicate, ncpus=4, parallel="snow" )
si_f1$fits[1:20]
wp(si_f1, ylim.all=1.5)

# now log them
si_z<-log(Silicate$Silicate)
si_f2 <- fitDist(si_z, ncpus=4, parallel="snow")
am_f2$fits[1:10]

wp(am_f2, ylim.all=1)
am_f3 <- histDist(am_z, family="BCPE", nbins=30)

gen.Family("BCPE", "log")
am_f4 <- histDist(Ammonium$Ammonium, family="logBCPE", nbins=30)
AIC(am_f1, am_f4)



# 1. Total number of observations
total_obs <- length(Silicate$Silicate)
print(paste("Total observations:", total_obs))

# 2. Number of observations in the top 10%
num_top_10_percent <- ceiling(total_obs * 0.10)
print(paste("Number of observations in the top 10%:", num_top_10_percent))




############# Back to RQ1
Silicate$Silicate_plus <- Silicate$Silicate + 1

si_f1 <- fitDist(Silicate$Silicate_plus, ncpus=4, parallel="snow")
si_f1$fits[1:20]
wp(si_f1, ylim.all=1)


## NOTE: based on the AIC alone, SEP2 is the best model as it has the lowest AIC

si_static_SEP2 <- gamlss(Silicate_plus ~ 1, family = SEP2(), data = Silicate,
                         #mu.start = mean(Silicate$Silicate_plus), sigma.start = sd(Silicate$Silicate_plus),
                         method = mixed(5,100),
                         control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
AIC(si_static_SEP2)

moment_bucket(si_static_SEP2)
wp(si_static_SEP2, ylim.all=1)


#### now let's look at them tails 
si_m4 <- fitTail(Silicate$Silicate, family=WEI, percentage=10)
si_m5 <- fitTail(Silicate$Silicate, family=LOGNO, percentage=10)
si_m6 <- fitTail(Silicate$Silicate, family=BCPE, percentage=10)
si_m7 <- fitTail(Silicate$Silicate, family=BCT, percentage=10) # got warning message 
si_m8 <- fitTail(Silicate$Silicate, family=JSU, percentage=10)
si_m9 <- fitTail(Silicate$Silicate, family=SHASHo2, percentage=10) # got warning message
si_m10 <- fitTail(Silicate$Silicate, family=GIG, percentage=10) 
si_m11 <- fitTail(Silicate$Silicate, family=SEP2, percentage=10)
si_m12 <- fitTail(Silicate$Silicate, family=BCCG, percentage=10) # got warning 
si_m13 <- fitTail(Silicate$Silicate, family=EXP, percentage=10)
si_m14 <- fitTail(Silicate$Silicate, family=GA, percentage=10)
si_m15 <- fitTail(Silicate$Silicate, family=GG, percentage=10)
si_m16 <- fitTail(Silicate$Silicate, family=IG, percentage=10)


AIC(si_m4, si_m5, si_m6, si_m7, si_m8, si_m9, si_m10, si_m11,
    si_m12, si_m13, si_m14, si_m15, si_m16)


si_static_WEI <- gamlss(Silicate_plus ~ 1, family = WEI(), data = Silicate,
                         #mu.start = mean(Silicate$Silicate_plus), sigma.start = sd(Silicate$Silicate_plus),
                         method = mixed(5,100),
                         control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
AIC(si_static_WEI)

moment_bucket(si_static_WEI)
wp(si_static_WEI, ylim.all=3)




# empirical estimate of moments
si_moments_summary <- Silicate %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Silicate, na.rm = TRUE),
    var = var(Silicate, na.rm = TRUE),
    skew = skewness(Silicate, na.rm = TRUE),
    kurt = kurtosis(Silicate, na.rm = TRUE)
  )
print(si_moments_summary, n=33)



# Silicate (based on periods)
par(mfrow = c(1, 3))
si_m1 <- loglogSurv1(Silicate$Silicate, prob=0.80, title="(a) TYPE I")
si_m2 <- loglogSurv2(Silicate$Silicate, prob=0.80, title="(b) TYPE II")
si_m3 <- loglogSurv3(Silicate$Silicate, prob=0.80, title="(c) TYPE III")
par(mfrow = c(1, 1))

### A: Type III (light tails)
### B: Type II (heavier than normal)*
### C: Type II (heavier than normal)*
### D: Type II (heavier than normal)*
### E: Type II (heavier than normal)*
### F: Type II (heavier than normal)
### G: Type III (light tails)
### H: Type II (heavier than normal)*







si_all_GB2 <- gamlss(Silicate_plus ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ 1,
                         family = GB2(), data = Silicate,
                         mu.start = mean(Silicate$Silicate_plus), sigma.start = sd(Silicate$Silicate_plus), nu.start = 2, tau.start = 2,
                         method = mixed(5,300),
                         control = gamlss.control(n.cyc = 300, c.crit = 0.01, trace = TRUE))


library(gamlss.tr)
gen.trun(0,"SEP4",type="left")

si_all_SEP4tr <- gamlss(Silicate_plus ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                           family = SEP4tr(), data = Silicate,
                           mu.start = mean(Silicate$Silicate_plus), sigma.start = sd(Silicate$Silicate_plus),
                           method = mixed(5,500),
                           control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))



# get the parameters for the year periods ======
si_all_SEPtr_meanpoly_param <- best_model_param(si_all_SEP4tr, Silicate)

# 4. Get PDF of those three year combination, then compare it to non-parametric density estimate ------------------------------

# get pdf for am_all_BCPE_sigmapoly -----
x <- seq(min(si_all_GB2$y, na.rm = TRUE),
         max(si_all_GB2$y, na.rm = TRUE))

# natural scale (back-transform) for all periods 
si_all_GB2_meanpoly_A_pdf  <- dGB2(x, mu = 6.582839, sigma = 5.137672, nu = 6.111668, tau = 10.16478)
si_all_GB2_meanpoly_B_pdf  <- dGB2(x, mu = 6.575298, sigma = 5.227842, nu = 7.591288, tau = 11.21919)
si_all_GB2_meanpoly_C_pdf  <- dGB2(x, mu = 5.724009, sigma = 5.308134, nu = 7.114745, tau = 10.67988)
si_all_GB2_meanpoly_D_pdf  <- dGB2(x, mu = 5.416849, sigma = 5.282102, nu = 7.209079, tau = 10.03653)
si_all_GB2_meanpoly_E_pdf  <- dGB2(x, mu = 6.196346, sigma = 5.668848, nu = 7.828953, tau = 10.78032)
si_all_GB2_meanpoly_F_pdf  <- dGB2(x, mu = 7.020484, sigma = 5.825220, nu = 7.953717, tau = 11.48802)
si_all_GB2_meanpoly_G_pdf  <- dGB2(x, mu = 8.325362, sigma = 5.963357, nu = 8.130878, tau = 12.35281)
si_all_GB2_meanpoly_H_pdf  <- dGB2(x, mu = 10.056940, sigma = 6.194196, nu = 8.408597, tau = 11.78617)



# get the parameters for the year periods ======
si_all_GB2_meanpoly_param <- best_model_param(si_all_GB2, Silicate)





si_nokurt_SEP4tr <- gamlss(Silicate_plus ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ 1,
                        family = SEP4tr(), data = Silicate,
                        mu.start = mean(Silicate$Silicate_plus), sigma.start = sd(Silicate$Silicate_plus),
                        method = mixed(5,500),
                        control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))

si_nokurt_SEP4tr_meanpoly_param <- best_model_param(si_nokurt_SEP4tr, Silicate)















library(gamlss.tr)
gen.trun(0,"SEP4",type="left")
gen.trun(0,"SHASHo",type="left")

si_all_BCPEo_spline <- gamlss(Silicate + 1 ~ cs(year) + month, sigma.fo = ~ cs(year) + month, nu.fo = ~ year + month, tau.fo = ~ cs(year) + month,
                        family = SEP4tr(), data = Silicate,
                        #mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                        method = mixed(5,500),
                        control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))



si_intercept_SEP4tr<- gamlss(Silicate ~ 1,
                               family = SEP4tr(), data = Silicate,
                               #mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                               method = mixed(5,500),
                               control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))
summary(si_intercept_SEP4tr)


si_intercept_SEP2<- gamlss(Silicate ~ 1,
                             family = SEP2(), data = Silicate,
                             #mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                             method = mixed(5,500),
                             control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))
summary(si_intercept_SEP2)




gen.trun(0,"SEP1",type="left")
gen.trun(0,"SEP2",type="left")
gen.trun(0,"SEP3",type="left")
gen.trun(0,"SEP4",type="left")
gen.trun(0,"SHASHo",type="left") # not working 
gen.trun(0,"SHASHo2",type="left") # not working 
gen.trun(0,"SHASH",type="left")

##### use bins in the model 
Silicate$period <- as.numeric(Silicate$period)

si_all_BCPEo_bins_spline <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ pb(period) + month, nu.fo = ~ pb(period) + month, tau.fo = ~ pb(period) + month,
                                   family = BCPEo(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))

si_all_BCPEo_bins_spline_param <- best_model_param(si_all_BCPEo_bins_spline, Silicate)

set.seed(123)
si_all_BCPEo_bins_spline_predict <- list(
  "1" = rBCPEo(122, mu = 4.371509, sigma = 37.88017, nu = 3.227595, tau = 0.3016460),
  "2" = rBCPEo(466, mu = 3.398992, sigma = 29.41983, nu = 2.654444, tau = 0.3334140),
  "3" = rBCPEo(547, mu = 2.048156, sigma = 25.51359, nu = 3.066303, tau = 0.3494071),
  "4" = rBCPEo(755, mu = 1.543980, sigma = 22.85957, nu = 2.606879, tau = 0.3775562),
  "5" = rBCPEo(961, mu = 1.717384, sigma = 19.48771, nu = 3.147447, tau = 0.4021955),
  "6" = rBCPEo(980, mu = 1.200094, sigma = 17.52896, nu = 3.016291, tau = 0.4252636),
  "7" = rBCPEo(973, mu = 6.213807, sigma = 15.18648, nu = 2.847544, tau = 0.4522932),
  "8" = rBCPEo(968, mu = 5.744028, sigma = 13.39220, nu = 4.352010, tau = 0.4821415)
)


sii_all_BCPEo_bins_spline_predict_long <- stack(si_all_BCPEo_bins_spline_predict)
names(sii_all_BCPEo_bins_spline_predict_long) <- c("pred", "period")
sii_all_BCPEo_bins_spline_predict_long$period <- as.numeric(sii_all_BCPEo_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Silicate, aes(x = Silicate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = sii_all_BCPEo_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Silicate (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1966", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  









# bins time-varying parameters 
si_all_BCPEo_bins_spline <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ period + month, nu.fo = ~ period + month, tau.fo = ~ period + month,
                                   family = BCPEo(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

si_nokurt_BCPEo_bins_spline <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ period + month, nu.fo = ~ period + month, tau.fo = ~ 1,
                                   family = BCPEo(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

si_noskew_BCPEo_bins_spline <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ period + month, nu.fo = ~ 1, tau.fo = ~ period + month,
                                   family = BCPEo(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))

si_mean_sigma_BCPEo_bins_spline <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ period + month, nu.fo = ~ 1, tau.fo = ~ 1,
                                   family = BCPEo(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 500, c.crit = 0.01, trace = TRUE))

si_meanonly_BCPEo_bins_spline <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                                   family = BCPEo(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

si_intercpet_BCPEo <- gamlss(Silicate + 1  ~ 1,
                                        family = BCPEo(), data = Silicate,
                                        mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                        #method = mixed(5,100),
                                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))






AIC(si_all_BCPEo_bins_spline, si_nokurt_BCPEo_bins_spline,
    si_noskew_BCPEo_bins_spline,si_mean_sigma_BCPEo_bins_spline,
    si_meanonly_BCPEo_bins_spline)


### RESULTS 
# si_nokurt_BCPEo_bins_spline      28978.09
# si_mean_sigma_BCPEo_bins_spline  29099.71 # convergence 
# si_all_BCPEo_bins_spline         29202.39
# si_noskew_BCPEo_bins_spline      29241.84
# si_meanonly_BCPEo_bins_spline    29728.66



## since nokurt has lowest AIC, let's redo the pdfs 
si_nokurt_BCPEo_bins_spline_param <- best_model_param(si_nokurt_BCPEo_bins_spline, Silicate)

set.seed(123)
si_nokurt_BCPEo_bins_spline_predict <- list(
  "1" = rBCPEo(122, mu = 2.083577, sigma = 33.389928, nu = 2.157591, tau = 0.4302056),
  "2" = rBCPEo(466, mu = 1.895755, sigma = 24.912906, nu = 2.306390, tau = 0.4302056),
  "3" = rBCPEo(547, mu = 1.587141, sigma = 18.062264, nu = 2.434490, tau = 0.4302056),
  "4" = rBCPEo(755, mu = 1.651377, sigma = 13.349290, nu = 2.572304, tau = 0.4302056),
  "5" = rBCPEo(961, mu = 1.852747, sigma = 10.937368, nu = 2.705547, tau = 0.4302056),
  "6" = rBCPEo(980, mu = 1.587014, sigma = 8.295265, nu = 2.834858, tau = 0.4302056),
  "7" = rBCPEo(973, mu = 7.979859, sigma = 6.145398, nu = 2.970139, tau = 0.4302056),
  "8" = rBCPEo(968, mu = 6.399850, sigma = 4.798561, nu = 3.097104, tau = 0.4302056)
)


sii_nokurt_BCPEo_bins_spline_predict_long <- stack(si_nokurt_BCPEo_bins_spline_predict)
names(sii_nokurt_BCPEo_bins_spline_predict_long) <- c("pred", "period")
sii_nokurt_BCPEo_bins_spline_predict_long$period <- as.numeric(sii_nokurt_BCPEo_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Silicate, aes(x = Silicate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = sii_nokurt_BCPEo_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Silicate (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1966", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  




si_intercpet_SEP2 <- gamlss(Silicate  ~ 1,
                             family = SEP2(), data = Silicate,
                             mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                             #method = mixed(5,100),
                             control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

si_interpt_models <- fitDist(Silicate$Silicate)
si_interpt_models$fits[1:20]


moment_bucket(si_intercpet_SEP2) +
resid_wp(si_intercpet_SEP2) + theme(plot.title = element_blank())





si_nokurt_BCPEo_bins_spline_2 <- gamlss(Silicate + 1  ~ pb(period) + month, sigma.fo = ~ pb(period) + month, nu.fo = ~ pb(period) + month, tau.fo = ~ 1,
                                      family = BCPEo(), data = Silicate,
                                      mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                      #method = mixed(5,100),
                                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


si_nokurt_BCPEo_bins_spline_3 <- gamlss(Silicate + 1  ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ 1,
                                        family = BCPEo(), data = Silicate,
                                        mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                        #method = mixed(5,100),
                                        control = gamlss.control(n.cyc = 300, c.crit = 0.01, trace = TRUE))








# monte carlo simulation 


si_all_BCPEo_bins_spline_predict <- list(
  "1" = rBCPEo(122, mu = 4.371509, sigma = 37.88017, nu = 3.227595, tau = 0.3016460),
  "2" = rBCPEo(466, mu = 3.398992, sigma = 29.41983, nu = 2.654444, tau = 0.3334140),
  "3" = rBCPEo(547, mu = 2.048156, sigma = 25.51359, nu = 3.066303, tau = 0.3494071),
  "4" = rBCPEo(755, mu = 1.543980, sigma = 22.85957, nu = 2.606879, tau = 0.3775562),
  "5" = rBCPEo(961, mu = 1.717384, sigma = 19.48771, nu = 3.147447, tau = 0.4021955),
  "6" = rBCPEo(980, mu = 1.200094, sigma = 17.52896, nu = 3.016291, tau = 0.4252636),
  "7" = rBCPEo(973, mu = 6.213807, sigma = 15.18648, nu = 2.847544, tau = 0.4522932),
  "8" = rBCPEo(968, mu = 5.744028, sigma = 13.39220, nu = 4.352010, tau = 0.4821415)
)

library(gamlss.dist)
library(moments)

############################## bin 1
mu_si_bin1   <-     4.371509
sigma_si_bin1 <-    37.88017
nu_si_bin1  <-    3.227595
tau_si_bin1  <- 0.3016460

set.seed(1)
si_bin1_sim <- rBCPEo(N, mu = mu_si_bin1, sigma = sigma_si_bin1, nu = nu_si_bin1, tau = tau_si_bin1)

# get skewness and kurtosis 
kurtosis(si_bin1_sim) - 3 # 2.29653
skewness(si_bin1_sim) # 1.284017
var(si_bin1_sim) # 35.82415
mean(si_bin1_sim) # 10.55564


############################# bin 2
mu_si_bin2   <-   3.398992  
sigma_si_bin2 <-    29.41983
nu_si_bin2  <-    2.654444
tau_si_bin2  <- 0.3334140
  
set.seed(2)
si_bin2_sim <- rBCPEo(N, mu = mu_si_bin2, sigma = sigma_si_bin2, nu = nu_si_bin2, tau = tau_si_bin2)

# get skewness and kurtosis 
kurtosis(si_bin2_sim) - 3 # 3.443164
skewness(si_bin2_sim) # 1.52059
var(si_bin2_sim) # 37.03401
mean(si_bin2_sim) # 9.270082

 
############################# bin 3
mu_si_bin3   <-    2.048156 
sigma_si_bin3 <-    25.51359
nu_si_bin3  <-    3.066303
tau_si_bin3  <- 0.3494071
  
set.seed(3)
si_bin3_sim <- rBCPEo(N, mu = mu_si_bin3, sigma = sigma_si_bin3, nu = nu_si_bin3, tau = tau_si_bin3)

# get skewness and kurtosis 
kurtosis(si_bin3_sim) - 3 # 1.898579
skewness(si_bin3_sim) # 1.188436
var(si_bin3_sim) # 7.277168
mean(si_bin3_sim) # 4.887771


############################# bin 4
mu_si_bin4   <-     1.543980
sigma_si_bin4 <-    22.85957
nu_si_bin4  <-    2.606879
tau_si_bin4  <- 0.3775562
  
set.seed(4)
si_bin4_sim <- rBCPEo(N, mu = mu_si_bin4, sigma = sigma_si_bin4, nu = nu_si_bin4, tau = tau_si_bin4)

# get skewness and kurtosis 
kurtosis(si_bin4_sim) - 3 # 2.733718
skewness(si_bin4_sim) # 1.376977
var(si_bin4_sim) # 6.757346
mean(si_bin4_sim) # 4.153336


############################# bin 5
mu_si_bin5   <-    1.717384 
sigma_si_bin5 <-    19.48771
nu_si_bin5  <-    3.147447
tau_si_bin5  <- 0.4021955
  
set.seed(5)
si_bin5_sim <- rBCPEo(N, mu = mu_si_bin5, sigma = sigma_si_bin5, nu = nu_si_bin5, tau = tau_si_bin5)

# get skewness and kurtosis 
kurtosis(si_bin5_sim) - 3 # 1.242842
skewness(si_bin5_sim) # 0.9962532
var(si_bin5_sim) # 3.954322
mean(si_bin5_sim) # 3.947018


############################# bin 6
mu_si_bin6   <-    1.200094 
sigma_si_bin6 <-    17.52896
nu_si_bin6  <-    3.016291
tau_si_bin6  <- 0.4252636
  
set.seed(6)
si_bin6_sim <- rBCPEo(N, mu = mu_si_bin6, sigma = sigma_si_bin6, nu = nu_si_bin6, tau = tau_si_bin6)

# get skewness and kurtosis 
kurtosis(si_bin6_sim) - 3 # 1.309219
skewness(si_bin6_sim) # 1.006845
var(si_bin6_sim) # 2.049662
mean(si_bin6_sim) # 2.807751


############################# bin 7
mu_si_bin7  <-    6.213807 
sigma_si_bin7 <-    15.18648
nu_si_bin7  <-    2.847544
tau_si_bin7  <- 0.4522932
  
set.seed(7)
si_bin7_sim <- rBCPEo(N, mu = mu_si_bin7, sigma = sigma_si_bin7, nu = nu_si_bin7, tau = tau_si_bin7)

# get skewness and kurtosis 
kurtosis(si_bin7_sim) - 3 # 1.301179
skewness(si_bin7_sim) # 1.015378
var(si_bin7_sim) # 59.05768
mean(si_bin7_sim) # 14.69842


############################# bin 8
mu_si_bin8   <-    5.744028 
sigma_si_bin8 <-    13.39220
nu_si_bin8  <-    4.352010
tau_si_bin8  <- 0.4821415
  
set.seed(8)
si_bin8_sim <- rBCPEo(N, mu = mu_si_bin8, sigma = sigma_si_bin8, nu = nu_si_bin8, tau = tau_si_bin8)

# get skewness and kurtosis 
kurtosis(si_bin8_sim) - 3 # -0.005948124
skewness(si_bin8_sim) # 0.4295109
var(si_bin8_sim) # 13.16864
mean(si_bin8_sim) # 10.68833





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
si_bin1_sim_se <- bootstrap_moments_sim_only(si_bin1_sim) 
si_bin2_sim_se <- bootstrap_moments_sim_only(si_bin2_sim) 
si_bin3_sim_se <- bootstrap_moments_sim_only(si_bin3_sim) 
si_bin4_sim_se <- bootstrap_moments_sim_only(si_bin4_sim) 
si_bin5_sim_se <- bootstrap_moments_sim_only(si_bin5_sim) 
si_bin6_sim_se <- bootstrap_moments_sim_only(si_bin6_sim) 
si_bin7_sim_se <- bootstrap_moments_sim_only(si_bin7_sim) 
si_bin8_sim_se <- bootstrap_moments_sim_only(si_bin8_sim) 

si_bin8_sim_se$se

si_bins_df <- tibble::tibble(
  year_mid = c(1966, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  mean = c(10.556,9.270,4.888,4.153,3.947,2.808,14.698,10.688),
  variance = c(35.824,37.034,7.277,6.757,3.954,2.050,59.058,13.169),
  skewness = c(1.284,1.521,1.188,1.377,0.996,1.007,1.015,0.430),
  kurtosis = c(2.297,3.443,1.899,2.733,1.243,1.309,1.301,-0.006)
)

si_se_df <- tibble::tibble(
  year_mid = si_bins_df$year_mid,
  se_mean = c(0.008327026,0.008557848,0.003864204,0.003659303,
              0.002793657,0.002072796,0.010651489,0.004965832),
  se_variance = c(0.104648253,0.120793546,0.020797956,0.020533033,
                  0.009885407,0.005250138,0.150520472,0.026137703),
  se_skew = c(0.007000402,0.008962337,0.005896449,0.007918784,
              0.005110358,0.005213796,0.005027148,0.003357494),
  se_kurt = c(0.051966837,0.078445843,0.036832837,0.068050130,
              0.027521371,0.028117714,0.025937383,0.009033781)
)

# join and pivot to long form
si_df_long <- si_bins_df %>%
  left_join(si_se_df, by = "year_mid") %>%
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


ggplot(si_df_long, aes(x = year_mid, y = value)) +
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







####====== combine moments of all nutrients ======#####

#### kurtosis #####################################################################

kurt_all <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_k = c(2.297,3.443,1.899,2.733,1.243,1.309,1.301,-0.006),
  Nitrite_k = c(1.371,1.342,1.075,0.997,2.100,2.161,0.939,0.105),
  Phosphate_k = c(3.716,2.716,0.664,0.407,0.750,1.202,2.115,3.937),
  Nitrate_k = c(1.468,2.659,2.186,1.601,0.685,-0.012,0.095,1.687),
  DIN_k = c(-0.283,0.376,1.038,2.001,2.122,1.044,0.395,0.419),
  Ammonium_k = c(0.370,0.398,0.273,0.586,1.636,3.264,5.046,4.623)
)

kurt_all_se <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_k_se = c(0.051966837,0.078445843,0.036832837,0.068050130,
               0.027521371,0.028117714,0.025937383,0.009033781),
  Nitrite_k_se = c(0.0285856634,0.0301624652,0.0236317279,0.0273073230,
                0.0500960830,0.0477571198,0.0206259617,0.0103299206),
  Phosphate_k_se = c(0.0910671947,0.0556927519,0.0187235385,0.0169571335,
               0.0214152375,0.0287815614,0.048765570,0.0922350562),
  Nitrate_k_se = c(0.027114858,0.052560952,0.037299409,0.030475993,
               0.015516501,0.008719303,0.009396691,0.023431861),
  DIN_k_se = c(0.006601600,0.01536872,0.025791753,0.069681021,
                0.056682017,0.02646837,0.014115509,0.01174154),
  Ammonium_k_se = c(0.013622697,0.015933674,0.012790616,0.018134935,
               0.034477138,0.067954178,0.108433020,0.093706115)
)

kurt_all_df <- merge(kurt_all, kurt_all_se, by = "year_mid")

### plot
# ggplot(kurt_all_df, aes(x = year_mid)) +
#   geom_line(aes(y = si_k), color = "azure4") +
#   geom_point(aes(y = si_k), color = "azure4") +
#   geom_errorbar(aes(ymin = si_k - 1.96 * si_k_se, ymax = si_k + 1.96 * si_k_se), 
#                 width = 0.7, color = "azure4") +
#   geom_line(aes(y = nii_k), color = "steelblue") +
#   geom_point(aes(y = nii_k), color = "steelblue") +
#   geom_errorbar(aes(ymin = nii_k - 1.96 * nii_k_se, ymax = nii_k + 1.96 * nii_k_se), 
#                 width = 0.7, color = "steelblue") +
#   geom_line(aes(y = ph_k), color = "darkolivegreen") +
#   geom_point(aes(y = ph_k), color = "darkolivegreen") +
#   geom_errorbar(aes(ymin = ph_k - 1.96 * ph_k_se, ymax = ph_k + 1.96 * ph_k_se), 
#                 width = 0.7, color = "darkolivegreen") +
#   geom_line(aes(y = am_k), color = "black") +
#   geom_point(aes(y = am_k), color = "black") +
#   geom_errorbar(aes(ymin = am_k - 1.96 * am_k_se, ymax = am_k + 1.96 * am_k_se), 
#                 width = 0.7, color = "black") +
#   geom_line(aes(y = DIN_k), color = "chocolate") +
#   geom_point(aes(y = DIN_k), color = "chocolate") +
#   geom_errorbar(aes(ymin = DIN_k - 1.96 * DIN_k_se, ymax = DIN_k + 1.96 * DIN_k_se), 
#                 width = 0.7, color = "chocolate") +
#   geom_line(aes(y = ni_k), color = "goldenrod") +
#   geom_point(aes(y = ni_k), color = "goldenrod") +
#   geom_errorbar(aes(ymin = ni_k - 1.96 * ni_k_se, ymax = ni_k + 1.96 * ni_k_se), 
#                 width = 0.7, color = "goldenrod") +
#   geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
#   xlab("Year bins") +
#   ylab("Excess kurtosis") +
#   theme_minimal()

kurt_long <- kurt_all_df %>%
  pivot_longer(
    cols = -year_mid,
    names_to = c("nutrient", "metric"),
    names_pattern = "^(.*)_(k(?:_se)?)$",
    values_drop_na = FALSE
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(kurtosis = k, se = k_se) %>%
  mutate(nutrient = str_remove(nutrient, "_$")) %>%
  arrange(nutrient, year_mid)

# plot
ggplot(kurt_long, aes(x = year_mid, y = kurtosis, color = nutrient, group = nutrient)) +
  geom_errorbar(aes(ymin = kurtosis - 1.96 * se, ymax = kurtosis + 1.96 * se),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Silicate" = "azure4", "Nitrite" = "steelblue", "Phosphate" = "darkolivegreen",
                                "Ammonium" = "black", "DIN" = "chocolate", "Nitrate" = "goldenrod")) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
  labs(x = "Year bins", y = "Excess kurtosis", color = "Nutrients") +
  scale_x_continuous(breaks = unique(kurt_long$year_mid)) +
  theme_bw() +
  theme(legend.position = "right")



### skewness #####################################################################
skew_all <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_sk = c(1.284,1.521,1.188,1.377,0.996,1.007,1.015,0.430),
  Nitrite_sk = c(0.975,1.036,0.995,0.909,1.084,1.130,0.970,0.732),
  Phosphate_sk = c(1.351,1.052,0.722,0.657,0.723,0.770,0.933,1.314),
  Nitrate_sk = c(1.059,1.330,1.213,1.084,0.870,0.646,0.737,1.281),
  DIN_sk = c(0.436,0.650,0.828,1.039,1.030,0.781,0.735,0.967),
  Ammonium_sk = c(0.610,0.633,0.563,0.652,0.968,1.338,1.640,1.340)
)

skew_all_se <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_sk_se = c(0.007000402,0.008962337,0.005896449,0.007918784,
              0.005110358,0.005213796,0.005027148,0.003357494),
  Nitrite_sk_se = c(0.0053590092,0.0052255988,0.0046475068,0.0049345745,
               0.0073008419,0.0070670580,0.0043976990,0.0031003094),
  Phosphate_sk_se = c(0.0108495722,0.0076980327,0.0043515345,0.0040797363,
               0.0046396849,0.0056142459,0.007746112,0.0111093414),
  Nitrate_sk_se = c(0.005153090,0.007194377,0.006110511,0.005324222,
               0.003817211,0.002978047,0.003031812,0.004510697),
  DIN_sk_se = c(0.002754898,0.00394651,0.005058384,0.008070872,
                0.008002874,0.00535313,0.003689457,0.00320749),
  Ammonium_sk_se = c(0.003876006,0.004041878,0.003716181,0.004287422,
               0.006093326,0.009030545,0.011670575,0.011669906)
)

skew_all_df <- merge(skew_all, skew_all_se, by = "year_mid")

### plot
# ggplot(skew_all_df, aes(x = year_mid)) +
#   geom_line(aes(y = si_sk), color = "azure4") +
#   geom_point(aes(y = si_sk), color = "azure4") +
#   geom_errorbar(aes(ymin = si_sk - 1.96 * si_sk_se, ymax = si_sk + 1.96 * si_sk_se), 
#                 width = 0.7, color = "azure4") +
#   geom_line(aes(y = nii_sk), color = "steelblue") +
#   geom_point(aes(y = nii_sk), color = "steelblue") +
#   geom_errorbar(aes(ymin = nii_sk - 1.96 * nii_sk_se, ymax = nii_sk + 1.96 * nii_sk_se), 
#                 width = 0.7, color = "steelblue") +
#   geom_line(aes(y = ph_sk), color = "darkolivegreen") +
#   geom_point(aes(y = ph_sk), color = "darkolivegreen") +
#   geom_errorbar(aes(ymin = ph_sk - 1.96 * ph_sk_se, ymax = ph_sk + 1.96 * ph_sk_se), 
#                 width = 0.7, color = "darkolivegreen") +
#   geom_line(aes(y = am_sk), color = "black") +
#   geom_point(aes(y = am_sk), color = "black") +
#   geom_errorbar(aes(ymin = am_sk - 1.96 * am_sk_se, ymax = am_sk + 1.96 * am_sk_se), 
#                 width = 0.7, color = "black") +
#   geom_line(aes(y = DIN_sk), color = "chocolate") +
#   geom_point(aes(y = DIN_sk), color = "chocolate") +
#   geom_errorbar(aes(ymin = DIN_sk - 1.96 * DIN_sk_se, ymax = DIN_sk + 1.96 * DIN_sk_se), 
#                 width = 0.7, color = "chocolate") +
#   geom_line(aes(y = ni_sk), color = "goldenrod") +
#   geom_point(aes(y = ni_sk), color = "goldenrod") +
#   geom_errorbar(aes(ymin = ni_sk - 1.96 * ni_sk_se, ymax = ni_sk + 1.96 * ni_sk_se), 
#                 width = 0.7, color = "goldenrod") +
#   geom_hline(yintercept = 0, color = "darkred", linetype = "dashed") +
#   xlab("Year bins") +
#   ylab("Skewness") +
#   theme_minimal()

skew_long <- skew_all_df %>%
  pivot_longer(
    cols = -year_mid,
    names_to = c("nutrient", "metric"),
    names_pattern = "^(.*)_(sk(?:_se)?)$",
    values_drop_na = FALSE
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(Skewness = sk, se = sk_se) %>%
  mutate(nutrient = str_remove(nutrient, "_$")) %>%
  arrange(nutrient, year_mid)

# plot
ggplot(skew_long, aes(x = year_mid, y = Skewness, color = nutrient, group = nutrient)) +
  geom_errorbar(aes(ymin = Skewness - 1.96 * se, ymax = Skewness + 1.96 * se),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Silicate" = "azure4", "Nitrite" = "steelblue", "Phosphate" = "darkolivegreen",
                                "Ammonium" = "black", "DIN" = "chocolate", "Nitrate" = "goldenrod")) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(x = "Year bins", y = "Skewness", color = "Nutrients") +
  scale_x_continuous(breaks = unique(skew_long$year_mid)) +
  theme_classic() +
  theme(legend.position = "none") +
  ggplot(kurt_long, aes(x = year_mid, y = kurtosis, color = nutrient, group = nutrient)) +
  geom_errorbar(aes(ymin = kurtosis - 1.96 * se, ymax = kurtosis + 1.96 * se),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Silicate" = "azure4", "Nitrite" = "steelblue", "Phosphate" = "darkolivegreen",
                                "Ammonium" = "black", "DIN" = "chocolate", "Nitrate" = "goldenrod")) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  labs(x = "Year bins", y = "Kurtosis", color = "Nutrients") +
  scale_x_continuous(breaks = unique(kurt_long$year_mid)) +
  theme_classic() +
  theme(legend.position = "none")




### Variance #####################################################################
var_all <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_var = c(35.824,37.034,7.277,6.757,3.954,2.050,59.058,13.169),
  Nitrite_var = c(0.201,0.235,0.222,0.196,0.222,0.246,0.246,0.148),
  Phosphate_var = c(0.109,0.105,0.089,0.090,0.102,0.095,0.080,0.097),
  Nitrate_var = c(26.601,27.231,25.931,34.721,51.815,72.015,139.563,364.277),
  DIN_var = c(53.377,60.253,60.318,58.975,61.107,80.593,158.502,474.165),
  Ammonium_var = c(9.160,10.878,10.671,7.757,5.270,2.868,2.340,1.952)
)

var_all_se <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_var_se = c(0.104648253,0.120793546,0.020797956,0.020533033,
               0.009885407,0.005250138,0.150520472,0.026137703),
  Nitrite_var_se = c(0.0005252876,0.0006202958,0.0005477940,0.0004864865,
                 0.0006394150,0.0007024575,0.0006084682,0.0003040004),
  Phosphate_var_se = c(0.0003669594,0.0002975767,0.0002083750,0.0002004666,
                0.0002367871,0.0002410548,0.000227958,0.0003337661),
  Nitrate_var_se = c(0.069718673,0.085008131,0.073731347,0.094024152,
                0.120707255,0.138510137,0.286956477,0.995877110),
  DIN_var_se = c(0.098832581,0.13149763,0.149360467,0.167194544,
                 0.177430045,0.19541194,0.349876593,1.03674250),
  Ammonium_var_se = c(0.020178283,0.024048321,0.022752320,0.017471896,
                0.014151229,0.009448144,0.008731005,0.007057966)
)

var_all_df <- merge(var_all, var_all_se, by = "year_mid")

var_long <- var_all_df %>%
  pivot_longer(
    cols = -year_mid,
    # capture nutrient name and the metric (either "var" or "var_se")
    names_to = c("nutrient", "metric"),
    names_pattern = "^(.*)_(var(?:_se)?)$",
    values_drop_na = FALSE
  ) %>%
  # now metric will be "var" or "var_se"; widen so we have columns var and var_se
  pivot_wider(names_from = metric, values_from = value) %>%
  # rename to nice names
  rename(variance = var, se = var_se) %>%
  # optional: clean nutrient names if they contain trailing underscores
  mutate(nutrient = str_remove(nutrient, "_$")) %>%
  arrange(nutrient, year_mid)

var_long <- var_long %>%
  group_by(nutrient) %>%
  mutate(
    var_mean = mean(variance, na.rm = TRUE),
    var_sd   = sd(variance, na.rm = TRUE),
    # protect against zero sd (if var_sd == 0, z and se_z will be NA)
    z = if_else(var_sd > 0, (variance - var_mean) / var_sd, NA_real_),
    se_z = if_else(var_sd > 0, se / var_sd, NA_real_),
    z_lo = z - 1.96 * se_z,
    z_hi = z + 1.96 * se_z
  ) %>%
  ungroup()

# plot
ggplot(var_long, aes(x = year_mid, y = z, color = nutrient, group = nutrient)) +
  geom_errorbar(aes(ymin = z_lo, ymax = z_hi),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Silicate" = "azure4", "Nitrite" = "steelblue", "Phosphate" = "darkolivegreen",
                                "Ammonium" = "black", "DIN" = "chocolate", "Nitrate" = "goldenrod")) +
  labs(x = "Year bins", y = "Z (Variance)", color = "Nutrients") +
  scale_x_continuous(breaks = unique(var_long$year_mid)) +
  theme_bw() +
  theme(legend.position = "right")


### mean #####################################################################
mean_all <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_m = c(10.556,9.270,4.888,4.153,3.947,2.808,14.698,10.688),
  Nitrite_m = c(0.687,0.775,0.831,0.842,0.886,0.857,0.788,0.657),
  Phosphate_m = c(0.647,0.723,0.816,0.883,0.894,0.815,0.694,0.640),
  Nitrate_m = c(7.597,7.333,7.640,9.587,13.576,18.168,23.334,26.575),
  DIN_m = c(17.871,17.951,18.141,19.421,21.283,24.211,28.322,33.889),
  Ammonium_m = c(9.409,8.744,8.032,7.117,5.811,4.674,3.771,3.225)
)

mean_all_se <- tibble::tibble(
  year_mid = c(1966, 1970, 1974, 1978, 1982, 1986, 1990, 1994),
  Silicate_m_se = c(0.008327026,0.008557848,0.003864204,0.003659303,
              0.002793657,0.002072796,0.010651489,0.004965832),
  Nitrite_m_se = c(0.0006448352,0.0006854193,0.0006687207,0.0006316300,
               0.0006792097,0.0007030647,0.0007034426,0.0005531132),
  Phosphate_m_se = c(0.0004739638,0.0004553520,0.0004190247,0.0004307061,
              0.0004425233,0.0004387041,0.000402048,0.0004405365),
  Nitrate_m_se = c(0.007424578,0.007710453,0.006965149,0.008256720,
              0.010187013,0.011729352,0.016204944,0.026905917),
  DIN_m_se = c(0.010306205,0.01093010,0.010975780,0.010889395,
               0.011115839,0.01278470,0.018202419,0.03081060),
  Ammonium_m_se = c(0.004436985,0.004564469,0.004453601,0.003951382,
              0.003271781,0.002427279,0.002168838,0.001981363)
)

mean_all_df <- merge(mean_all, mean_all_se, by = "year_mid")

mean_long <- mean_all_df %>%
  pivot_longer(
    cols = -year_mid,
    names_to = c("nutrient", "metric"),
    names_pattern = "^(.*)_(m(?:_se)?)$",
    values_drop_na = FALSE
  ) %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  rename(mean = m, se = m_se) %>%
  mutate(nutrient = str_remove(nutrient, "_$")) %>%
  arrange(nutrient, year_mid)

mean_long <- mean_long %>%
  group_by(nutrient) %>%
  mutate(
    mean_mean = mean(mean, na.rm = TRUE),
    mean_sd   = sd(mean, na.rm = TRUE),
    z = if_else(mean_sd > 0, (mean - mean_mean) / mean_sd, NA_real_),
    se_z = if_else(mean_sd > 0, se / mean_sd, NA_real_),
    z_lo = z - 1.96 * se_z,
    z_hi = z + 1.96 * se_z
  ) %>%
  ungroup()

# plot
ggplot(mean_long, aes(x = year_mid, y = z, color = nutrient, group = nutrient)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = z_lo, ymax = z_hi),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  scale_color_manual(values = c("Silicate" = "azure4", "Nitrite" = "steelblue", "Phosphate" = "darkolivegreen",
                               "Ammonium" = "black", "DIN" = "chocolate", "Nitrate" = "goldenrod")) +
  labs(y = "Z (Mean)", color = "Nutrients") +
  scale_x_continuous(breaks = unique(mean_long$year_mid)) +
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.title.x = element_blank()) +
  ggplot(var_long, aes(x = year_mid, y = z, color = nutrient, group = nutrient)) +
  geom_errorbar(aes(ymin = z_lo, ymax = z_hi),
                position = position_dodge(width = 0),
                width = 4, size = 0.6, alpha = 1) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_color_manual(values = c("Silicate" = "azure4", "Nitrite" = "steelblue", "Phosphate" = "darkolivegreen",
                                "Ammonium" = "black", "DIN" = "chocolate", "Nitrate" = "goldenrod")) +
  labs(y = "Z (Variance)", color = "Nutrients") +
  scale_x_continuous(breaks = unique(var_long$year_mid)) +
  theme_classic() +
  theme(legend.position = "right") +
  theme(axis.title.x = element_blank())




