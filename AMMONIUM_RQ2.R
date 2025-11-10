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
library(ggpubr)

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


am_time_AIC <- data.frame(
  model = c("mean only(time)", "mean & sigma(time)", 
            "constant skewness", "constant kurtosis", "all parameter(time)"),
  AIC = c("29746.79", "29397.90", "29369.81", "29156.90",
               "29124.04")
)
# change column names 
colnames(am_time_AIC) <- c("time-varying parameters", "AIC")

# make them table 
am_time_AIC_table <- ggtexttable(am_time_AIC, rows = NULL, theme = ttheme("light", base_size = 12))



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
  theme_classic()

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
am_all_BCPE_sigmapoly_param <- best_model_param(am_all_BCPE_sigmapoly, Ammonium)
# for when higher-moments are not changing through time 
am_meansigma_BCPE_param <- best_model_param(am_mean_and_sigma_BCPE, Ammonium)



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


am_BCPE_param <- data.frame(
  years = c("1962 - 1966", "1967 - 1970", "1971 - 1974", "1975 - 1978", 
            "1979 - 1982", "1983 - 1986", "1987 - 1990", "1991 - 1994"),
  mu = c("9.05", "8.02", "7.15", "6.26",
          "5.40", "4.53", "3.67", "2.76"),
  sigma = c("0.403", "0.448", "0.479", "0.512",
            "0.532", "0.545", "0.546", "0.546"),
  nu = c("0.925", "0.802", "0.729", "0.627",
         "0.565", "0.477", "0.392", "0.303"),
  tau = c("1.950330", "1.884630", "1.805774", "1.739688",
          "1.669930", "1.606992", "1.541246", "1.488819")
)

# make them table 
am_BCPE_param_table <- ggtexttable(am_BCPE_param, rows = NULL, theme = ttheme("light", base_size = 12)) 

am_BCPE_param_table  <- table_cell_font(am_BCPE_param_table , row = 2:9, column = 1,
                       face = "bold")

am_BCPE_param_table %>% 
  tab_add_vline(at.column = 1, column.side = "right") 
  

##### for mean and sigma only 
x <- seq(min(am_mean_and_sigma_BCPE$y, na.rm = TRUE),
         max(am_mean_and_sigma_BCPE$y, na.rm = TRUE))

am_meansigma_BCPE_A_pdf  <- dBCPE(x, mu = 9.088002, sigma = 0.4219519, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_B_pdf  <- dBCPE(x, mu = 8.061705, sigma = 0.4417238, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_C_pdf  <- dBCPE(x, mu = 7.195985, sigma = 0.4522838, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_D_pdf  <- dBCPE(x, mu = 6.310598, sigma = 0.4713902, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_E_pdf  <- dBCPE(x, mu = 5.449665, sigma = 0.4843436, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_F_pdf  <- dBCPE(x, mu = 4.584880, sigma = 0.5007778, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_G_pdf  <- dBCPE(x, mu = 3.729521, sigma = 0.5150814, nu = 0.3808586, tau = 1.536581)
am_meansigma_BCPE_H_pdf  <- dBCPE(x, mu = 2.826814, sigma = 0.5356297, nu = 0.3808586, tau = 1.536581)


# 5. plot that sucka to compare model prediction and non-parametric estimate ---------------------------------------------------

par(mfrow = c(2, 2))
# plot density plots for period A ------------------------------
plot(Ammonium_A_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     main = "A", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_A_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_A_pdf, type = "l", lwd = 2, col = "steelblue") 

# plot density plots for period B ------------------------------
plot(Ammonium_B_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     main = "B", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_B_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_B_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period C ------------------------------
plot(Ammonium_C_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.142), 
     main = "C", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_C_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_C_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period D ------------------------------
plot(Ammonium_D_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.16),
     main = "D", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_D_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_D_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period E ------------------------------
plot(Ammonium_E_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.175),
     main = "E", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_E_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_E_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period F ------------------------------
plot(Ammonium_F_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.2),
     main = "F", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_F_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_F_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period G ------------------------------
plot(Ammonium_G_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.25),
     main = "G", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_G_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_G_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period H ------------------------------
plot(Ammonium_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", ylim = c(0,0.3), 
     main = "H", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_H_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_G_pdf, type = "l", lwd = 2, col = "steelblue") 

# reset 
par(mfrow = c(1, 1))


###### tail emphasis 
plot(Ammonium_A_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     xlim = c(16,22), ylim = c(0,0.05),
     main = "A", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_A_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_A_pdf, type = "l", lwd = 2, col = "steelblue") 

# plot density plots for period B ------------------------------
plot(Ammonium_B_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     xlim = c(16,35), ylim = c(0,0.01),
     main = "B", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_B_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_B_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period C ------------------------------
plot(Ammonium_C_pdf, lwd = 2, lty = 2, col = "black", type = "l",
     xlim = c(16,30), ylim = c(0,0.01),
     main = "C", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_C_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_C_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period D ------------------------------
plot(Ammonium_D_pdf, lwd = 2, lty = 2, col = "black", type = "l",
     xlim = c(14,18), ylim = c(0,0.04),
     main = "D", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_D_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_D_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period E ------------------------------
plot(Ammonium_E_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     xlim = c(13,18), ylim = c(0,0.015),
     main = "E", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_E_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_E_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period F ------------------------------
plot(Ammonium_F_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     xlim = c(10,25), ylim = c(0,0.03),
     main = "F", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_F_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_F_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period G ------------------------------
plot(Ammonium_G_pdf, lwd = 2, lty = 2, col = "black", type = "l",
     xlim = c(10,20), ylim = c(0,0.02),
     main = "G", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_G_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_G_pdf, type = "l", lwd = 2, col = "steelblue") 


# plot density plots for period H ------------------------------
plot(Ammonium_H_pdf, lwd = 2, lty = 2, col = "black", type = "l", 
     xlim = c(8,11), ylim = c(0,0.04),
     main = "H", xlab = "Ammonium (µmol/l)")
lines(am_all_BCPE_sigmapoly_H_pdf, type = "l", lwd = 2, col = "goldenrod") 
lines(am_meansigma_BCPE_G_pdf, type = "l", lwd = 2, col = "steelblue") 

########


# METHOD 1 ===========
logSurv(Ammonium$Ammonium, prob=0.90, tail="right")


# METHOD 2 ===========
par(mfrow = c(1, 3))
am_m1 <- loglogSurv1(Ammonium$Ammonium, prob=0.90, title="(a) TYPE I")
am_m2 <- loglogSurv2(Ammonium$Ammonium, prob=0.90, title="(b) TYPE II")
am_m3 <- loglogSurv3(Ammonium$Ammonium, prob=0.90, title="(c) TYPE III")
par(mfrow = c(1, 1))


# METHOD 3 ===========

# use some truncated family distribution 
am_m4 <- fitTail(Ammonium$Ammonium, family=WEI, percentage=10)
am_m5 <- fitTail(Ammonium$Ammonium, family=LOGNO, percentage=10)
am_m6 <- fitTail(Ammonium$Ammonium, family=BCPE, percentage=10)
am_m7 <- fitTail(Ammonium$Ammonium, family=BCT, percentage=10)

# check AIC 
AIC(am_m4, am_m5, am_m6, am_m7)
## NOTE: BCPE does the best

wp(am_m6, ylim.all = 1)
am_m6_2 <- fitTailAll(Ammonium$Ammonium, family=BCPE) # got like 50 warnings
plot(am_m6_2)


# METHOD 4 ===========
am_f1 <- fitDist(Ammonium$Ammonium, ncpus=4, parallel="snow" )
am_f1$fits[1:5]
wp(am_f1, ylim.all=1.5)

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
total_obs <- length(Ammonium$Ammonium)
print(paste("Total observations:", total_obs))

# 2. Number of observations in the top 10%
num_top_10_percent <- ceiling(total_obs * 0.10)
print(paste("Number of observations in the top 10%:", num_top_10_percent))




# empirical estimate of moments
am_moments_summary <- Ammonium %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Ammonium, na.rm = TRUE),
    var = var(Ammonium, na.rm = TRUE),
    skew = skewness(Ammonium, na.rm = TRUE),
    kurt = kurtosis(Ammonium, na.rm = TRUE)
  )
print(am_moments_summary, n=33)










library(ggridges)
library(gamlss)
# adding smoothing spline 

am_all_BCPE_spline <- gamlss(Ammonium ~ pb(year) + month, sigma.fo = ~ pb(year) + month, nu.fo = ~ pb(year) + month, tau.fo = ~ pb(year) + month,
                               family = BCPE(), data = Ammonium,
                               method = mixed(5,1000),
                             control = gamlss.control(n.cyc = 1000, c.crit = 0.01, trace = TRUE))

am_all_BCPE_spline_param <- best_model_param(am_all_BCPE_spline, Ammonium)

set.seed(123)
am_spline_predict <- list(
  A = rBCPE(656, mu = 9.272481, sigma = 0.3048066, nu = 0.8616693, tau = 2.307871),
  B = rBCPE(532, mu = 8.167273, sigma = 0.3562662, nu = 0.7774261, tau = 2.187535),
  C = rBCPE(562, mu = 6.608230, sigma = 0.5546833, nu = 0.7297309, tau = 2.062586),
  D = rBCPE(862, mu = 7.935265, sigma = 0.3463028, nu = 0.6610458, tau = 1.925822),
  E = rBCPE(959, mu = 5.181854, sigma = 0.4927180, nu = 0.6189171, tau = 1.833270),
  F = rBCPE(980, mu = 4.467627, sigma = 0.3929425, nu = 0.5622363, tau = 1.739978),
  G = rBCPE(973, mu = 3.443021, sigma = 0.5122292, nu = 0.5092856, tau = 1.640085),
  H = rBCPE(968, mu = 2.993934, sigma = 0.4600001, nu = 0.4433367, tau = 1.557013)
)


am_spline_predict_long <- stack(am_spline_predict)
names(am_spline_predict_long) <- c("pred", "period")


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Ammonium, aes(x = Ammonium, y = period), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = am_spline_predict_long, aes(x = pred, y = period), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  




library(gamlss.tr)
gen.trun(0,"SEP4",type="left")
gen.trun(0,"SHASH",type="left")

##### use bins in the model 
Ammonium$period <- as.numeric(Ammonium$period)

am_all_GB2_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                  family = GB2(), data = Ammonium,
                                  mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                  method = mixed(5,100),
                                  control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


am_all_SEP4tr_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                 family = SEP4tr(), data = Ammonium,
                                 mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                 method = mixed(5,100),
                                 control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


am_all_GB2_bins_spline_param <- best_model_param(am_all_GB2_bins_spline, Ammonium)

am_all_SEP4tr_bins_spline_param <- best_model_param(am_all_SEP4tr_bins_spline, Ammonium)


set.seed(123)
am_all_SEP4tr_bins_spline_predict <- list(
  "1" = rSEP4tr(656, mu = 9.121592, sigma = 4.2968953, nu = 3.819410, tau = 1.6051666),
  "2" = rSEP4tr(532, mu = 8.412779, sigma = 4.6619553, nu = 3.945809, tau = 1.5906259),
  "3" = rSEP4tr(562, mu = 7.740065, sigma = 4.6636639, nu = 3.549853, tau = 1.6358919),
  "4" = rSEP4tr(862, mu = 6.815223, sigma = 3.7858888, nu = 3.208364, tau = 1.5080143),
  "5" = rSEP4tr(959, mu = 5.431378, sigma = 2.7543580, nu = 2.996246, tau = 1.2678445),
  "6" = rSEP4tr(980, mu = 4.288864, sigma = 1.7498774, nu = 2.752936, tau = 1.0919705),
  "7" = rSEP4tr(973, mu = 3.352023, sigma = 1.2482654, nu = 1.942267, tau = 0.9321908),
  "8" = rSEP4tr(968, mu = 2.966985, sigma = 0.9349161, nu = 1.084364, tau = 0.8724551)
)


am_all_SEP4tr_bins_spline_predict_long <- stack(am_all_SEP4tr_bins_spline_predict)
names(am_all_SEP4tr_bins_spline_predict_long) <- c("pred", "period")
am_all_SEP4tr_bins_spline_predict_long$period <- as.numeric(am_all_SEP4tr_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Ammonium, aes(x = Ammonium, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = am_all_SEP4tr_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Ammonium (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  
  



## bins, time-varying parameters 

am_all_SEP4tr_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                    family = SEP4tr(), data = Ammonium,
                                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                    method = mixed(5,100),
                                    control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


am_nokurt_SEP4tr_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ 1,
                                    family = SEP4tr(), data = Ammonium,
                                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                    method = mixed(5,100),
                                    control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


am_noskew_SEP4tr_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ cs(period) + month,
                                    family = SEP4tr(), data = Ammonium,
                                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                    method = mixed(5,100),
                                    control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


am_mean_sigma_SEP4tr_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ 1,
                                    family = SEP4tr(), data = Ammonium,
                                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                    method = mixed(5,100),
                                    control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


am_meanonly_SEP4tr_bins_spline <- gamlss(Ammonium ~ cs(period) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                                    family = SEP4tr(), data = Ammonium,
                                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                    method = mixed(5,100),
                                    control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


AIC(am_all_SEP4tr_bins_spline, am_nokurt_SEP4tr_bins_spline,
    am_noskew_SEP4tr_bins_spline, am_mean_sigma_SEP4tr_bins_spline,
    am_meanonly_SEP4tr_bins_spline)


### RESULTS 
# am_all_SEP4tr_bins_spline         28821.53
# am_noskew_SEP4tr_bins_spline      28871.96
# am_nokurt_SEP4tr_bins_spline      29064.76
# am_mean_sigma_SEP4tr_bins_spline  29113.06
# am_meanonly_SEP4tr_bins_spline    29756.22




am_intercept_SEP2 <- gamlss(Ammonium ~ 1,
                                         family = SEP2(), data = Ammonium,
                                         mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                         method = mixed(5,100),
                                         control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

set.seed(52)
moment_bucket(am_intercept_SEP2) +
resid_wp(am_intercept_SEP2) + theme(plot.title = element_blank())














# monte carlo simulation 

am_all_SEP4tr_bins_spline_predict <- list(
  "1" = rSEP4tr(656, mu = 9.121592, sigma = 4.2968953, nu = 3.819410, tau = 1.6051666),
  "2" = rSEP4tr(532, mu = 8.412779, sigma = 4.6619553, nu = 3.945809, tau = 1.5906259),
  "3" = rSEP4tr(562, mu = 7.740065, sigma = 4.6636639, nu = 3.549853, tau = 1.6358919),
  "4" = rSEP4tr(862, mu = 6.815223, sigma = 3.7858888, nu = 3.208364, tau = 1.5080143),
  "5" = rSEP4tr(959, mu = 5.431378, sigma = 2.7543580, nu = 2.996246, tau = 1.2678445),
  "6" = rSEP4tr(980, mu = 4.288864, sigma = 1.7498774, nu = 2.752936, tau = 1.0919705),
  "7" = rSEP4tr(973, mu = 3.352023, sigma = 1.2482654, nu = 1.942267, tau = 0.9321908),
  "8" = rSEP4tr(968, mu = 2.966985, sigma = 0.9349161, nu = 1.084364, tau = 0.8724551)
)

library(gamlss.dist)
library(moments)

############################## bin 1
mu_am_bin1   <-     9.121592
sigma_am_bin1 <-    4.2968953
nu_am_bin1  <-    3.819410
tau_am_bin1  <- 1.6051666

set.seed(1)
am_bin1_sim <- rSEP4tr(N, mu = mu_am_bin1, sigma = sigma_am_bin1, nu = nu_am_bin1, tau = tau_am_bin1)

# get skewness and kurtosis 
kurtosis(am_bin1_sim) - 3 # 0.3695847
skewness(am_bin1_sim) # 0.6147768
var(am_bin1_sim) # 9.160423
mean(am_bin1_sim) # 9.40857


############################# bin 2
mu_am_bin2   <-     8.412779
sigma_am_bin2 <-    4.6619553
nu_am_bin2  <-    3.945809
tau_am_bin2  <- 1.5906259
  
set.seed(2)
am_bin2_sim <- rSEP4tr(N, mu = mu_am_bin2, sigma = sigma_am_bin2, nu = nu_am_bin2, tau = tau_am_bin2)

# get skewness and kurtosis 
kurtosis(am_bin2_sim) - 3 # 0.3976929
skewness(am_bin2_sim) # 0.6325068
var(am_bin2_sim) # 10.87824
mean(am_bin2_sim) # 8.744228


############################# bin 3
mu_am_bin3   <-     7.740065
sigma_am_bin3 <-    4.6636639
nu_am_bin3  <-    3.549853
tau_am_bin3  <- 1.6358919
  
set.seed(3)
am_bin3_sim <- rSEP4tr(N, mu = mu_am_bin3, sigma = sigma_am_bin3, nu = nu_am_bin3, tau = tau_am_bin3)

# get skewness and kurtosis 
kurtosis(am_bin3_sim) - 3 # 0.2732945
skewness(am_bin3_sim) # 0.5629908
var(am_bin3_sim) # 10.6712
mean(am_bin3_sim) # 8.031593


############################# bin 4
mu_am_bin4   <-  6.815223   
sigma_am_bin4 <-    3.7858888
nu_am_bin4  <-    3.208364
tau_am_bin4  <- 1.5080143
  
set.seed(4)
am_bin4_sim <- rSEP4tr(N, mu = mu_am_bin4, sigma = sigma_am_bin4, nu = nu_am_bin4, tau = tau_am_bin4)

# get skewness and kurtosis 
kurtosis(am_bin4_sim) - 3 # 0.5864298
skewness(am_bin4_sim) # 0.6518984
var(am_bin4_sim) # 7.757278
mean(am_bin4_sim) # 7.117409


############################# bin 5
mu_am_bin5   <-     5.431378
sigma_am_bin5 <-    2.7543580
nu_am_bin5  <-    2.996246
tau_am_bin5  <- 1.2678445
  
set.seed(5)
am_bin5_sim <- rSEP4tr(N, mu = mu_am_bin5, sigma = sigma_am_bin5, nu = nu_am_bin5, tau = tau_am_bin5)

# get skewness and kurtosis 
kurtosis(am_bin5_sim) - 3 # 1.635682
skewness(am_bin5_sim) # 0.968011
var(am_bin5_sim) # 5.269674
mean(am_bin5_sim) # 5.811096


############################# bin 6
mu_am_bin6   <-     4.288864
sigma_am_bin6 <-    1.7498774
nu_am_bin6  <-    2.752936
tau_am_bin6  <- 1.0919705
  
set.seed(6)
am_bin6_sim <- rSEP4tr(N, mu = mu_am_bin6, sigma = sigma_am_bin6, nu = nu_am_bin6, tau = tau_am_bin6)

# get skewness and kurtosis 
kurtosis(am_bin6_sim) - 3 # 3.264097
skewness(am_bin6_sim) # 1.337953
var(am_bin6_sim) # 2.868245
mean(am_bin6_sim) # 4.674122


############################# bin 7
mu_am_bin7   <-     3.352023
sigma_am_bin7 <-   1.2482654 
nu_am_bin7  <-    1.942267
tau_am_bin7  <- 0.9321908
  
set.seed(7)
am_bin7_sim <- rSEP4tr(N, mu = mu_am_bin7, sigma = sigma_am_bin7, nu = nu_am_bin7, tau = tau_am_bin7)

# get skewness and kurtosis 
kurtosis(am_bin7_sim) - 3 # 5.046296
skewness(am_bin7_sim) # 1.640126
var(am_bin7_sim) # 2.340453
mean(am_bin7_sim) # 3.770817


############################# bin 8
mu_am_bin8   <-     2.966985
sigma_am_bin8 <-    0.9349161
nu_am_bin8  <-    1.084364
tau_am_bin8  <- 0.8724551
  
set.seed(8)
am_bin8_sim <- rSEP4tr(N, mu = mu_am_bin8, sigma = sigma_am_bin8, nu = nu_am_bin8, tau = tau_am_bin8)

# get skewness and kurtosis 
kurtosis(am_bin8_sim) - 3 # 4.623925
skewness(am_bin8_sim) # 1.339615
var(am_bin8_sim) # 1.951983
mean(am_bin8_sim) # 3.225493



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
am_bin1_sim_se <- bootstrap_moments_sim_only(am_bin1_sim) 
am_bin2_sim_se <- bootstrap_moments_sim_only(am_bin2_sim)
am_bin3_sim_se <- bootstrap_moments_sim_only(am_bin3_sim)
am_bin4_sim_se <- bootstrap_moments_sim_only(am_bin4_sim)
am_bin5_sim_se <- bootstrap_moments_sim_only(am_bin5_sim)
am_bin6_sim_se <- bootstrap_moments_sim_only(am_bin6_sim)
am_bin7_sim_se <- bootstrap_moments_sim_only(am_bin7_sim)
am_bin8_sim_se <- bootstrap_moments_sim_only(am_bin8_sim)

am_bin8_sim_se$se

am_bins_df <- tibble::tibble(
  year_mid = c(1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  mean = c(9.409,8.744,8.032,7.117,5.811,4.674,3.771,3.225),
  variance = c(9.160,10.878,10.671,7.757,5.270,2.868,2.340,1.952),
  skewness = c(0.610,0.633,0.563,0.652,0.968,1.338,1.640,1.340),
  kurtosis = c(0.370,0.398,0.273,0.586,1.636,3.264,5.046,4.623)
)

# Suppose you already have standard errors (one per year per stat).
# Example small SE values — replace with your actual SE or bootstrap-based SE
# If you have CI directly, skip this step and create lo/hi columns instead.
am_se_df <- tibble::tibble(
  year_mid = am_bins_df$year_mid,
  se_mean = c(0.004436985,0.004564469,0.004453601,0.003951382,
              0.003271781,0.002427279,0.002168838,0.001981363),
  se_variance = c(0.020178283,0.024048321,0.022752320,0.017471896,
                  0.014151229,0.009448144,0.008731005,0.007057966),
  se_skew = c(0.003876006,0.004041878,0.003716181,0.004287422,
              0.006093326,0.009030545,0.011670575,0.011669906),
  se_kurt = c(0.013622697,0.015933674,0.012790616,0.018134935,
              0.034477138,0.067954178,0.108433020,0.093706115)
)

# join and pivot to long form
am_df_long <- am_bins_df %>%
  left_join(am_se_df, by = "year_mid") %>%
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


ggplot(am_df_long, aes(x = year_mid, y = value)) +
  geom_line(size = 0.6) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.6) +   # vertical error bars
  facet_wrap(~stat, scales = "free_y", ncol = 1) +
  theme_minimal(base_size = 14) +
  labs(x = "Year (midpoint)", y = NULL, title = "Distribution moments through time") +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 12),
    plot.title = element_text(size = 18, hjust = 0.02)
  )





