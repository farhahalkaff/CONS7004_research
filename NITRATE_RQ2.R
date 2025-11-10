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
Nitrate$mu_hat_all_JSU_ni <- predict(ni_all_JSU, what = "mu", type = "response")


ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "grey", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_BCTo_ni), color = "goldenrod", linewidth = 1) +
  geom_line(aes(y = mu_hat_all_JSU_ni), color = "steelblue", linewidth = 1) +
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()

m <- fitTail(Nitrate$Nitrate, family=BCTo, percentage=10)
plot(m)
m1 <- fitTail(Nitrate$Nitrate, family=JSU, percentage=10)
plot(m1)

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



################

# METHOD 1 ======
logSurv(Nitrate$Nitrate, prob=0.90, tail="right")

# METHOD 2 ======
par(mfrow = c(1, 3))
ni_m1 <- loglogSurv1(Nitrate$Nitrate, prob=0.90, title="(a) TYPE I")
ni_m2 <- loglogSurv2(Nitrate$Nitrate, prob=0.90, title="(b) TYPE II")
ni_m3 <- loglogSurv3(Nitrate$Nitrate, prob=0.90, title="(c) TYPE III")
par(mfrow = c(1, 1))


# METHOD 3 ======

# use some truncated family distribution 
ni_m4 <- fitTail(Nitrate$Nitrate, family=WEI, percentage=10)
ni_m5 <- fitTail(Nitrate$Nitrate, family=LOGNO, percentage=10)
ni_m6 <- fitTail(Nitrate$Nitrate, family=BCPE, percentage=10)
ni_m7 <- fitTail(Nitrate$Nitrate, family=BCT, percentage=10)
ni_m8 <- fitTail(Nitrate$Nitrate, family=BCTo, percentage=10)

# check AIC 
AIC(ni_m4, ni_m5, ni_m6, ni_m7, ni_m8)
## NOTE: BCPE does the best

wp(ni_m8, ylim.all = 1)
am_m6_2 <- fitTailAll(Ammonium$Ammonium, family=BCPE) # got like 50 warnings
plot(am_m6_2)



# 1. Total number of observations
total_obs <- length(Nitrate$Nitrate)
print(paste("Total observations:", total_obs))

# 2. Number of observations in the top 10%
num_top_10_percent <- ceiling(total_obs * 0.10)
print(paste("Number of observations in the top 10%:", num_top_10_percent))





# empirical estimate of moments
ni_moments_summary <- Nitrate %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Nitrate, na.rm = TRUE),
    var = var(Nitrate, na.rm = TRUE),
    skew = skewness(Nitrate, na.rm = TRUE),
    kurt = kurtosis(Nitrate, na.rm = TRUE)
  )
print(ni_moments_summary, n=33)







ni_all_BCPE_meanpoly <- gamlss(Nitrate ~ pb(year) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCPE(), data = Nitrate,
                      method = mixed(5,100),
                      control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
 

ni_all_BCPEo_spline_param <- best_model_param(ni_all_BCPE_meanpoly, Nitrate)

set.seed(123)
ni_spline_predict <- list(
  A = rBCPE(656, mu = 9.010212, sigma = 0.6572113, nu = 0.2998671, tau = 2.567016),
  B = rBCPE(532, mu = 8.343500, sigma = 0.6339938, nu = 0.3357239, tau = 2.383719),
  C = rBCPE(562, mu = 9.573813, sigma = 0.6131986, nu = 0.3550464, tau = 2.233479),
  D = rBCPE(862, mu = 9.981417, sigma = 0.5985880, nu = 0.3776479, tau = 2.086386),
  E = rBCPE(959, mu = 13.956159, sigma = 0.5791233, nu = 0.3999346, tau = 1.946922),
  F = rBCPE(980, mu = 15.491823, sigma = 0.5619768, nu = 0.4250580, tau = 1.817000),
  G = rBCPE(973, mu = 24.138190, sigma = 0.5468654, nu = 0.4470378, tau = 1.707011),
  H = rBCPE(968, mu = 24.470387, sigma = 0.5302117, nu = 0.4717803, tau = 1.596139)
)


ni_spline_predict_long <- stack(ni_spline_predict)
names(ni_spline_predict_long) <- c("pred", "period")


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Nitrate, aes(x = Nitrate, y = period), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = ni_spline_predict_long, aes(x = pred, y = period), 
                      color = "blue4", alpha = 0.1) +
  theme_classic()  





##### use bins in the model 
Nitrate$period <- as.numeric(Nitrate$period)

library(gamlss.tr)
gen.trun(0,"SEP2",type="left")
gen.trun(0,"SHASHo2",type="left")

ni_all_SHASHo2tr_bins_spline <- gamlss(Nitrate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                   family = SHASHo2tr(), data = Nitrate,
                                   mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                   #method = mixed(5,100),
                                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ni_all_SHASHo2tr_bins_spline_noseason <- gamlss(Nitrate ~ cs(period), sigma.fo = ~ cs(period), nu.fo = ~ cs(period), tau.fo = ~ cs(period),
                                       family = SHASHo2tr(), data = Nitrate,
                                       mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                       #method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 300, c.crit = 0.01, trace = TRUE))


ni_all_SHASHo2tr_bins_spline_param <- best_model_param(ni_all_SHASHo2tr_bins_spline, Nitrate)
ni_all_SHASHo2tr_bins_noseason_spline_param <- best_model_param(ni_all_SHASHo2tr_bins_spline_noseason, Nitrate)

## with seasonaility 
set.seed(123)
ni_all_SHASHo2tr_bins_spline_predict <- list(
  "1" = rSHASHo2tr(656, mu = 5.174920, sigma = 5.522785, nu = 0.05077909, tau = 0.7664181),
  "2" = rSHASHo2tr(532, mu = 4.851424, sigma = 5.022379, nu = 0.05085468, tau = 0.6801695),
  "3" = rSHASHo2tr(562, mu = 5.422420, sigma = 4.644028, nu = 0.12158318, tau = 0.7085193),
  "4" = rSHASHo2tr(862, mu = 6.595169, sigma = 4.932488, nu = 0.32034797, tau = 0.8000461),
  "5" = rSHASHo2tr(959, mu = 8.285170, sigma = 5.982398, nu = 0.61810835, tau = 1.0212367),
  "6" = rSHASHo2tr(980, mu = 9.330573, sigma = 7.517602, nu = 0.92258305, tau = 1.3682099),
  "7" = rSHASHo2tr(973, mu = 8.958703, sigma = 9.321589, nu = 1.16271996, tau = 1.3711291),
  "8" = rSHASHo2tr(968, mu = 6.611717, sigma = 10.001079, nu = 1.15990392, tau = 0.9788938)
)


ni_all_SHASHo2tr_bins_spline_predict_long <- stack(ni_all_SHASHo2tr_bins_spline_predict)
names(ni_all_SHASHo2tr_bins_spline_predict_long) <- c("pred", "period")
ni_all_SHASHo2tr_bins_spline_predict_long$period <- as.numeric(ni_all_SHASHo2tr_bins_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Nitrate, aes(x = Nitrate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = ni_all_SHASHo2tr_bins_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Nitrate (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  
  

## without seasonaility 
set.seed(123)
ni_all_SHASHo2tr_bins_noseason_spline_predict <- list(
  "1" = rSHASHo2tr(656, mu = , sigma = , nu = , tau = ),
  "2" = rSHASHo2tr(532, mu = , sigma = , nu = , tau = ),
  "3" = rSHASHo2tr(562, mu = , sigma = , nu = , tau = ),
  "4" = rSHASHo2tr(862, mu = , sigma = , nu = , tau = ),
  "5" = rSHASHo2tr(959, mu = , sigma = , nu = , tau = ),
  "6" = rSHASHo2tr(980, mu = , sigma = , nu = , tau = ),
  "7" = rSHASHo2tr(973, mu = , sigma = , nu = , tau = ),
  "8" = rSHASHo2tr(968, mu = , sigma = , nu = , tau = )
)


ni_all_SHASHo2tr_bins_noseason_spline_predict_long <- stack(ni_all_SHASHo2tr_bins_noseason_spline_predict)
names(ni_all_SHASHo2tr_bins_noseason_spline_predict_long) <- c("pred", "period")
ni_all_SHASHo2tr_bins_noseason_spline_predict_long$period <- as.numeric(ni_all_SHASHo2tr_bins_noseason_spline_predict_long$period)


#### create ridge plot to compare non-parametric density plot and model's prediction
ggplot() +
  geom_density_ridges(data = Nitrate, aes(x = Nitrate, y = as.factor(period)), 
                      color = "azure4", alpha = 0.35) +
  geom_density_ridges(data = ni_all_SHASHo2tr_bins_noseason_spline_predict_long, aes(x = pred, y = as.factor(period)), 
                      color = "blue4", alpha = 0.1) +
  labs(x = "Nitrate (µmol/l)", y = "Years") +
  scale_y_discrete(labels = c("1962 - '66", "1967 - '70", "1971 - '74", "1975 - '78",
                              "1979 - '82", "1983 - '86", "1987 - '90", "1991 - '94")) +
  theme_classic()  










## Bins, time-varying parameters 

ni_all_SHASHo2tr_bins_spline <- gamlss(Nitrate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ cs(period) + month,
                                       family = SHASHo2tr(), data = Nitrate,
                                       mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                       #method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


ni_nokurt_SHASHo2tr_bins_spline <- gamlss(Nitrate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ cs(period) + month, tau.fo = ~ 1,
                                       family = SHASHo2tr(), data = Nitrate,
                                       mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                       #method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ni_noskew_SHASHo2tr_bins_spline <- gamlss(Nitrate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ cs(period) + month,
                                       family = SHASHo2tr(), data = Nitrate,
                                       mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                       #method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ni_mean_sigma_SHASHo2tr_bins_spline <- gamlss(Nitrate ~ cs(period) + month, sigma.fo = ~ cs(period) + month, nu.fo = ~ 1, tau.fo = ~ 1,
                                       family = SHASHo2tr(), data = Nitrate,
                                       mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                       #method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = TRUE))


ni_meanonly_SHASHo2tr_bins_spline <- gamlss(Nitrate ~ cs(period) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                                       family = SHASHo2tr(), data = Nitrate,
                                       mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                       #method = mixed(5,100),
                                       control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


AIC(ni_all_SHASHo2tr_bins_spline, ni_nokurt_SHASHo2tr_bins_spline,
    ni_noskew_SHASHo2tr_bins_spline, ni_mean_sigma_SHASHo2tr_bins_spline,
    ni_meanonly_SHASHo2tr_bins_spline)



## RESULTS 
# ni_all_SHASHo2tr_bins_spline         42846.14
# ni_noskew_SHASHo2tr_bins_spline      43135.68
# ni_nokurt_SHASHo2tr_bins_spline      43146.61
# ni_mean_sigma_SHASHo2tr_bins_spline  43293.13
# ni_meanonly_SHASHo2tr_bins_spline    44707.54





library(gamlss.tr)
gen.trun(family = "SEP2", type = "left")

ni_intercept_SEP2 <- gamlss(Nitrate ~  1,
                                            family = SEP2(), data = Nitrate,
                                            mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                                            #method = mixed(5,100),
                                            control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))

ni_intercept_SEP2tr <- gamlss(Nitrate ~  1,
                            family = SEP2tr(), data = Nitrate,
                            mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                            #method = mixed(5,100),
                            control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


set.seed(52)
moment_bucket(ni_intercept_SEP2) +
resid_wp(ni_intercept_SEP2) + theme(plot.title = element_blank())







#### monte carlo simulation 

ni_all_SHASHo2tr_bins_spline_predict <- list(
  "1" = rSHASHo2tr(656, mu = 5.174920, sigma = 5.522785, nu = 0.05077909, tau = 0.7664181),
  "2" = rSHASHo2tr(532, mu = 4.851424, sigma = 5.022379, nu = 0.05085468, tau = 0.6801695),
  "3" = rSHASHo2tr(562, mu = 5.422420, sigma = 4.644028, nu = 0.12158318, tau = 0.7085193),
  "4" = rSHASHo2tr(862, mu = 6.595169, sigma = 4.932488, nu = 0.32034797, tau = 0.8000461),
  "5" = rSHASHo2tr(959, mu = 8.285170, sigma = 5.982398, nu = 0.61810835, tau = 1.0212367),
  "6" = rSHASHo2tr(980, mu = 9.330573, sigma = 7.517602, nu = 0.92258305, tau = 1.3682099),
  "7" = rSHASHo2tr(973, mu = 8.958703, sigma = 9.321589, nu = 1.16271996, tau = 1.3711291),
  "8" = rSHASHo2tr(968, mu = 6.611717, sigma = 10.001079, nu = 1.15990392, tau = 0.9788938)
)


library(gamlss.dist)
library(moments)

############################## bin 1
mu_ni_bin1   <- 5.174920    
sigma_ni_bin1 <- 5.522785   
nu_ni_bin1  <- 0.05077909   
tau_ni_bin1  <- 0.7664181

set.seed(1)
N <- 5e5                    
ni_bin1_sim <- rSHASHo2tr(N, mu = mu_ni_bin1, sigma = sigma_ni_bin1, nu = nu_ni_bin1, tau = tau_ni_bin1)

# get skewness and kurtosis 
kurtosis(ni_bin1_sim) - 3 # 1.46827
skewness(ni_bin1_sim) # 1.05859
var(ni_bin1_sim) # 26.60115
mean(ni_bin1_sim) # 7.596501


############################# bin 2
mu_ni_bin2   <- 4.851424    
sigma_ni_bin2 <- 5.022379   
nu_ni_bin2  <- 0.05085468   
tau_ni_bin2  <- 0.6801695

set.seed(2)
ni_bin2_sim <- rSHASHo2tr(N, mu = mu_ni_bin2, sigma = sigma_ni_bin2, nu = nu_ni_bin2, tau = tau_ni_bin2)

# get skewness and kurtosis 
kurtosis(ni_bin2_sim) - 3 # 2.659685
skewness(ni_bin2_sim) # 1.329857
var(ni_bin2_sim) # 27.23116
mean(ni_bin2_sim) # 7.332844


############################# bin 3
mu_ni_bin3   <- 5.422420    
sigma_ni_bin3 <- 4.644028   
nu_ni_bin3  <- 0.12158318   
tau_ni_bin3  <- 0.7085193

set.seed(3)
ni_bin3_sim <- rSHASHo2tr(N, mu = mu_ni_bin3, sigma = sigma_ni_bin3, nu = nu_ni_bin3, tau = tau_ni_bin3)

# get skewness and kurtosis 
kurtosis(ni_bin3_sim) - 3 # 2.185784
skewness(ni_bin3_sim) # 1.212528
var(ni_bin3_sim) # 25.93133
mean(ni_bin3_sim) # 7.640072


############################# bin 4
mu_ni_bin4   <- 6.595169    
sigma_ni_bin4 <- 4.932488   
nu_ni_bin4  <- 0.32034797   
tau_ni_bin4  <- 0.8000461

set.seed(4)
ni_bin4_sim <- rSHASHo2tr(N, mu = mu_ni_bin4, sigma = sigma_ni_bin4, nu = nu_ni_bin4, tau = tau_ni_bin4)

# get skewness and kurtosis 
kurtosis(ni_bin4_sim) - 3 # 1.600815
skewness(ni_bin4_sim) # 1.083545
var(ni_bin4_sim) # 34.72103
mean(ni_bin4_sim) # 9.586785


############################# bin 5
mu_ni_bin5   <- 8.285170 
sigma_ni_bin5 <- 5.982398   
nu_ni_bin5  <-  0.61810835  
tau_ni_bin5  <- 1.0212367

set.seed(5)
ni_bin5_sim <- rSHASHo2tr(N, mu = mu_ni_bin5, sigma = sigma_ni_bin5, nu = nu_ni_bin5, tau = tau_ni_bin5)

# get skewness and kurtosis 
kurtosis(ni_bin5_sim) - 3 # 0.6845814
skewness(ni_bin5_sim) # 0.8698108
var(ni_bin5_sim) # 51.81536
mean(ni_bin5_sim) # 13.57627


############################# bin 6
mu_ni_bin6   <-  9.330573   
sigma_ni_bin6 <-    7.517602
nu_ni_bin6  <-  0.92258305  
tau_ni_bin6  <- 1.3682099

set.seed(6)
ni_bin6_sim <- rSHASHo2tr(N, mu = mu_ni_bin6, sigma = sigma_ni_bin6, nu = nu_ni_bin6, tau = tau_ni_bin6)

# get skewness and kurtosis 
kurtosis(ni_bin6_sim) - 3 # -0.01217206
skewness(ni_bin6_sim) # 0.6464722
var(ni_bin6_sim) # 72.01458
mean(ni_bin6_sim) # 18.16769



############################# bin 7
mu_ni_bin7   <-   8.958703  
sigma_ni_bin7 <-  9.321589  
nu_ni_bin7  <-    1.16271996
tau_ni_bin7  <- 1.3711291

set.seed(7)
ni_bin7_sim <- rSHASHo2tr(N, mu = mu_ni_bin7, sigma = sigma_ni_bin7, nu = nu_ni_bin7, tau = tau_ni_bin7)

# get skewness and kurtosis 
kurtosis(ni_bin7_sim) - 3 # 0.09462743
skewness(ni_bin7_sim) # 0.7370327
var(ni_bin7_sim) # 139.5628
mean(ni_bin7_sim) # 23.33412


############################# bin 8
mu_ni_bin8   <-     6.611717
sigma_ni_bin8 <-    10.001079
nu_ni_bin8  <-    1.15990392
tau_ni_bin8  <- 0.9788938

set.seed(8)
ni_bin8_sim <- rSHASHo2tr(N, mu = mu_ni_bin8, sigma = sigma_ni_bin8, nu = nu_ni_bin8, tau = tau_ni_bin8)

# get skewness and kurtosis 
kurtosis(ni_bin8_sim) - 3 # 1.687261
skewness(ni_bin8_sim) # 1.280998
var(ni_bin8_sim) # 364.2776
mean(ni_bin8_sim) # 26.57479




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
ni_bin1_sim_se <- bootstrap_moments_sim_only(ni_bin1_sim) 
ni_bin2_sim_se <- bootstrap_moments_sim_only(ni_bin2_sim) 
ni_bin3_sim_se <- bootstrap_moments_sim_only(ni_bin3_sim) 
ni_bin4_sim_se <- bootstrap_moments_sim_only(ni_bin4_sim) 
ni_bin5_sim_se <- bootstrap_moments_sim_only(ni_bin5_sim) 
ni_bin6_sim_se <- bootstrap_moments_sim_only(ni_bin6_sim) 
ni_bin7_sim_se <- bootstrap_moments_sim_only(ni_bin7_sim) 
ni_bin8_sim_se <- bootstrap_moments_sim_only(ni_bin8_sim) 

ni_bin8_sim_se$se

ni_bins_df <- tibble::tibble(
  year_mid = c(1964, 1968, 1972, 1976, 1980, 1984, 1988, 1992),
  mean = c(7.597,7.333,7.640,9.587,13.576,18.168,23.334,26.575),
  variance = c(26.601,27.231,25.931,34.721,51.815,72.015,139.563,364.277),
  skewness = c(1.059,1.330,1.213,1.084,0.870,0.646,0.737,1.281),
  kurtosis = c(1.468,2.659,2.186,1.601,0.685,-0.012,0.095,1.687)
)

ni_se_df <- tibble::tibble(
  year_mid = ni_bins_df$year_mid,
  se_mean = c(0.007424578,0.007710453,0.006965149,0.008256720,
              0.010187013,0.011729352,0.016204944,0.026905917),
  se_variance = c(0.069718673,0.085008131,0.073731347,0.094024152,
                  0.120707255,0.138510137,0.286956477,0.995877110),
  se_skew = c(0.005153090,0.007194377,0.006110511,0.005324222,
              0.003817211,0.002978047,0.003031812,0.004510697),
  se_kurt = c(0.027114858,0.052560952,0.037299409,0.030475993,
              0.015516501,0.008719303,0.009396691,0.023431861)
)

# join and pivot to long form
ni_df_long <- ni_bins_df %>%
  left_join(ni_se_df, by = "year_mid") %>%
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


ggplot(ni_df_long, aes(x = year_mid, y = value)) +
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



