library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(gamlss)
library(DHARMa)

# Set the directory where your files are located
data_directory <- "dataset"

files <- list.files(path = data_directory, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("Helgoland_2002.csv", files)]
files <- files[!grepl("Helgoland_1998.csv", files)]

list_of_dataframes <- lapply(files, read.csv)

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

# create own df to remove NAs 
Nitrate <- df %>% filter(!is.na(Nitrate)) %>% 
  select(Date, Nitrate, year, month)

# density plot 
ggplot(Nitrate, aes(y = Nitrate)) +
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(Nitrate$Nitrate), linetype = "dashed", color = "azure4") +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


#============================================================================================ SST
# Intercept model 
niSSTm <- gamlss(Nitrate ~ 1, family = SST(), data = Nitrate,
                 mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                 method = mixed(10,200),
                 control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm)


# Only mean changing through time linearly 
niSSTm2 <- gamlss(Nitrate ~ Date, family = SST(), data = Nitrate,
                  #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                  method = mixed(10,200),
                  control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm2)


# mean, sigma and nu changing through time linearly 
niSSTm3 <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrate,
                  #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                  method = mixed(10,200),
                  control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm3)
AIC(niSSTm3)

# polynomial
niSSTpolym1 <- gamlss(Nitrate ~ poly(Date, 2), sigma.fo = ~ poly(Date, 2), nu.fo = ~ poly(Date, 2), family = SST(), data = Nitrate,
                     #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTpolym1)
AIC(niSSTpolym1)

# mu sig and nu changing with seasonality
niSSTm3_sea <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, family = SST(), data = Nitrate,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm3_sea)
AIC(niSSTm3_sea)

# only mean seasonaility
niSSTm2_sea <- gamlss(Nitrate ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrate,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm2_sea)
AIC(niSSTm2_sea)
#============================================================================================

# predict mu based on seasonaility model
Nitrate$mu_hat <- fitted(niSSTm3_sea, what = "mu")
Nitrate$sigma_hat <- fitted(niSSTm3_sea, what = "sigma", type = "response")
Nitrate$nu_hat <- fitted(niSSTm3_sea, what = "nu", type = "response")
Nitrate$tau_hat <- fitted(niSSTm3_sea, what = "tau", type = "response")


# plot that sucka (intercept model)
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "grey", na.rm = TRUE) +
  geom_line(aes(y = mu_hat), color = "darkolivegreen", linewidth = 1) +
  geom_abline(intercept = niSSTm2$mu.coefficients[1], slope = niSSTm2$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # only mean as function of time 
  geom_abline(intercept = niSSTm3$mu.coefficients[1], slope = niSSTm3$mu.coefficients[2], 
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu as function of time
  #geom_abline(intercept = niNOm1$mu.coefficients[1], slope = niNOm1$mu.coefficients[2], 
  #color = "darkolivegreen", linewidth = 1) + # mean, sigma and nu as function of time
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "salmon") +
  geom_hline(yintercept = 16.718, linetype = "dashed", color = "black") + # intercept mean
  geom_vline(xintercept = as.Date("1963-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1980-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1994-04-13"), linetype = "dashed", color = "darkred") +
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_classic()

#============================================================================================
# function for pdf 
pdf_SST_at_date <- function(model, data, date_str, x = NULL, n = 200,
                            time_var = "Date", date_var = "Date") {
  # 1. Find the row in df matching the date of interest
  t_val <- data[[time_var]][as.Date(data[[date_var]]) == as.Date(date_str)]
  if (length(t_val) == 0) stop("Date not found in data frame.")
  
  # 2. Build newdata with the correct covariate
  newdat <- data.frame(setNames(list(t_val), time_var))
  
  # 3. Predict distribution parameters on natural scale
  mu    <- predict(model, "mu",    newdata = newdat, type = "response")
  sigma <- predict(model, "sigma", newdata = newdat, type = "response")
  nu    <- predict(model, "nu",    newdata = newdat, type = "response")
  tau   <- predict(model, "tau",   newdata = newdat, type = "response")
  
  # 4. Create x grid if none supplied
  if (is.null(x)) {
    x <- seq(min(model$y, na.rm = TRUE),
             max(model$y, na.rm = TRUE),
             length.out = n)
  }
  
  # 5. Evaluate PDF
  dens <- dSST(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
  
  list(x = x, density = dens,
       params = c(mu = mu, sigma = sigma, nu = nu, tau = tau))
}

# get the pdf and param for three dates
# m1 (intercept model)
resm1 <- pdf_SST_at_date(niSSTm, Nitrate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm1b <- pdf_SST_at_date(niSSTm, Nitrate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm1c <- pdf_SST_at_date(niSSTm, Nitrate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m2 (mean only)
resm2 <- pdf_SST_at_date(niSSTm2, Nitrate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm2b <- pdf_SST_at_date(niSSTm2, Nitrate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm2c <- pdf_SST_at_date(niSSTm2, Nitrate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m3 (mean sigma and nu)
resm3 <- pdf_SST_at_date(niSSTm3, Nitrate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm3b <- pdf_SST_at_date(niSSTm3, Nitrate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm3c <- pdf_SST_at_date(niSSTm3, Nitrate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# pdf for seasonaility model
x <- seq(min(niSSTm3_sea$y, na.rm = TRUE),
         max(niSSTm3_sea$y, na.rm = TRUE),
         length.out = 200)

pdf_1963 <- dSST(x, mu = 1.2156707, sigma = 3.290147, nu = 9.663093, tau = 6.874041)
pdf_1980 <- dSST(x, mu = 8.671901, sigma = 6.657359, nu = 7.422281, tau = 6.874041)
pdf_1994 <- dSST(x, mu = 31.41314, sigma = 24.87263, nu = 0.9332162, tau = 6.874041)

#============================================================================================
# plot all three dates 

par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.15),
     ylab = "Density", xlab = "Value",
     main = "1963-09-23") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3$x, resm3$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(x, pdf_1963, col = "chocolate", lwd = 2, lty = 1)# seasonaility model
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.15),
     ylab = "Density", xlab = "Value",
     main = "1980-09-03") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3b$x, resm3b$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma an nu 
lines(x, pdf_1980, col = "chocolate", lwd = 2, lty = 1)# seasonaility model
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
# plot the models together for date 3: 1994-04-13
plot(resm1c$x, resm1c$density, col = "black", type = "l", lwd = 2, lty =2, ylim=c(0,0.15),
     ylab = "Density", xlab = "Value",
     main = "1994-04-13") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3c$x, resm3c$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(x, pdf_1994, col = "chocolate", lwd = 2, lty = 1)# seasonaility model
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
par(mfrow = c(1,1)) # reset
#============================================================================================

# m1 (intercept model)
resm1$params #1963
resm1b$params #1980
resm1c$params #1994
# m2 (mean only)
resm2$params #1963
resm2b$params #1980
resm2c$params #1994
# m3 (mean sigma and nu)
resm3$params #1963
resm3b$params #1980
resm3c$params #1994
# normal 
resm4$params #1963
resm4b$params #1980
resm4c$params #1994
# seasonaility 
param_1963 <- c(Nitrate[242,5], Nitrate[242,6], Nitrate[242,7], Nitrate[242,8])
param_1980 <- c(Nitrate[3025,5], Nitrate[3025,6], Nitrate[3025,7], Nitrate[3025,8])
param_1994 <- c(Nitrate[6323,5], Nitrate[6323,6], Nitrate[6323,7], Nitrate[6323,8])
param_1994 <- c(Nitrate[6435,5], Nitrate[6435,6], Nitrate[6435,7], Nitrate[6435,8])

param_table <- data.frame(
  models = c("Intercept", "Mean(time)", "All_param(time)", "Seasonality"),
  mu = c(resm1$params[1], resm2$params[1],resm3$params[1],  Nitrate[242,5]),
  sigma = c(resm1$params[2], resm2$params[2], resm3$params[2], Nitrate[242,6]),
  nu = c(resm1$params[3], resm2$params[3], resm3$params[3], Nitrate[242,7]),
  tau = c(resm1$params[4], resm2$params[4], resm3$params[4], Nitrate[242,8])
)


# look at seasoanility of same months different year 
# 1963
pdf_1963_Apr <- dSST(x, mu = 17.8164808, sigma = 6.879606, nu = 1.509802, tau = 6.874041) # April
pdf_1963_Aug <- dSST(x, mu = 0.8980581, sigma = 3.253148, nu = 9.926294, tau = 6.874041) # August
pdf_1963_Dec <- dSST(x, mu = 5.6005733, sigma = 4.221979, nu = 2.678465, tau = 6.874041) # December
# 1980
pdf_1980_Apr <- dSST(x, mu = 25.272711, sigma = 13.920353, nu = 1.159688, tau = 6.874041) # April
pdf_1980_Aug <- dSST(x, mu = 8.354289, sigma = 6.582495, nu = 7.624447, tau = 6.874041) # August
pdf_1980_Dec <- dSST(x, mu = 13.056804, sigma = 8.542849, nu = 2.057345, tau = 6.874041) # December
## 1994
pdf_1994_Apr <- dSST(x, mu = 31.41314, sigma = 24.87263, nu = 0.9332162, tau = 6.874041) # April
pdf_1994_Aug <- dSST(x, mu = 14.49471, sigma = 11.76148, nu = 6.1354915, tau = 6.874041) # August
pdf_1994_Dec <- dSST(x, mu = 19.19723, sigma = 15.26420, nu = 1.6555721, tau = 6.874041) # December


par(mfrow = c(1, 3)) 
# plot the models together for date month of April
plot(x, pdf_1963_Apr, col = "chocolate", type = "l", lwd = 2, ylim=c(0,0.14),
     ylab = "Density", xlab = "Value",
     main = "April") # 1963
lines(x, pdf_1980_Apr, col = "cornflowerblue", lwd = 2) # 1980
lines(x, pdf_1994_Apr, col = "darkseagreen4",lwd = 2) # 1994
legend("topright", legend = c("1963", "1980", "1994"),
       col = c("chocolate", "cornflowerblue", "darkseagreen4"), cex = 0.8, lty = c(1,1,1), bty = "n")
# plot the models together for month August
plot(x, pdf_1963_Aug, col = "chocolate", type = "l", lwd = 2, ylim=c(0,0.14),
      ylab = "Density", xlab = "Value",
      main = "August") # 1963
lines(x, pdf_1980_Aug, col = "cornflowerblue", lwd = 2) # 1980
lines(x, pdf_1994_Aug, col = "darkseagreen4",lwd = 2) # 1994
legend("topright", legend = c("1963", "1980", "1994"),
       col = c("chocolate", "cornflowerblue", "darkseagreen4"), cex = 0.8, lty = c(1,1,1), bty = "n")
# plot the models together for month December
plot(x, pdf_1963_Dec, col = "chocolate", type = "l", lwd = 2, ylim=c(0,0.14),
      ylab = "Density", xlab = "Value",
      main = "December") # 1963
lines(x, pdf_1980_Dec, col = "cornflowerblue", lwd = 2) # 1980
lines(x, pdf_1994_Dec, col = "darkseagreen4",lwd = 2) # 1994
legend("topright", legend = c("1963", "1980", "1994"),
       col = c("chocolate", "cornflowerblue", "darkseagreen4"), cex = 0.8, lty = c(1,1,1), bty = "n")
par(mfrow = c(1,1)) # reset


# look at parameters 
# 1963
param_1963_Apr <- c(Nitrate[242,5], Nitrate[242,6], Nitrate[242,7], Nitrate[242,8])
param_1963_Aug <- c(Nitrate[242,5], Nitrate[242,6], Nitrate[242,7], Nitrate[242,8])
param_1963_Dec <- c(Nitrate[242,5], Nitrate[242,6], Nitrate[242,7], Nitrate[242,8])
# 1980
param_1980_Apr <- c(Nitrate[3025,5], Nitrate[3025,6], Nitrate[3025,7], Nitrate[3025,8])
param_1980_Aug <- c(Nitrate[3025,5], Nitrate[3025,6], Nitrate[3025,7], Nitrate[3025,8])
param_1980_Dec <- c(Nitrate[6323,5], Nitrate[6323,6], Nitrate[6323,7], Nitrate[6323,8])
# 1994
param_1994_Apr <- c(Nitrate[3025,5], Nitrate[3025,6], Nitrate[3025,7], Nitrate[3025,8])
param_1994_Aug <- c(Nitrate[3025,5], Nitrate[3025,6], Nitrate[3025,7], Nitrate[3025,8])
param_1994_Dec <- c(Nitrate[6323,5], Nitrate[6323,6], Nitrate[6323,7], Nitrate[6323,8])




# create df for phytoplankton count
Phytopl <- df %>% filter(!is.na(Phytopl)) %>% 
  select(Date, Phytopl, year, month)

# create df for phytoplankton biomass
Phytopl_C <- df %>% filter(!is.na(Phytopl_C)) %>% 
  select(Date, Phytopl_C, year, month)

# create day of the year for nitrate and phytopl
Nitrate <- Nitrate %>%
  mutate(DOY  = yday(Date)) 

Phytopl <- Phytopl %>%
  mutate(DOY  = yday(Date)) 

Phytopl_C <- Phytopl_C %>%
  mutate(DOY  = yday(Date)) 


# look at nitrate patterns for different years 
years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)
years_to_plot <- c(1963, 1980, 1994)

ggplot(Nitrate %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Nitrate, color = factor(year))) +
  geom_line() +
  scale_colour_manual(values = c("Other" = "lightgrey",
                                "1963" = "chocolate",
                                 "1980" = "cornflowerblue",
                                 "1994" = "darkseagreen4")) + 
  geom_vline(xintercept = 355, linetype = "dashed") + # start of winter
  geom_vline(xintercept = 80, linetype = "dashed") + # start of spring
  geom_vline(xintercept = 173, linetype = "dashed") + # start of summer
  geom_vline(xintercept = 266, linetype = "dashed") + # start of fall
  annotate("text", x = 37, y = 145, label = "Winter", size = 4) +
  annotate("text", x = 130, y = 145, label = "Spring", size = 4) +
  annotate("text", x = 220, y = 145, label = "Summer", size = 4) +
  annotate("text", x = 315, y = 145, label = "Fall", size = 4) +
  labs(x = "Day of Year", y = "Nitrate (µmol/l)", colour = "Year") +
  theme_classic()


ggplot(Phytopl %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Phytopl, color = factor(year))) +
  geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Day of Year", y = "Phytoplankton count", colour = "Year") +
  theme_classic()

ggplot(Phytopl_C %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Phytopl_C, color = factor(year))) +
  geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Day of Year", y = "Phytoplankton biomass", colour = "Year") +
  theme_classic()



#PCA 

# create new df of interested columns and omit NA
nutrients <- df %>% 
  select(Nitrate, Nitrite, Phosphate, Silicate, Ammonium, DIN, Phytopl, Phytopl_C )
# remove NAs
nutrients <- na.omit(nutrients)

# standardize all nutrients 
std_function<-function(x) ((x-mean(x, na.rm=T)) / sd(x, na.rm=T))
std_data<-nutrients 

std_data$Nitrate <- std_function(x=std_data$Nitrate)
std_data$DIN <- std_function(x=std_data$DIN)
std_data$Nitrite <- std_function(x=std_data$Nitrite)
std_data$Silicate <- std_function(x=std_data$Silicate)
std_data$Ammonium<- std_function(x=std_data$Ammonium)
std_data$Phosphate <- std_function(x=std_data$Phosphate)
std_data$Phytopl<- std_function(x=std_data$Phytopl)
std_data$Phytopl_C <- std_function(x=std_data$Phytopl_C)

# look at correlation between nutrients
cov_mat <- round(var(std_data))
cov_mat



## SIMULATIONS

# seasonality pattern
set.seed(1234)
n <- 1825
t_index  <- seq_len(n)
mu0       <- 16
sigma0   <- 15
nu0      <- 10
tau0     <- 10

# changing rate
beta1 <- 0.6
beta2 <- 0.4

cov1 <- sin(2*pi*(t_index)/12)
cov2 <- cos(2*pi*(t_index)/12)

yt <- beta1*cov1 + beta2*cov2
plot(yt, type="l", xlab="t")



n <- 365*2             
beta1 <- 0.6; beta2 <- 0.4
period <- 365   
t_index <- 1:n

# mean changing through time
yt <- beta1 * sin(2 * pi * t_index / period) +
  beta2 * cos(2 * pi * t_index / period)


plot(t_index, yt, type="l", xlab="t (days)", ylab="yt")

# sigma and nu changing though time
sigma_t <- sigma0 + 4.146e-02
nu_t <- nu0 + -0.015507


y <- rSST(n, mu = yt, sigma = sigma_t, nu = nu0, tau = tau0)

dat <- data.frame(
  t = t_index,
  y = y,
  sd_true = sigma_t,
  mu_true = yt,
  nu_true = nu0,
  tau_true = tau0,
  dist = "SST"
)

ggplot(dat, aes(x = t)) +
  geom_line(
    aes(y = y),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)






# Simulate 

simulate_nitrate <- function(years = 33,
                             # MEAN (μ): baseline + trend + shaped seasonality
                             mu_base   = 8,            # baseline level
                             mu_trend  = 0.20,         # ↑ per year (additive)
                             peak_day  = 120,          # spring peak ~ day-of-year
                             A1        = 14,           # annual amplitude
                             A2        = -0.25,        # small 2nd harmonic to sharpen spring (keep |A2| <= ~0.33)
                             # VARIANCE (σ): log-link trend + seasonal bump in winter
                             sigma0    = 3.5,          # starting SD
                             sigma_trend_mult = 3.0,   # by the last year SD ≈ sigma0 * sigma_trend_mult
                             sigma_seas = 0.35,        # wintery ↑ in spread (0..~0.6)
                             winter_shift = 0.15,      # shifts “winter spread” timing
                             # Distribution (Gaussian or SST)
                             family = c("GAUS", "SST"),
                             nu = 0.9, tau = 7,        # SST shape params (ignored if GAUS)
                             seed = 42) {
  
  set.seed(seed)
  family <- match.arg(family)
  
  # timeline
  T  <- 365 * years
  t  <- 1:T
  u  <- (t - 1) %% 365 / 365                 # within-year position [0,1)
  yr <- floor((t - 1) / 365)                 # 0,1,...,years-1
  
  # ---- MEAN: trend + single spring peak ----
  phi <- 2*pi*peak_day/365                   # place the peak in spring
  # annual wave + small 2nd harmonic to skew shape (kept small so only 1 peak)
  seas <- sin(2*pi*u - phi) + A2 * sin(4*pi*u - 2*phi)
  mu_t <- mu_base + mu_trend * (yr) + A1 * seas
  
  # ---- SIGMA: multiplicative growth + wintery spread ----
  # log-sigma rises linearly from log(sigma0) to log(sigma0*sigma_trend_mult)
  log_sigma_t <- log(sigma0) +
    log(sigma_trend_mult) * (yr / (years - 1)) +
    sigma_seas * cos(2*pi*(u + winter_shift))  # bigger spread in winter
  sigma_t <- pmax(0.1, exp(log_sigma_t))
  
  # ---- draw series ----
  if (family == "GAUS") {
    y <- rnorm(T, mean = mu_t, sd = sigma_t)
  } else {
    # skew + heavy tails like nitrate extremes
    if (!requireNamespace("gamlss.dist", quietly = TRUE))
      stop("Please install.packages('gamlss.dist') for SST draws.")
    y <- gamlss.dist::rSST(T, mu = mu_t, sigma = sigma_t, nu = nu, tau = tau)
  }
  
  # keep nonnegative (nitrate units)
  y <- pmax(0, y)
  
  data.frame(t = t, year = yr, u = u, mu_t = mu_t, sigma_t = sigma_t, y = y)
}

start_date <- as.Date("1962-01-01")


# add Date, Year, Month (no need to change how you simulated y)
sim <- sim %>%
  mutate(
    Date  = start_date + (t - 1),
    Year  = year(Date),
    Month = factor(month(Date), levels = 1:12, labels = month.abb)
  )

# plot that sucka
ggplot(sim, aes(Date, y)) +
  geom_line(color = "#7A0C0C") +
  labs(x = "Date", y = "Value") +
  theme_classic()

# center ghe year and date 
sim$Year_C <- sim$Year - min(sim$Year)
sim$Date_C <- sim$Date - min(sim$Date)

#===========================================================================================
# seasonality for just the mean
m1 <- gamlss(y ~ Year + Month, family = SST(), data = sim,
             method = mixed(10,200),
             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(m1)

# seasonality for other param
m2 <- gamlss(y ~ Year_C + Month, sigma.fo = ~ Year_C + Month, nu.fo = ~ Date, family = SST(), data = sim,
             mu.start = mean(sim$mu_t), sigma.start = sd(sim$sigma_t),
             method = mixed(10,200),
             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(m2)

# mean only model 
mean_m <- gamlss(y ~ Date_C, family = SST(), data = sim,
             mu.start = mean(sim$y), sigma.start = sd(sim$y),
             method = mixed(10,200),
             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


# all param (time)
param_m <- gamlss(y ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date,  family = SST(), data = sim,
                  mu.start = mean(sim$mu_t), sigma.start = sd(sim$sigma_t),
                 method = mixed(10,200),
                 control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
#===========================================================================================


# predict mu based on seasonality model
sim$mu_hat <- fitted(m2, what = "mu")
sim$sigma_hat <- fitted(m2, what = "sigma")
sim$nu_hat <- fitted(m2, what = "nu")
sim$tau_hat <- fitted(m2, what = "tau")

# plot best of fit line to data
ggplot(sim, aes(Date, y)) +
  geom_line(color = "grey") +
  geom_line(aes(y = mu_hat), color = "darkolivegreen", linewidth = 1) +
  geom_hline(yintercept = mean(sim$y), linetype = "dashed", color = "black") +
  labs(x = "Date", y = "Value") +
  theme_classic()


