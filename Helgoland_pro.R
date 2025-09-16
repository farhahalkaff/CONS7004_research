# Set the directory where your files are located
data_directory <- "dataset"

# Get a list of all files with a specific extension (e.g., .csv)
files <- list.files(path = data_directory, pattern = "\\.csv$", full.names = TRUE)
files <- files[!grepl("Helgoland_2002.csv", files)]
files <- files[!grepl("Helgoland_1998.csv", files)]

# Use lapply to read each file into a list of data frames
# Replace read.csv with the appropriate function for your file type (e.g., read.delim, read_excel from 'readxl' package)
list_of_dataframes <- lapply(files, read.csv)

df <- files %>%
  lapply(read.csv) %>%
  bind_rows()

library(dplyr)
library(ggplot2)
library(patchwork)
library(lubridate)
library(tidyr)

library(mvgam)           # Fit, interrogate and forecast DGAMs
library(forecast)        # Construct fourier terms for time series
library(gratia)          # Graceful plotting of smooth terms
library(marginaleffects) # Interrogating regression models
library(janitor)         # Creating clean, tidy variable names

# add year as a column
df$year <- year(ymd_hm(df$Date))
df <- df %>%
  mutate(Date = as.Date(ymd_hm(Date)))



# =========================
# changes in temperature through the years
# =========================

## Calculate yearly average temperature
yearly_avg <- df %>%
  group_by(year) %>%
  summarise(avg_temp = mean(Temp, na.rm = TRUE))
### plot the average temp for each year
ggplot(yearly_avg, aes(x = year, y = avg_temp)) +
  geom_line() +
  geom_point() +
  labs(title = "Average Temperature per Year",
       x = "Year",
       y = "Average Temperature")



# =========================
# How higher-moments in resources changes over time
# =========================

# see how nitrite differ over the years 
library(dplyr)
library(e1071)  # for skewness and kurtosis
## variance
yearly_stats <- df %>%
  group_by(year) %>%
  summarise(
    mean_nitrite = mean(Nitrite, na.rm = TRUE),
    var_nitrite = var(Nitrite, na.rm = TRUE),
    skew_nitrite = skewness(Nitrite, na.rm = TRUE),
    kurt_nitrite = kurtosis(Nitrite, na.rm = TRUE),
    mean_nitrate = mean(Nitrate, na.rm = TRUE),
    var_nitrate = var(Nitrate, na.rm = TRUE),
    skew_nitrate = skewness(Nitrate, na.rm = TRUE),
    kurt_nitrate = kurtosis(Nitrate, na.rm = TRUE),
    mean_phospate = mean(Phosphate, na.rm = TRUE),
    var_phosphate = var(Phosphate, na.rm = TRUE),
    skew_phosphate = skewness(Phosphate, na.rm = TRUE),
    kurt_phosphate = kurtosis(Phosphate, na.rm = TRUE)
  )
##visualize for variance 
ni_v <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = var_nitrite), color = 'blue') +
  labs(title = "Variance of Nitrite Over Years", y = "Variance")
##visualize for skewness 
ni_s <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = skew_nitrite), color = 'blue') +
  labs(title = "Skewness of Nitrite Over Years", y = "Skewness")
##visualize for kurtosis 
ni_k <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = kurt_nitrite), color = 'blue') +
  labs(title = "Kurtosis of Nitrite Over Years", y = "Kurtosis")
###stack them together 
(ni_v / ni_s / ni_k) + plot_layout(ncol = 1)


# see how nitrate differ over the years 
##visualize for variance 
na_v <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = var_nitrate), color = 'goldenrod') +
  labs(title = "Variance of Nitrate Over Years", y = "Variance")
##visualize for skewness 
na_s <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = skew_nitrate), color = 'goldenrod') +
  labs(title = "Skewness of Nitrate Over Years", y = "Skewness")
##visualize for kurtosis 
na_k <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = kurt_nitrate), color = 'goldenrod') +
  labs(title = "Kurtosis of Nitrate Over Years", y = "Kurtosis")
###stack them together 
(na_v / na_s / na_k) + plot_layout(ncol = 1)


# see how phosphate differ over the years 
##visualize for variance 
p_v <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = var_phosphate), color = 'darkolivegreen') +
  labs(title = "Variance of Phosphate Over Years", y = "Variance")
##visualize for skewness 
p_s <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = skew_phosphate), color = 'darkolivegreen') +
  labs(title = "Skewness of Phosphate Over Years", y = "Skewness")
##visualize for kurtosis 
p_k <- ggplot(yearly_stats, aes(x = year)) +
  geom_line(aes(y = kurt_phosphate), color = 'darkolivegreen') +
  labs(title = "Kurtosis of Phosphate Over Years", y = "Kurtosis")
###stack them together 
(p_v / p_s / p_k) + plot_layout(ncol = 1)




# =========================
# compare nutrients and phytoplankton
# =========================

nutrients <- c("Nitrate", "Nitrite", "Silicate", "Phosphate", "DIN", "Silicate", "Ammonium") 
nitrate <- "Nitrate"
phyto_col <- "Phytopl"
phyto_biom <- "Phytopl_C"

df_year <- df %>%
  group_by(year) %>%
  summarise(across(all_of(c(nutrients, phyto_col)), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup()

ggplot(df_year, aes(x = year, y = z, color = variable)) +
  # geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_smooth(se = FALSE, method = "loess", span = 0.3, linewidth = 0.6) +
  labs(x = "Year", y = "Standardized (z-score)",
       color = "Variable") +
  theme_minimal()


##################################
### Looking at trends per year ###
##################################

target_year <- 1968

df_one <- df %>%
  mutate(year = year(Date), month = month(Date, label = TRUE)) %>%
  filter(year == target_year) %>%
  group_by(month) %>%
  summarise(across(all_of(c(nitrate, phyto_col)), ~mean(., na.rm = TRUE)),
            .groups = "drop")

df_long <- df_one %>%
  pivot_longer(-month, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup()

ggplot(df_long, aes(month, z, color = variable, group = variable)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  labs(x = "Month", y = "Z-score", color = "Variable") +
  theme_minimal()



# =========================
# Raw plots
# =========================

####################################################
### Nitrate and fit models of different families ###
####################################################

# create own df to remove NAs 
Nitrate <- df %>% filter(!is.na(Nitrate)) %>% 
  select(Date, Nitrate, year)

# plot that sucka (intercept model)
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "grey", na.rm = TRUE) +
  geom_abline(intercept = niSSTm2$mu.coefficients[1], slope = niSSTm2$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # only mean as function of time 
  geom_abline(intercept = niSSTm3$mu.coefficients[1], slope = niSSTm3$mu.coefficients[2], 
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu as function of time
  geom_hline(yintercept = 16.718, linetype = "dashed", color = "black") + # intercept mean
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()


# plot for a specific year 
ggplot(dplyr::filter(Nitrate, year(Date) == 1994),
       aes(Date, Nitrate)) +
  geom_line(color = "black")

stl(Nitrate$Nitrate, s.window = 9) %>%
  autoplot()

head(Nitrate, 24)



### CHECK DISTRIBUTION ###

# histogram
hist(Nitrate$Nitrate)

# density plot 
ggplot(Nitrate, aes(y = Nitrate)) +
  geom_density(
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()

# PDF
pdf <- dSST(x, mu = 16.718, sigma = exp(2.75160), 
            nu = exp(2.2996), tau = exp(0.9038) + 2)
x <- seq(min(niSSTm$y, na.rm = TRUE),
         max(niSSTm$y, na.rm = TRUE),
         length.out = 200)
plot(x, pdf, type = "l")




### DIFFERENT FAMILY MODELS ###

#============================================================================================
# Intercept model 
niSSTm <- gamlss(Nitrate ~ 1, family = SST(), data = Nitrate,
                 mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                 method = mixed(10,200),
                 control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm)

# Only mean changing through time 
niSSTm2 <- gamlss(Nitrate ~ Date, family = SST(), data = Nitrate,
                 #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                 method = mixed(10,200),
                 control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm2)

# mean, sigma and nu changing through time 
niSSTm3 <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrate,
                  #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                  method = mixed(10,200),
                  control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm3)
#============================================================================================


niTF2m <- gamlss(Nitrate ~ 1, family = TF2(), data = Nitrate, 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niSN1m <- gamlss(Nitrate ~ 1, family = SN1(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niNOm <- gamlss(Nitrate ~ 1, family = NO(), data = Nitrate, 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

summary(Nitrate$Nitrate)
sum(is.nan(Nitrate$Nitrate))

summary(niSSTm)
summary(niTF2m)
AIC(niTF2m)
AIC(niTF2mtest)


### SST FAMILY MODELS WITH TIME-VARYING PARAMETERS ###

niSSTm_1 <- gamlss(Nitrate ~ Date, family = SST(), data = Nitrate, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_1)
with(Nitrate, plot(Nitrate ~ Date))


curve(cbind(1,x)%*%coef(niSSTm_1), add =T, col = "red", lwd=2)


niSSTm_sig <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, family = SST(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_sig)

niSSTm_nu <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_nu)

niSSTm_tau <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), 
                     data = Nitrate, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_tau)

#omit <- na.omit(df)


# summary
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  niSSTm_1 = c(coefAll(niSSTm_1)$mu[1],coefAll(niSSTm_1)$mu[2] ,exp(coefAll(niSSTm_1)$sigma),"NA", 
                exp(coefAll(niSSTm_1)$nu),"NA", (exp(coefAll(niSSTm_1)$tau) + 2), "NA", AIC(niSSTm_1)),
  niSSTm_sig = c(coefAll(niSSTm_sig)$mu[1],coefAll(niSSTm_sig)$mu[2] ,exp(coefAll(niSSTm_sig)$sigma[1]),
                 exp(coefAll(niSSTm_sig)$sigma[2]), exp(coefAll(niSSTm_sig)$nu),"NA", 
                 (exp(coefAll(niSSTm_sig)$tau) + 2), "NA", AIC(niSSTm_sig)),
  niSSTm_nu = c(coefAll(niSSTm_nu)$mu[1],coefAll(niSSTm_nu)$mu[2], exp(coefAll(niSSTm_nu)$sigma[1]),
                 exp(coefAll(niSSTm_nu)$sigma[2]), exp(coefAll(niSSTm_nu)$nu[1]),exp(coefAll(niSSTm_nu)$nu[2]), 
                 (exp(coefAll(niSSTm_nu)$tau) + 2), "NA", AIC(niSSTm_nu))
)
print(param_summary)


# V = 1 --> symmetrical 
# V > 1 --> positively skewed 
# 0 < V < 1 --> negatively skewed 
# Decreasing T increases the heaviness of the tail
# When V = 1, = TF2 distribution 
# As V reaches infinity, it tends to be half t distribution TF
# As T reaches infinity, it tends to Skew normal type 2, SN2
# If both V and T reaches infinity, it tends to a half normal distribution


#############################################
############################################# LOG

niSSTm_1b <- gamlss(log(Nitrate) ~ Date, family = SST(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_1b)
with(Nitrate, plot(Nitrate ~ Date))
curve((cbind(1,x)%*%coef(niSSTm_1b)), add =T, col = "red", lwd=2)


niSSTm_sigb <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, family = SST(), data = Nitrate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_sigb)

niSSTm_nub <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrate, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_nub)

niSSTm_taub <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), 
                     data = Nitrate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_taub)


# logging takes away heavy tails #
# pdf changing over time 
# ecology 
# mean follow season vs not 
# mean centering

# full intercept model for all nutrients 
# fit a line (monotonic for now - linear) -- also non-monotonic to see what the data is doing 
# intercept - linear - quadratic - spline
# stick with one fam distribution

# smaller the bar the more important it is
# consequences of assumptions 

param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  niSSTm_1b = c(coefAll(niSSTm_1b)$mu[1],coefAll(niSSTm_1b)$mu[2] ,exp(coefAll(niSSTm_1b)$sigma),"NA", 
               exp(coefAll(niSSTm_1b)$nu),"NA", (exp(coefAll(niSSTm_1b)$tau) + 2), "NA", AIC(niSSTm_1b)),
  niSSTm_sigb = c(coefAll(niSSTm_sigb)$mu[1],coefAll(niSSTm_sigb)$mu[2] ,exp(coefAll(niSSTm_sigb)$sigma[1]),
                 exp(coefAll(niSSTm_sigb)$sigma[2]), exp(coefAll(niSSTm_sigb)$nu),"NA", 
                 (exp(coefAll(niSSTm_sigb)$tau) + 2), "NA", AIC(niSSTm_sigb)),
  niSSTm_nub = c(coefAll(niSSTm_nub)$mu[1],coefAll(niSSTm_nub)$mu[2], exp(coefAll(niSSTm_nub)$sigma[1]),
                exp(coefAll(niSSTm_nub)$sigma[2]), exp(coefAll(niSSTm_nub)$nu[1]),exp(coefAll(niSSTm_nub)$nu[2]), 
                (exp(coefAll(niSSTm_nub)$tau) + 2), "NA", AIC(niSSTm_nub))
)
print(param_summary)

#############

# scale by df and center the mean
Nitrate$t_index <- seq(0, nrow(Nitrate) - 1)

Nitrate <- Nitrate %>%
  mutate(
    # keep Date as is
    Date = as.Date(Date),              
    t_num = as.numeric(Date),          
    t_scaled = (t_num - min(t_num, na.rm = TRUE)) /
      (max(t_num, na.rm = TRUE) - min(t_num, na.rm = TRUE)),
    t_centered = t_index - mean(t_index)
  )

max(Nitrate$t_centered)
min(Nitrate$t_centered)

# center the mean
Nitrate <- Nitrate %>%
  mutate(t_centered = t_index - mean(t_index))


# new model with t_scaled
niSSTm_nu <- gamlss(Nitrate ~ t_scaled, sigma.fo = ~ t_scaled, nu.fo = ~ t_scaled, family = SST(), data = Nitrate, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_nu)

# keep model unscaled but use t_index
niSSTm_nu <- gamlss(Nitrate ~ t_index, sigma.fo = ~ t_index, nu.fo = ~ t_index, family = SST(), data = Nitrate, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_nu)


# function for pdf
pdf_SST_at_date <- function(niSSTm_nu, Nitrate, date_str, x = NULL, n = 200,
                            time_var = "t_index", date_var = "Date") {
  # 1. Find the row in df matching the date of interest
  t_val <- Nitrate[[time_var]][as.Date(Nitrate[[date_var]]) == as.Date(date_str)]
  if (length(t_val) == 0) stop("Date not found in data frame.")
  
  # 2. Build newdata with the correct covariate
  newdat <- data.frame(setNames(list(t_val), time_var))
  
  # 3. Predict distribution parameters on natural scale
  mu    <- predict(niSSTm_nu, "mu",    newdata = newdat, type = "response")
  sigma <- predict(niSSTm_nu, "sigma", newdata = newdat, type = "response")
  nu    <- predict(niSSTm_nu, "nu",    newdata = newdat, type = "response")
  tau   <- predict(niSSTm_nu, "tau",   newdata = newdat, type = "response")
  
  # 4. Create x grid if none supplied
  if (is.null(x)) {
    x <- seq(min(niSSTm_nu$y, na.rm = TRUE),
             max(niSSTm_nu$y, na.rm = TRUE),
             length.out = n)
  }
  
  # 5. Evaluate PDF
  dens <- dSST(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
  
  list(x = x, density = dens,
       params = c(mu = mu, sigma = sigma, nu = nu, tau = tau))
}

res <- pdf_SST_at_date(niSSTm_nu, Nitrate, "1962-01-15",
                       time_var = "t_index", date_var = "Date")

plot(res$x, res$density, type = "l")
res$params



# manually
x <- seq(min(niSSTm_nu$y, na.rm = TRUE),
         max(niSSTm_nu$y, na.rm = TRUE),
         length.out = 200)
pdf <- dSST(x, mu = (3022*1.59e+01)+2.80e-03, sigma = exp((3022*2.47e+00) + 1.64e-04), 
            nu = exp((3022*1.83e+00) + 4.61e-05), tau = exp(3022*2.2234) + 2)
plot(x, pdf, type = "l")

# mu 
3022*1.59e+01+2.80e-03

# sigma exp after
exp((4971*1.93e+00) + 1.64e-04)
# sigma exp before 
4971*exp(1.93e+00) + exp(1.64e-04)

# nu exp after 
exp((0.765*1.679) + 0.299)
# nu exp before 
0.765*exp(1.679) + exp(0.299)

# tau exp after 
exp(0.765*2.2234) + 2
# tau exp before 
(exp(2.2234) + 2)*0.765



#############################################
############################################# SHASHo

niSSTm_1c <- gamlss(log(Nitrate) ~ Date, family = SHASHo(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_1c)


niSSTm_sigc <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, family = SHASHo(), data = Nitrate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_sigc)

niSSTm_nuc <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SHASHo(), data = Nitrate, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_nu)

niSSTm_tauc <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = BCT(), 
                     data = Nitrate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_tau)


#####================
# niSSTm_1c <- gamlss(Nitrate ~ Date, family = BCT(), data = na.omit(df), 
#                     method = mixed(5, 200),
#                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_1c)
# with(Nitrate, plot(Nitrate ~ Date))
# curve(cbind(1,x)%*%coef(niSSTm_1), add =T, col = "red", lwd=2)
# 
# 
# niSSTm_sigc <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, family = BCT(), data = na.omit(df), 
#                       method = mixed(5, 200),
#                       control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_sig)
# 
# niSSTm_nuc <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = BCT(), data = na.omit(df), 
#                      method = mixed(5, 200),
#                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_nu)
# 
# niSSTm_tauc <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = BCT(), 
#                       data = na.omit(df), 
#                       method = mixed(5, 200),
#                       control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_tau)
#######==================

# summary
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  niSSTm_1c = c(coefAll(niSSTm_1c)$mu[1],coefAll(niSSTm_1c)$mu[2] ,exp(coefAll(niSSTm_1c)$sigma),"NA", 
               coefAll(niSSTm_1c)$nu,"NA", (exp(coefAll(niSSTm_1c)$tau)), "NA", AIC(niSSTm_1c)),
  niSSTm_sigc = c(coefAll(niSSTm_sigc)$mu[1],coefAll(niSSTm_sigc)$mu[2] ,exp(coefAll(niSSTm_sigc)$sigma[1]),
                 exp(coefAll(niSSTm_sigc)$sigma[2]), coefAll(niSSTm_sigc)$nu,"NA", 
                 (exp(coefAll(niSSTm_sigc)$tau)), "NA", AIC(niSSTm_sigc)),
  niSSTm_nuc = c(coefAll(niSSTm_nuc)$mu[1],coefAll(niSSTm_nuc)$mu[2], exp(coefAll(niSSTm_nuc)$sigma[1]),
                exp(coefAll(niSSTm_nuc)$sigma[2]), coefAll(niSSTm_nuc)$nu[1], coefAll(niSSTm_nuc)$nu[2], 
                (exp(coefAll(niSSTm_nuc)$tau)), "NA", AIC(niSSTm_nuc))
)
print(param_summary)

# V = 0 --> symmetrical 
# V > 0 --> positively skewed 
# V < 0 --> negatively skewed 
# T < 1 --> heavier tails than normal
# T > 1 --> lighter tails than normal
# centile skewness and kurtosis


#############################################
############################################# Skew exponential 

niSEPm_1 <- gamlss(log(Nitrate) ~ Date, family = SEP3(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSEPm_1)


niSEPm_sig <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, family = SEP3(), data = Nitrate, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niSEPm_nu <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SEP3(), data = Nitrate, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSEPm_nu)

# summary 
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  niSEPm_1 = c(coefAll(niSEPm_1)$mu[1],coefAll(niSEPm_1)$mu[2] ,exp(coefAll(niSEPm_1)$sigma),"NA", 
                exp(coefAll(niSEPm_1)$nu),"NA", (exp(coefAll(niSEPm_1)$tau)), "NA", AIC(niSEPm_1)),
  niSEPm_sig = c(coefAll(niSEPm_sig)$mu[1],coefAll(niSEPm_sig)$mu[2] ,exp(coefAll(niSEPm_sig)$sigma[1]),
                  exp(coefAll(niSEPm_sig)$sigma[2]), exp(coefAll(niSEPm_sig)$nu),"NA", 
                  (exp(coefAll(niSEPm_sig)$tau)), "NA", AIC(niSEPm_sig)),
  niSEPm_nu = c(coefAll(niSEPm_nu)$mu[1],coefAll(niSEPm_nu)$mu[2], exp(coefAll(niSEPm_nu)$sigma[1]),
                 exp(coefAll(niSEPm_nu)$sigma[2]), exp(coefAll(niSEPm_nu)$nu[1]),exp(coefAll(niSEPm_nu)$nu[2]), 
                 (exp(coefAll(niSEPm_nu)$tau)), "NA", AIC(niSEPm_nu))
)
print(param_summary)

# V = 0 --> symmetrical 
# V > 1 --> positively skewed 
# 0 < V < 1 --> negatively skewed 
# T < 2 --> moment leptokurtic (heavy tails)
# T > 2 --> moment platykurtic (lighter tials) or leptokurtic




#############################################
############################################# SST left-truncated

library(gamlss.tr)
gen.trun(0,"SHASHo",type="left")

niSSTm_1tr <- gamlss(Nitrate ~ Date, family = SSTtr(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 900, c.crit = 0.01, trace = FALSE))


niSHASHom <- gamlss(Nitrate + 100 ~ 1, family = SHASHotr(), data = Nitrate,
                 method = mixed(10,200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSHASHom)


#======================
# SEASONALITY
#======================
 
# See missing values for each month of each year
Nitrate_NA <- Nitrate %>%
  mutate(Year = year(Date),
         Month = factor(month(Date), levels = 1:12, labels = month.abb)) %>%
  group_by(Year, Month) %>%
  summarize(Nitrate = mean(Nitrate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month, values_from = Nitrate) %>%
  arrange(Year)

# Average values of nitrate per months each year 
Nitrate_monthly <- Nitrate %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(Nitrate = mean(Nitrate), .groups = "drop") %>%
  arrange(Year, Month)

# Add missing data (May and June of 1965, Jan of 1970)
Nitrate_monthly <- Nitrate_monthly %>% 
  add_row(Year = 1965, Month = 5, Nitrate = 25.82, .after = 40) %>% 
  add_row(Year = 1965, Month = 6, Nitrate = 19.265, .after = 41) %>% 
  add_row(Year = 1970, Month = 1, Nitrate = 8.29, .after = 96)

# turn dataset into vector 
Nitrate_monthly_v <- ts(
  Nitrate_monthly$Nitrate,
  start = c(min(Nitrate_monthly$Year), min(Nitrate_monthly$Month)),
  frequency = 12
)

# check that shit out 
head(Nitrate_monthly_v, 72)

# STL decomposition (Seasonal and Trend decomposition using Loess)
stl(Nitrate_monthly_v, s.window = 12) %>%  #== LOOK INTO S.WINDOW ==#
  autoplot()
?stl

# split data into training and testing set 
NM <- series_to_mvgam(Nitrate_monthly_v, freq = frequency(Nitrate_monthly_v))
dplyr::glimpse(NM$data_train)

# plot them 
plot_mvgam_series(data = NM$data_train,
                  newdata = NM$data_test)




######################################################
### Phosphate and fit models of different families ###
######################################################

# create own df to remove NAs
Phosphate <- df %>% filter(!is.na(Phosphate)) %>% 
  select(Date, Phosphate, year)

# plot that sucka
ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "grey") +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0.82, linetype = "dashed", color = "black") + 
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()

# plot a for a specific year 
ggplot(dplyr::filter(Phosphate, year(Date) == 1970),
       aes(Date, Phosphate)) +
  geom_line(color = "black")

# rescaled phosphate
Phosphate$y_rescaled <- Phosphate$Phosphate * 100

### CHECK DISTRIBUTION ###

# histogram 
hist(Phosphate$Phosphate)

# density plot 
ggplot(Phosphate, aes(y = Phosphate)) +
  geom_density(
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### MODELS ###

# Intercept model (different family tho ;( ) => mu is mode
phSEPm1 <- gamlss(Phosphate ~ 1, family = SEP3(), data = Phosphate,
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phSEPm1)
mfv(Phosphate$Phosphate)


phSSTm <- gamlss(Phosphate ~ poly(Date,2), family = SST(), data = Phosphate, 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phTF2m <- gamlss(Phosphate ~ Date, family = TF2(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phSN1m <- gamlss(Phosphate ~ Date, family = SN1(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phNOm <- gamlss(Phosphate ~ Date, family = NO(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))



AIC(phTF2m)
AIC(phSN1m)
AIC(phNOm)


#======================
# SEASONALITY
#======================

# Check missing data for each month  
Phosphate_NA <- Phosphate %>%
  mutate(Year = year(Date),
         Month = factor(month(Date), levels = 1:12, labels = month.abb)) %>%
  group_by(Year, Month) %>%
  summarize(Phosphate = mean(Phosphate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month, values_from = Phosphate) %>%
  arrange(Year)

# Average values of nitrate per months each year 
Phosphate_monthly <- Phosphate %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(Phosphate = mean(Phosphate, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year, Month)

# Add missing data (May and June of 1965)
Phosphate_monthly <- Phosphate_monthly %>% 
  add_row(Year = 1965, Month = 5, Phosphate = 0.655, .after = 40) %>% 
  add_row(Year = 1965, Month = 6, Phosphate = 0.228, .after = 41) 

# turn dataset into vector 
Phosphate_monthly_v <- ts(
  Phosphate_monthly$Phosphate,
  start = c(min(Phosphate_monthly$Year), min(Phosphate_monthly$Month)),
  frequency = 12
)

# check that shit out 
head(Phosphate_monthly_v, 72)

# STL decomposition (Seasonal and Trend decomposition using Loess)
stl(Phosphate_monthly_v, s.window = 7) %>%  #== LOOK INTO S.WINDOE ==#
  autoplot()
?stl

# split data into training and testing set 
PM <- series_to_mvgam(Phosphate_monthly_v, freq = frequency(Phosphate_monthly_v))
dplyr::glimpse(PM$data_train)

# plot them 
plot_mvgam_series(data = PM$data_train,
                  newdata = PM$data_test)




####################################################
### Nitrite and fit models of different families ###
####################################################

# Create own df to remove NAs
Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite, year)

# plot that sucka
ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "grey") + 
  geom_abline(intercept = niiSSTm$mu.coefficients[1], slope = niiSSTm$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) +
  geom_hline(yintercept = 0.81659, linetype = "dashed", color = "black") + 
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()

# plot a for a specific year 
ggplot(dplyr::filter(Nitrite, year(Date) == 1970),
       aes(Date, Nitrite)) +
  geom_line(color = "black")


### CHECK DISTRIBUTION ###

# histogram
hist(Nitrite$Nitrite)

# density plot 
ggplot(Nitrite, aes(y = Nitrite)) +
  geom_density(
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### MODELS ###

# intercept model
niiSSTm1 <- gamlss(Nitrite ~ 1, family = SST(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niiSSTm)


niiSSTm <- gamlss(Nitrite ~ Date, family = SST(), data = Nitrite, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niiTF2m <- gamlss(Nitrite ~ Date, family = TF2(), data = Nitrite, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niiSN1m <- gamlss(Nitrite ~ Date, family = SN1(), data = Nitrite, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niiNOm <- gamlss(Nitrite ~ Date, family = NO(), data = Nitrite, 
                method = mixed(5, 200),
                control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


AIC(niiSSTm)
AIC(niiTF2m)
AIC(niiSN1m)
AIC(niiNOm)


#### CHANGING PARAMTER AS FUNCTION OF TIME WITH SST MODEL ####

niiSSTm_1 <- gamlss(Nitrite ~ Date, family = SST(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niiSSTm_1poly <- gamlss(Nitrite ~ poly(Date,3), family = SST(), data = Nitrite, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

summary(niiSSTm_1poly)

niiSSTm_sig <- gamlss(Nitrite ~ Date, sigma.fo = ~ Date, family = SST(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niiSSTm_nu <- gamlss(Nitrite ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrite, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niiSSTm_tau <- gamlss(Nitrite ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), data = Nitrite, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


# summary
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  niiSSTm_1 = c(coefAll(niiSSTm_1)$mu[1],coefAll(niiSSTm_1)$mu[2] ,exp(coefAll(niiSSTm_1)$sigma),"NA", 
               exp(coefAll(niiSSTm_1)$nu),"NA", (exp(coefAll(niiSSTm_1)$tau) + 2), "NA", AIC(niiSSTm_1)),
  niiSSTm_sig = c(coefAll(niiSSTm_sig)$mu[1],coefAll(niiSSTm_sig)$mu[2] ,exp(coefAll(niiSSTm_sig)$sigma[1]),
                 exp(coefAll(niiSSTm_sig)$sigma[2]), exp(coefAll(niiSSTm_sig)$nu),"NA", 
                 (exp(coefAll(niiSSTm_sig)$tau) + 2), "NA", AIC(niiSSTm_sig)),
  niiSSTm_nu = c(coefAll(niiSSTm_nu)$mu[1],coefAll(niiSSTm_nu)$mu[2], exp(coefAll(niiSSTm_nu)$sigma[1]),
                exp(coefAll(niiSSTm_nu)$sigma[2]), exp(coefAll(niiSSTm_nu)$nu[1]),exp(coefAll(niiSSTm_nu)$nu[2]), 
                (exp(coefAll(niiSSTm_nu)$tau) + 2), "NA", AIC(niiSSTm_nu)),
  niiSSTm_tau = c(coefAll(niiSSTm_tau)$mu[1],coefAll(niiSSTm_tau)$mu[2], exp(coefAll(niiSSTm_tau)$sigma[1]),
                 exp(coefAll(niiSSTm_tau)$sigma[2]), exp(coefAll(niiSSTm_tau)$nu[1]),exp(coefAll(niiSSTm_tau)$nu[2]), 
                 (exp(coefAll(niiSSTm_tau)$tau[1]) + 2), (exp(coefAll(niiSSTm_tau)$tau[2]) + 2), AIC(niiSSTm_tau))
)
print(param_summary)



#======================
# SEASONALITY
#======================

# Check missing data for each month  
Nitrite_NA <- Nitrite %>%
  mutate(Year = year(Date),
         Month = factor(month(Date), levels = 1:12, labels = month.abb)) %>%
  group_by(Year, Month) %>%
  summarize(Nitrite = mean(Nitrite, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month, values_from = Nitrite) %>%
  arrange(Year)

# Average values of nitrate per months each year 
Nitrite_monthly <- Nitrite %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(Nitrite = mean(Nitrite, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year, Month)

# Add missing data (May and June of 1965)
Nitrite_monthly <- Nitrite_monthly %>% 
  add_row(Year = 1965, Month = 5, Nitrite = 1.23, .after = 40) %>% 
  add_row(Year = 1965, Month = 6, Nitrite = 0.752, .after = 41) 

# turn dataset into vector 
Nitrite_monthly_v <- ts(
  Nitrite_monthly$Nitrite,
  start = c(min(Nitrite_monthly$Year), min(Nitrite_monthly$Month)),
  frequency = 12
)

# check that shit out 
head(Nitrite_monthly_v, 72)

# STL decomposition (Seasonal and Trend decomposition using Loess)
stl(Nitrite_monthly_v, s.window = 7) %>%  #== LOOK INTO S.WINDOW ==#
  autoplot()
?stl

# split data into training and testing set 
NaM <- series_to_mvgam(Nitrite_monthly_v, freq = frequency(Nitrite_monthly_v))
dplyr::glimpse(NaM$data_train)

# plot them 
plot_mvgam_series(data = NaM$data_train,
                  newdata = NaM$data_test)


################################################
### DIN and fit models of different families ###
################################################

# create it's own df to remove NAs
DIN <- df %>% filter(!is.na(DIN)) %>% 
  select(Date, DIN, year)

# plot that sucka
ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "grey") +
  geom_abline(intercept = DINSSTm$mu.coefficients[1], slope = DINSSTm$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # only mean as function of time 
  geom_abline(intercept = DINSSTm_nu$mu.coefficients[1], slope = DINSSTm_nu$mu.coefficients[2], 
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu as function of time
  geom_hline(yintercept = 23.182, linetype = "dashed", color = "black") + # mean intercept
  labs(x = "Time", y = "DIN (µmol/l)") +
  theme_minimal()

# plot a for a specific year 
ggplot(dplyr::filter(DIN, year(Date) == 1970),
       aes(Date, DIN)) +
  geom_line(color = "black")


### CHECK DISTRIBUTION ###

# histogram 
hist(DIN$DIN)

# density plot 
ggplot(DIN, aes(y = DIN)) +
  geom_density(
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### MODELS ###

# intercept model 
DINSSTm1 <- gamlss(DIN ~ 1, family = SST(), data = DIN,
                  method = mixed(10, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(DINSSTm1)


DINSSTm <- gamlss(DIN ~ Date, family = SST(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

DINTF2m <- gamlss(DIN ~ Date, family = TF2(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

DINSN1m <- gamlss(DIN ~ Date, family = SN1(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

DINNOm <- gamlss(DIN ~ Date, family = NO(), data = DIN, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


AIC(DINSSTm)
AIC(DINTF2m)
AIC(DINSN1m)
AIC(DINNOm)


#### CHANGING PARAMTER AS FUNCTION OF TIME WITH SST MODEL ####

DINSSTm_1 <- gamlss(log(DIN) ~ Date, family = SST(), data = DIN, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


DINSSTm_sig <- gamlss(log(DIN) ~ Date, sigma.fo = ~ Date, family = SST(), data = DIN, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


DINSSTm_nu <- gamlss(DIN ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = DIN, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


DINSSTm_tau <- gamlss(log(DIN) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), data = DIN, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


# summary
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  DINSSTm_1 = c(coefAll(DINSSTm_1)$mu[1],coefAll(DINSSTm_1)$mu[2] ,exp(coefAll(DINSSTm_1)$sigma),"NA", 
                exp(coefAll(DINSSTm_1)$nu),"NA", (exp(coefAll(DINSSTm_1)$tau) + 2), "NA", AIC(DINSSTm_1)),
  DINSSTm_sig = c(coefAll(DINSSTm_sig)$mu[1],coefAll(DINSSTm_sig)$mu[2] ,exp(coefAll(DINSSTm_sig)$sigma[1]),
                  exp(coefAll(DINSSTm_sig)$sigma[2]), exp(coefAll(DINSSTm_sig)$nu),"NA", 
                  (exp(coefAll(DINSSTm_sig)$tau) + 2), "NA", AIC(DINSSTm_sig)),
  DINSSTm_nu = c(coefAll(DINSSTm_nu)$mu[1],coefAll(DINSSTm_nu)$mu[2], exp(coefAll(DINSSTm_nu)$sigma[1]),
                 exp(coefAll(DINSSTm_nu)$sigma[2]), exp(coefAll(DINSSTm_nu)$nu[1]),exp(coefAll(DINSSTm_nu)$nu[2]), 
                 (exp(coefAll(DINSSTm_nu)$tau) + 2), "NA", AIC(DINSSTm_nu)),
  DINSSTm_tau = c(coefAll(DINSSTm_tau)$mu[1],coefAll(DINSSTm_tau)$mu[2], exp(coefAll(DINSSTm_tau)$sigma[1]),
                 exp(coefAll(DINSSTm_tau)$sigma[2]), exp(coefAll(DINSSTm_tau)$nu[1]),exp(coefAll(DINSSTm_tau)$nu[2]), 
                 (exp(coefAll(DINSSTm_tau)$tau[1]) + 2), (exp(coefAll(DINSSTm_tau)$tau[2]) + 2), AIC(DINSSTm_tau))
)
print(param_summary)



#======================
# SEASONALITY
#======================

# Check missing data for each month  
DIN_NA <- DIN %>%
  mutate(Year = year(Date),
         Month = factor(month(Date), levels = 1:12, labels = month.abb)) %>%
  group_by(Year, Month) %>%
  summarize(DIN = mean(DIN, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month, values_from = DIN) %>%
  arrange(Year)

# Average values of DIN per months each year 
DIN_monthly <- DIN %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(DIN = mean(DIN, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year, Month)

# Add missing data (Jan,Feb,Mar of 1962; May,Jun of 1965, Jan of 1970, Oct,Nov,Dec of 1973)
DIN_monthly <- DIN_monthly %>% 
  add_row(Year = 1962, Month = 1, DIN = 15, .before = 1) %>% 
  add_row(Year = 1962, Month = 2, DIN = 22.6, .after = 1) %>% 
  add_row(Year = 1962, Month = 3, DIN = 28.2, .after = 2) %>% 
  add_row(Year = 1965, Month = 5, DIN = 32.6, .after = 40) %>%
  add_row(Year = 1965, Month = 6, DIN = 27.9, .after = 41) %>% 
  add_row(Year = 1970, Month = 1, DIN = 17, .after = 96) %>% 
  add_row(Year = 1973, Month = 10, DIN = 14.62, .after = 141) %>% 
  add_row(Year = 1973, Month = 11, DIN = 15.22, .after = 142) %>% 
  add_row(Year = 1973, Month = 12, DIN = 15.89, .after = 143) 

# turn dataset into vector 
DIN_monthly_v <- ts(
  DIN_monthly$DIN,
  start = c(min(DIN_monthly$Year), min(DIN_monthly$Month)),
  frequency = 12
)

# check that shit out 
head(DIN_monthly_v, 396)

# STL decomposition (Seasonal and Trend decomposition using Loess)
stl(DIN_monthly_v, s.window = 7) %>%  #== LOOK INTO S.WINDOW ==#
  autoplot()

# split data into training and testing set 
DINM <- series_to_mvgam(DIN_monthly_v, freq = frequency(DIN_monthly_v))
dplyr::glimpse(DINM$data_train)

# plot them 
plot_mvgam_series(data = DINM$data_train,
                  newdata = DINM$data_test)


#####################################################
### Silicate and fit models of different families ###
#####################################################

# create own df to remove NAs 
Silicate <- df %>% filter(!is.na(Silicate)) %>% 
  select(Date, Silicate, year)

# plot that sucka 
ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "grey") +
  geom_smooth(method = NULL) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "black") +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()

# plot a for a specific year 
ggplot(dplyr::filter(Silicate, year(Date) == 1972),
       aes(Date, Silicate)) +
  geom_line(color = "black")

# rescale 
Silicate$y_rescaled <- Silicate$Silicate + 100

### CHECK DISTRIBUTION ###

# histogram 
hist(Silicate$Silicate)

# look at the density plot 
ggplot(Ammonium, aes(y = Ammonium)) +
  geom_density(
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### MODELS ###
library(modeest)
# intercept model (there are 0s in the data, that SST doesn't like)

# siSSTm1 <- gamlss(y_rescaled ~ 1, family = SST(), data = Silicate, 
#                  mu.start = mean(Silicate$y_rescaled), sigma.start = max(sd(Silicate$y_rescaled), 1),
#                  method = mixed(10, 200),
#                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(siSSTm1)

# Intercept model with SEP3 => mu is mode, not mean
siSEPm1 <- gamlss(Silicate ~ 1, family = SEP3(), data = Silicate, 
                  mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                 method = mixed(20, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(siSEPm1)
mfv(Silicate$Silicate)


siSASm <- gamlss(Silicate ~ 1, family = SHASHo(), data = Silicate, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(siSASm)


siSSTm <- gamlss(Silicate ~ Date, family = SST(), data = Silicate, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

siTF2m <- gamlss(Silicate ~ Date, family = TF2(), data = Silicate, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

siSN1m <- gamlss(Silicate ~ Date, family = SN1(), data = Silicate, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

siNOm <- gamlss(Silicate ~ Date, family = NO(), data = Silicate, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


summary(siSSTm)
AIC(siSSTm)
AIC(siTF2m)
AIC(siSN1m)
AIC(siNOm)


#### CHANGING PARAMTER AS FUNCTION OF TIME WITH SST MODEL ####

siSSTm_1 <- gamlss(Silicate ~ Date, family = SST(), data = Silicate, 
                    method = mixed(5, 200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


siSSTm_sig <- gamlss(Silicate ~ Date, sigma.fo = ~ Date, family = SST(), data = Silicate, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


siSSTm_nu <- gamlss(Silicate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Silicate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


siSSTm_tau <- gamlss(Silicate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), data = Silicate, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


# summary
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  siSSTm_1 = c(coefAll(siSSTm_1)$mu[1],coefAll(siSSTm_1)$mu[2] ,exp(coefAll(siSSTm_1)$sigma),"NA", 
                exp(coefAll(siSSTm_1)$nu),"NA", (exp(coefAll(siSSTm_1)$tau) + 2), "NA", AIC(siSSTm_1)),
  siSSTm_sig = c(coefAll(siSSTm_sig)$mu[1],coefAll(siSSTm_sig)$mu[2] ,exp(coefAll(siSSTm_sig)$sigma[1]),
                  exp(coefAll(siSSTm_sig)$sigma[2]), exp(coefAll(siSSTm_sig)$nu),"NA", 
                  (exp(coefAll(siSSTm_sig)$tau) + 2), "NA", AIC(siSSTm_sig)),
  siSSTm_nu = c(coefAll(siSSTm_nu)$mu[1],coefAll(siSSTm_nu)$mu[2], exp(coefAll(siSSTm_nu)$sigma[1]),
                 exp(coefAll(siSSTm_nu)$sigma[2]), exp(coefAll(siSSTm_nu)$nu[1]),exp(coefAll(siSSTm_nu)$nu[2]), 
                 (exp(coefAll(siSSTm_nu)$tau) + 2), "NA", AIC(siSSTm_nu))
)
print(param_summary)



#======================
# SEASONALITY
#======================

# Check missing data for each month  
Silicate_NA <- Silicate %>%
  mutate(Year = year(Date),
         Month = factor(month(Date), levels = 1:12, labels = month.abb)) %>%
  group_by(Year, Month) %>%
  summarize(Silicate = mean(Silicate, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month, values_from = Silicate) %>%
  arrange(Year)

# Average values of Silicate per months each year 
Silicate_monthly <- Silicate %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(Silicate = mean(Silicate, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year, Month)

# Add missing data (Jan,Feb,Mar of 1962; May,Jun of 1965, Jan of 1970, Oct,Nov,Dec of 1973)
Silicate_monthly <- Silicate_monthly %>% 
  add_row(Year = 1966, Month = 1, Silicate = 6.2, .before = 1) %>% 
  add_row(Year = 1969, Month = 6, Silicate = 2.21, .after = 41) %>% 
  add_row(Year = 1969, Month = 7, Silicate = 1.26, .after = 42) %>% 
  add_row(Year = 1969, Month = 8, Silicate = 2, .after = 43) %>% 
  add_row(Year = 1975, Month = 3, Silicate = 9.55, .after = 110) %>% 
  add_row(Year = 1976, Month = 1, Silicate = 16.25, .after = 120) %>% 
  add_row(Year = 1976, Month = 2, Silicate = 8.22, .after = 121) %>% 
  add_row(Year = 1976, Month = 3, Silicate = 9.63, .after = 122) 

# turn dataset into vector 
Silicate_monthly_v <- ts(
  Silicate_monthly$Silicate,
  start = c(min(Silicate_monthly$Year), min(Silicate_monthly$Month)),
  frequency = 12
)

# check that shit out 
head(Silicate_monthly_v, 396)

# STL decomposition (Seasonal and Trend decomposition using Loess)
stl(Silicate_monthly_v, s.window = 7) %>%  #== LOOK INTO S.WINDOW ==#
  autoplot()

# split data into training and testing set 
SM <- series_to_mvgam(Silicate_monthly_v, freq = frequency(Silicate_monthly_v))
dplyr::glimpse(SM$data_train)

# plot them 
plot_mvgam_series(data = SM$data_train,
                  newdata = SM$data_test)



#####################################################
### Ammonium and fit models of different families ###
#####################################################

# create own df to remove NAs 
Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium, year)

# plot that sucka
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey") +
  geom_abline(intercept = amSSTm2$mu.coefficients[1], slope = amSSTm2$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # only mean as function of time 
  geom_abline(intercept = amSSTm3$mu.coefficients[1], slope = amSSTm3$mu.coefficients[2], 
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu as function of time
  geom_hline(yintercept = 5.8704, linetype = "dashed", color = "black") +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()


# plot a for a specific year 
ggplot(dplyr::filter(Ammonium, year(Date) == 1970),
       aes(Date, Ammonium)) +
  geom_line(color = "black")


### CHECK DISTRIBUTION ###

# histogram
hist(Ammonium$Ammonium)

# look at the density plot 
ggplot(Ammonium, aes(y = Ammonium)) +
  geom_density(
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### MODELS ###

# intercept model 
amSSTm1 <- gamlss(Ammonium ~ 1, family = SST(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(amSSTm1)

# mean as a function of time 
amSSTm2 <- gamlss(Ammonium ~ Date, family = SST(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

# mean, sigma and nu as a function of time 
amSSTm3 <- gamlss(Ammonium ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))



amTF2m <- gamlss(Ammonium ~ Date, family = TF2(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

amSN1m <- gamlss(Ammonium ~ Date, family = SN1(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

amNOm <- gamlss(Ammonium ~ Date, family = NO(), data = Ammonium, 
                method = mixed(5, 200),
                control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


AIC(amSSTm)
AIC(amTF2m)
AIC(amSN1m)
AIC(amNOm)


#======================
# SEASONALITY
#======================

# Check missing data for each month  
Ammonium_NA <- Ammonium %>%
  mutate(Year = year(Date),
         Month = factor(month(Date), levels = 1:12, labels = month.abb)) %>%
  group_by(Year, Month) %>%
  summarize(Ammonium = mean(Ammonium, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Month, values_from = Ammonium) %>%
  arrange(Year)

# Average values of Ammonium per months each year 
Ammonium_monthly <- Ammonium %>%
  mutate(Year = year(Date), Month = month(Date)) %>%
  group_by(Year, Month) %>%
  summarize(Ammonium = mean(Ammonium, na.rm = TRUE), .groups = "drop") %>%
  arrange(Year, Month)

# Add missing data (May and June of 1965)
Ammonium_monthly <- Ammonium_monthly %>% 
  add_row(Year = 1965, Month = 5, Ammonium = 10.2, .after = 40) %>% 
  add_row(Year = 1965, Month = 6, Ammonium = 10.3, .after = 41) 

# turn dataset into vector 
Ammonium_monthly_v <- ts(
  Ammonium_monthly$Ammonium,
  start = c(min(Ammonium_monthly$Year), min(Ammonium_monthly$Month)),
  frequency = 12
)

# check that shit out 
head(Ammonium_monthly_v, 396)

# STL decomposition (Seasonal and Trend decomposition using Loess)
stl(Ammonium_monthly_v, s.window = 7) %>%  #== LOOK INTO S.WINDOW ==#
  autoplot()

# split data into training and testing set 
AM <- series_to_mvgam(Ammonium_monthly_v, freq = frequency(Ammonium_monthly_v))
dplyr::glimpse(AM$data_train)

# plot them 
plot_mvgam_series(data = AM$data_train,
                  newdata = AM$data_test)


