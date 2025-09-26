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
library(modeest)
library(e1071)
library(mvgam)           # Fit, interrogate and forecast DGAMs
library(forecast)        # Construct fourier terms for time series
library(gratia)          # Graceful plotting of smooth terms
library(marginaleffects) # Interrogating regression models
library(janitor)         # Creating clean, tidy variable names


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

target_year <- 1963

df_one <- df %>%
  mutate(year = year(Date), month = month(Date, label = TRUE)) %>%
  filter(year == target_year) %>%
  group_by(month) %>%
  summarise(across(all_of(c(Nitrate, phyto_col)), ~mean(., na.rm = TRUE)),
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
  select(Date, Nitrate, year, month)

# plot that sucka (intercept model)
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "azure4", na.rm = TRUE) +
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
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(Nitrate$Nitrate), linetype = "dashed", color = "azure4") +
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

niTF2m <- gamlss(Nitrate ~ 1, family = TF2(), data = Nitrate, 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niSN1m <- gamlss(Nitrate ~ 1, family = SN1(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niNOm <- gamlss(Nitrate ~ 1, family = NO(), data = Nitrate, 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))



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


#############

# V = 1 --> symmetrical 
# V > 1 --> positively skewed 
# 0 < V < 1 --> negatively skewed 
# Decreasing T increases the heaviness of the tail
# When V = 1, = TF2 distribution 
# As V reaches infinity, it tends to be half t distribution TF
# As T reaches infinity, it tends to Skew normal type 2, SN2
# If both V and T reaches infinity, it tends to a half normal distribution

#############

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

#############


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

# normal distribution
niNOm1 <- gamlss(Nitrate ~ Date, family = NO(), data = Nitrate, trace = FALSE)
summary(niNOm1)


niSSTpoly1 <- gamlss(Nitrate ~ poly(Date, 2), sigma.fo = ~ poly(Date, 2), nu.fo = ~ poly(Date, 2), family = SST(), data = Nitrate,
                  #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                  method = mixed(10,200),
                  control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTpoly1)



# mean, sigma and nu and tau changing through time linearly 
# niSSTm4 <- gamlss(Nitrate ~ t_scaled, sigma.fo = ~ t_scaled, nu.fo = ~ t_scaled, tau.fo = ~ t_scaled, family = SST(), data = Nitrate,
#                   mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
#                   method = mixed(10,200),
#                   control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))
# summary(niSSTm4)


## ADDING SEASONAILITY ##

# mean only
niSSTm2_sea <- gamlss(Nitrate ~ year + month, family = SST(), data = Nitrate,
                     #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm2_sea)

# mu sig and nu changing
niSSTm3_sea <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, family = SST(), data = Nitrate,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niSSTm3_sea)

#============================================================================================

# predict mu based on seasonaility model
Nitrate$mu_hat <- fitted(niSSTm3_sea, what = "mu")

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

# plot that sucka, with seasonaility 
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "azure4", na.rm = TRUE) +
  geom_line(aes(y = mu_hat), color = "darkred", linewidth = 1) +
  theme_minimal()

################################################################### PDF 

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


##################################################################################
# function for pdf
pdf_NO_at_date <- function(model, data, date_str, x = NULL, n = 200,
                            time_var = "Date", date_var = "Date") {
  # 1. Find the row in df matching the date of interest
  t_val <- data[[time_var]][as.Date(data[[date_var]]) == as.Date(date_str)]
  if (length(t_val) == 0) stop("Date not found in data frame.")
  
  # 2. Build newdata with the correct covariate
  newdat <- data.frame(setNames(list(t_val), time_var))
  
  # 3. Predict distribution parameters on natural scale
  mu    <- predict(model, "mu",    newdata = newdat, type = "response")
  sigma <- predict(model, "sigma", newdata = newdat, type = "response")

  
  # 4. Create x grid if none supplied
  if (is.null(x)) {
    x <- seq(min(model$y, na.rm = TRUE),
             max(model$y, na.rm = TRUE),
             length.out = n)
  }
  
  # 5. Evaluate PDF
  dens <- dNO(x, mu = mu, sigma = sigma)
  
  list(x = x, density = dens,
       params = c(mu = mu, sigma = sigma))
}


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
# Normal (mean)
resm4 <- pdf_NO_at_date(niNOm1, Nitrate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm4b <- pdf_NO_at_date(niNOm1, Nitrate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm4c <- pdf_NO_at_date(niNOm1, Nitrate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# poly (mean)
resm5 <- pdf_SST_at_date(niSSTpoly1, Nitrate, "1963-09-23",
                        time_var = "Date", date_var = "Date")
resm5b <- pdf_SST_at_date(niSSTpoly1, Nitrate, "1980-09-23",
                         time_var = "Date", date_var = "Date")
resm5c <- pdf_SST_at_date(niSSTpoly1, Nitrate, "1994-04-13",
                         time_var = "Date", date_var = "Date")
# season
resm6 <- pdf_SST_at_date(niSSTm2_sea, Nitrate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm6b <- pdf_SST_at_date(niSSTm2_sea, Nitrate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm6c <- pdf_SST_at_date(niSSTm2_sea, Nitrate, "1994-04-13",
                          time_var = "Date", date_var = "Date")


# pdf for seasonaility model
x <- seq(min(niSSTm3_sea$y, na.rm = TRUE),
         max(niSSTm3_sea$y, na.rm = TRUE),
         length.out = 200)

pdf_1963 <- dSST(x, mu = 1.2156707, sigma = 3.290147, nu = 9.663093, tau = 6.874041)
pdf_1980 <- dSST(x, mu = 8.671901, sigma = 6.657359, nu = 7.422281, tau = 6.874041)
pdf_1994 <- dSST(x, mu = 31.41314, sigma = 24.87263, nu = 0.9332162, tau = 6.874041)



# plot all three dates 

par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.15),
     ylab = "Density", xlab = "Value",
     main = "1963-09-23") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3$x, resm3$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(x, pdf_1963, col = "darkolivegreen", lwd = 2, lty = 1)# seasonaility model
# lines(resm4$x, resm4$density, col = "darkolivegreen",lwd = 2, lty = 1) # NORMAL 
#lines(resm5$x, resm5$density, col = "salmon",lwd = 2, lty = 1) # poly 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.15),
     ylab = "Density", xlab = "Value",
     main = "1980-09-03") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3b$x, resm3b$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma an nu 
lines(x, pdf_1980, col = "darkolivegreen", lwd = 2, lty = 1)# seasonaility model
#lines(resm4b$x, resm4b$density, col = "darkolivegreen",lwd = 2, lty = 1) # NORMAL 
#lines(resm5b$x, resm5b$density, col = "salmon",lwd = 2, lty = 1) # poly 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
# plot the models together for date 3: 1994-04-13
plot(resm1c$x, resm1c$density, col = "black", type = "l", lwd = 2, lty =2, ylim=c(0,0.15),
     ylab = "Density", xlab = "Value",
     main = "1994-04-13") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3c$x, resm3c$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(x, pdf_1994, col = "darkolivegreen", lwd = 2, lty = 1)# seasonaility model
#lines(resm4c$x, resm4c$density, col = "darkolivegreen",lwd = 2, lty = 1) # NORMAL 
#lines(resm5c$x, resm5c$density, col = "salmon",lwd = 2, lty = 1) # poly 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
par(mfrow = c(1,1)) # reset




### parameters for the three dates ###

# get the predicted param for seasonaility model
Nitrate$sigma_hat <- fitted(niSSTm3_sea, what = "sigma", type = "response")
Nitrate$nu_hat <- fitted(niSSTm3_sea, what = "nu", type = "response")
Nitrate$tau_hat <- fitted(niSSTm3_sea, what = "tau", type = "response")


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






################################################################### Specific year PDF

Nitrate_1962 <- Nitrate %>%
  filter(year(Date) == 1962)

Nitrate_1963 <- Nitrate %>%
  filter(year(Date) == 1963)

Nitrate_1980 <- Nitrate %>%
  filter(year(Date) == 1980)

Nitrate_1994 <- Nitrate %>%
  filter(year(Date) == 1994)


Nitrate <- Nitrate %>%
  mutate(DOY  = yday(Date))   # 1–365

years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)

ggplot(Nitrate %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Nitrate, color = factor(year))) +
  geom_line() +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Day of Year", y = "Nitrate (µmol/l)", colour = "Year") +
  theme_classic()


# look at nitrate trends per year 
ggplot(Nitrate_1963, aes(x = Date, y = Nitrate)) +
  geom_line(color = "azure4", na.rm = TRUE) +
  theme_classic()

ggplot() +
  geom_line(data = Nitrate_1963, aes(x = month, y = Nitrate), linetype = 2) +
  geom_line(data = Nitrate_1980, aes(x = month, y = Nitrate), color = "azure4") +
  geom_line(data = Nitrate_1994, aes(x = month, y = Nitrate), linetype = 4 , color = "grey") +
  theme_minimal()
  
# density plot
ggplot(Nitrate_1963, aes(y = Nitrate)) +
  geom_density(linewidth = 0.8) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  geom_hline(yintercept = mean(Nitrate_1963$Nitrate), linetype = "dashed", color = "azure4") +
  theme_minimal(base_size = 14) +
  coord_flip()

# what is geom_density doing 
m <- ggplot(Nitrate_1963, aes(x = Nitrate))
m <- m + geom_density()
p <- ggplot_build(m)
head(p$data[[1]], 3)

# pdf per year function
pdf_SST_at_year <- function(model, data, year, x = NULL, n = 200,
                            time_var = "Date", date_var = "Date") {
  # 1. Filter dates for the specified year
  all_dates <- unique(data[[date_var]][format(as.Date(data[[date_var]]), "%Y") == as.character(year)])
  if (length(all_dates) == 0) stop("Year not found in data frame.")
  
  # 2. Create common x grid
  if (is.null(x)) {
    x <- seq(min(model$y, na.rm = TRUE),
             max(model$y, na.rm = TRUE),
             length.out = n)
  }
  
  # 3. For each date, calculate params and PDF
  results <- lapply(all_dates, function(date_str) {
    t_val <- data[[time_var]][as.Date(data[[date_var]]) == as.Date(date_str)]
    if (length(t_val) == 0) return(NULL)
    newdat <- data.frame(setNames(list(t_val), time_var))
    mu    <- predict(model, "mu",    newdata = newdat, type = "response")
    sigma <- predict(model, "sigma", newdata = newdat, type = "response")
    nu    <- predict(model, "nu",    newdata = newdat, type = "response")
    tau   <- predict(model, "tau",   newdata = newdat, type = "response")
    dens  <- dSST(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
    list(params = c(mu = mu, sigma = sigma, nu = nu, tau = tau), density = dens)
  })
  
  # Remove NULLs if any
  results <- Filter(Negate(is.null), results)
  
  # 4. Average parameters and density
  params_mat <- do.call(rbind, lapply(results, function(res) res$params))
  avg_params <- colMeans(params_mat, na.rm = TRUE)
  
  dens_mat <- do.call(rbind, lapply(results, function(res) res$density))
  avg_density <- colMeans(dens_mat, na.rm = TRUE)
  
  list(
    x = x,
    avg_density = avg_density,
    avg_params = avg_params,
    all_params = params_mat
  )
}


# intercept model
ni1963_1 <- pdf_SST_at_year(niSSTm, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni1980_1 <- pdf_SST_at_year(niSSTm, Nitrate, "1980", time_var = "Date", date_var = "Date")
ni1994_1 <- pdf_SST_at_year(niSSTm, Nitrate, "1994", time_var = "Date", date_var = "Date")
# mean only model
ni1963_2 <- pdf_SST_at_year(niSSTm2, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni1980_2 <- pdf_SST_at_year(niSSTm2, Nitrate, "1980", time_var = "Date", date_var = "Date")
ni1994_2 <- pdf_SST_at_year(niSSTm2, Nitrate, "1994", time_var = "Date", date_var = "Date")
# all param changing model (not tau)
ni1963_3 <- pdf_SST_at_year(niSSTm3, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni1980_3 <- pdf_SST_at_year(niSSTm3, Nitrate, "1980", time_var = "Date", date_var = "Date")
ni1994_3 <- pdf_SST_at_year(niSSTm3, Nitrate, "1994", time_var = "Date", date_var = "Date")


par(mfrow = c(1, 3)) 
# plot the models together for year 1993
plot(ni1963_1$x, ni1963_1$avg_density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.08),
     ylab = "Density", xlab = "Value",
     main = "1963") # intercept model
lines(ni1963_2$x, ni1963_2$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(ni1963_3$x, ni1963_3$avg_density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
# plot the models together for year 1980
plot(ni1980_1$x, ni1980_1$avg_density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.08),
     ylab = "Density", xlab = "Value",
     main = "1980") # intercept model
lines(ni1980_2$x, ni1980_2$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(ni1980_3$x, ni1980_3$avg_density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
# plot the models together for year 1994
plot(ni1994_1$x, ni1994_1$avg_density, col = "black", type = "l", lwd = 2, lty =2, ylim=c(0,0.08),
     ylab = "Density", xlab = "Value",
     main = "1994") # intercept model
lines(ni1994_2$x, ni1994_2$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(ni1994_3$x, ni1994_3$avg_density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
par(mfrow = c(1,1)) # reset



# check heavy tails
mu0 <- 16.717847; sig0 <- 15.667711; nu0 <- 9.970027; tau0 <- 4.468998
mu1 <- 15.530967; sig1 <- 15.344068; nu1 <- 7.176139; tau1 <- 4.346268
mu2 <- 7.058096; sig2 <- 7.168460; nu2 <- 4.863863; tau2 <- 10.096813

k <- 3
pZ <- c(Intercept = 1 - pSST(mu0 + sig0*k, mu0, sig0, nu0, tau0),
        Meanonly = 1 - pSST(mu1 + sig1*k, mu1, sig1, nu1, tau1),
        TimeVarying = 1 - pSST(mu2 + sig2*k, mu2, sig2, nu2, tau2))
pZ
## timeVarying has higher probability, means extreme values are more often than other models

###################################################################

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

###################################################################


################################################## Moment bucket 
# intercept model 
moment_bucket(niSSTm) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean only 
moment_bucket(niSSTm2) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(niSSTm3) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

###################################################################


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


##################################################################################

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


##################################################################################





######################################################
### Phosphate and fit models of different families ###
######################################################

# create own df to remove NAs
Phosphate <- df %>% filter(!is.na(Phosphate)) %>% 
  select(Date, Phosphate, year, month)


ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "azure4") +
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
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(Phosphate$Phosphate), linetype = "dashed", color = "azure4") +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


density_estimate <- density(Phosphate$Phosphate)

# Plot the estimated density
plot(density_estimate, main = "Estimated Probability Density Function",
     xlab = "Value", ylab = "Density")





### MODELS ###

phSSTm <- gamlss(Phosphate ~ poly(Date,2), family = SST(), data = Phosphate, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phTF2m <- gamlss(Phosphate ~ Date, family = TF2(), data = Phosphate, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phSN1m <- gamlss(Phosphate ~ Date, family = SN1(), data = Phosphate, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phNOm <- gamlss(Phosphate ~ Date, family = NO(), data = Phosphate, 
                method = mixed(5, 200),
                control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


#========================================================================================
# # Intercept model (different family tho ;( ) => mu is mode
# phSEPm1 <- gamlss(Phosphate ~ 1, family = SEP3(), data = Phosphate,
#                  method = mixed(5, 200),
#                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(phSEPm1)
# mfv(Phosphate$Phosphate) # find the mode
# 
# # mean only
# phSEPm2 <- gamlss(Phosphate ~ Date, family = SEP3(), data = Phosphate,
#                   method = mixed(5, 200),
#                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(phSEPm2)
# 
# # mean sigma and nu changing through time
# phSEPm3 <- gamlss(Phosphate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SEP3(), data = Phosphate,
#                   method = mixed(5, 200),
#                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(phSEPm3)
#========================================================================================


# TEST
#============================================================================================
# intercept model 

### LINEAR ###

## just mean
#  ~ Date
# ~ year 
# ~ year + month 
# ~ year + month + DOY

## all param basic
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
# ~ year, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date


# intercept 
ph_Intercept <- gamlss(Phosphate ~ 1, family = JSU(), data = Phosphate,
                    mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

## just mean
#  ~ Date
ph_mean_date <- gamlss(Phosphate ~ Date, family = JSU(), data = Phosphate,
                    mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year 
ph_mean_year <- gamlss(Phosphate ~ year, family = JSU(), data = Phosphate,
                    mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month
ph_mean_year_month <- gamlss(Phosphate ~ year + month, family = JSU(), data = Phosphate,
                          mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY
ph_mean_year_month_DOY <- gamlss(Phosphate ~ year + month + DOY, family = JSU(), data = Phosphate,
                              mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(ph_mean_date)
AIC(ph_mean_year)
AIC(ph_mean_year_month)
AIC(ph_mean_year_month_DOY)


## all param basic 
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
ph_param_date <- gamlss(Phosphate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date,  tau.fo = ~ Date, family = JSU(), data = Phosphate,
                     #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
ph_param_year <- gamlss(Phosphate ~ year, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, family = JSU(), data = Phosphate,
                     #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
ph_param_year_month <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                              family = JSU(), data = Phosphate,
                           #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY
ph_param_year_month_DOY <- gamlss(Phosphate ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY, tau.fo = ~ year + month + DOY,
                                  family = JSU(), data = Phosphate,
                                     mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                     method = mixed(10,200),
                                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(ph_param_date)
AIC(ph_param_year)
AIC(ph_param_year_month)
AIC(ph_mean_year_month_DOY) 

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
ph_param_date_sw_year <- gamlss(Phosphate ~ Date, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year,
                                family = JSU(), data = Phosphate,
                             mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
ph_param_year_month_sw_year <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year,
                                      family = JSU(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
ph_param_year_month_sw_date <- gamlss(Phosphate ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                      family = JSU(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
ph_param_year_month_DOY_sw_year <- gamlss(Phosphate ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year,
                                          family = JSU(), data = Phosphate,
                                       #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
ph_param_year_month_DOY_sw_year_month <- gamlss(Phosphate ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                                                family = JSU(), data = Phosphate,
                                             #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                             method = mixed(10,200),
                                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date
ph_param_year_month_DOY_sw_date <- gamlss(Phosphate ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                          family = JSU(), data = Phosphate,
                                       #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


AIC(ph_param_date_sw_year) 
AIC(ph_param_year_month_sw_year)
AIC(ph_param_year_month_sw_date)
AIC(ph_param_year_month_DOY_sw_year)
AIC(ph_param_year_month_DOY_sw_year_month)
AIC(ph_param_year_month_DOY_sw_date)

summary(ph_param_year_month) # BEST MODEL
summary(ph_param_year_month_DOY_sw_year_month) # SECOND BEST MODEL
summary(ph_param_year_month_sw_year) # THIRD BEST MODEL

# add polynomial to best model 
poly_ph_param_year_month <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                              family = JSU(), data = Phosphate,
                              #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
AIC(poly_ph_param_year_month) # Better?? negative AIC
summary(poly_ph_param_year_month)

# compare different families for the best model 
SHASHo_poly_ph_param_year_month <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                   family = SHASHo(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

NET_poly_ph_param_year_month <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                   family = NET(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

SEP3_poly_ph_param_year_month <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                   family = SEP3(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

SST_poly_ph_param_year_month <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = SST(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(SHASHo_poly_ph_param_year_month)
AIC(NET_poly_ph_param_year_month)
AIC(SEP3_poly_ph_param_year_month)
AIC(SST_poly_ph_param_year_month) # tau is only ~ year
#=======================================================================================

#=======================================================================================
# intercept model
phJSUm1 <- gamlss(Phosphate ~ 1, family = JSU(), data = Phosphate,
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phJSUm1)

# mean only
phJSUm2 <- gamlss(Phosphate ~ Date, family = JSU(), data = Phosphate,
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phJSUm2)

# mean sigma and nu changing through time
phJSUm3 <- gamlss(Phosphate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = JSU(), data = Phosphate,
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phJSUm3)

# mean sigma nu and tau changing through time
phJSUm4 <- gamlss(Phosphate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = JSU(), data = Phosphate,
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phJSUm4)

# mean sigma nu and tau changing through time with poly
phJSUpolym1 <- gamlss(Phosphate ~ poly(Date,2), sigma.fo = ~ poly(Date,2), nu.fo = ~ poly(Date,2), tau.fo = ~ poly(Date,2), family = JSU(), data = Phosphate,
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phJSUpolym1)


# only mean poly
phJSUpolym2 <- gamlss(Phosphate ~ poly(Date,2), sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = JSU(), data = Phosphate,
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(phJSUpolym2)


# mu sig and nu changing with seasonality
phJSUm4_sea <- gamlss(Phosphate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                      family = SST(), data = Phosphate,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(phJSUm4_sea)
#========================================================================================
# predict mu based on seasonality model
Phosphate$mu_hat <- fitted(ph_param_year_month, what = "mu")
Phosphate$mu_hat2 <- fitted(poly_ph_param_year_month, what = "mu")

# plot that sucka
ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
  geom_line(color = "grey") +
  geom_abline(intercept = ph_mean_date$mu.coefficients[1], slope = ph_mean_date$mu.coefficients[2],
              color = "steelblue", linewidth = 1) + #mean only
  #geom_abline(intercept = phJSUm3$mu.coefficients[1], slope = phJSUm3$mu.coefficients[2],
  #color = "goldenrod", linewidth = 1) + # mean sigma and nu changing through time
  geom_abline(intercept = ph_param_date$mu.coefficients[1], slope = ph_param_date$mu.coefficients[2],
              color = "goldenrod", linewidth = 1) + # mean sigma nu and tau changing through time
  geom_line(aes(y = mu_hat2), color = "darkolivegreen", linewidth = 1) + 
  geom_hline(yintercept = 0.75230, linetype = "dashed", color = "black") + 
  geom_vline(xintercept = as.Date("1963-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1980-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1994-04-13"), linetype = "dashed", color = "darkred") +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()

############################################ PDF

# get the pdf and param for three dates
# m1 (intercept model)
resm1 <- pdf_SST_at_date(phJSUm1, Phosphate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm1b <- pdf_SST_at_date(phJSUm1, Phosphate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm1c <- pdf_SST_at_date(phJSUm1, Phosphate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m2 (mean only)
resm2 <- pdf_SST_at_date(phJSUm2, Phosphate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm2b <- pdf_SST_at_date(phJSUm2, Phosphate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm2c <- pdf_SST_at_date(phJSUm2, Phosphate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m3 (mean sigma nu and tau)
resm3 <- pdf_SST_at_date(phJSUm4, Phosphate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm3b <- pdf_SST_at_date(phJSUm4, Phosphate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm3c <- pdf_SST_at_date(phJSUm4, Phosphate, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# poly (mean sigma nu and tau)
resm4 <- pdf_SST_at_date(phJSUpolym2, Phosphate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm4b <- pdf_SST_at_date(phJSUpolym2, Phosphate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm4c <- pdf_SST_at_date(phJSUpolym2, Phosphate, "1994-04-13",
                          time_var = "Date", date_var = "Date")



par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,2.5),
     ylab = "Density", xlab = "Value",
     main = "1963") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3$x, resm3$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(resm4$x, resm4$density, col = "darkolivegreen",lwd = 2, lty = 1) # poly
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model", "poly"),
       col = c("black", "steelblue", "goldenrod", "darkolivegreen"), cex = 0.5, lty = c(2,1,1,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,2.3),
     ylab = "Density", xlab = "Value",
     main = "1980") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3b$x, resm3b$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(resm4b$x, resm4b$density, col = "darkolivegreen",lwd = 2, lty = 1) # poly
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model", "poly"),
       col = c("black", "steelblue", "goldenrod", "darkolivegreen"), cex = 0.5, lty = c(2,1,1,1), bty = "n")
# plot the models together for date 3: 1994-04-13
plot(resm1c$x, resm1c$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,2.3),
     ylab = "Density", xlab = "Value",
     main = "1994") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3c$x, resm3c$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
lines(resm4c$x, resm4c$density, col = "darkolivegreen",lwd = 2, lty = 1) # poly
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model", "poly"),
       col = c("black", "steelblue", "goldenrod", "darkolivegreen"), cex = 0.5, lty = c(2,1,1,1), bty = "n")
par(mfrow = c(1, 1)) # reset


# # plot all three dates 
# par(mfrow = c(1, 3))
# plot(res1$x, res1$density, type = "l", ylim=c(0,2.3)) # 1.5 for m1 and m2
# plot(res2$x, res2$density, type = "l", ylim=c(0,2.3)) # 1.55 for m3
# plot(res3$x, res3$density, type = "l", ylim=c(0,2.3)) # 2.3 for m4
# par(mfrow = c(1, 1)) # reset


### parameters for the three dates ###

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
# poly (mean sigma and nu)
resm4$params #1963
resm4b$params #1980
resm4c$params #1994

#######################################################

###### Look at moment changes every n years ########


summary <- Phosphate %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Phosphate, na.rm = TRUE),
    var = var(Phosphate, na.rm = TRUE),
    skew = skewness(Phosphate, na.rm = TRUE),
    kurt = kurtosis(Phosphate, na.rm = TRUE)
  )

print(summary, n=33)

###### Look at moment changes average months ########

month_moments <- Phosphate %>%
  mutate(
    month_num = month(Date),
    month_lab = month(Date, label = TRUE, abbr = TRUE)
  ) %>%
  group_by(month_num, month_lab) %>%
  summarise(
    n     = sum(!is.na(Phosphate)),
    mean  = mean(Phosphate, na.rm = TRUE),
    var   = var(Phosphate, na.rm = TRUE),
    skew  = skewness(Phosphate, na.rm = TRUE),
    kurt  = kurtosis(Phosphate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_num)

month_moments

# compare the param with the predicted param of the best model 
poly_ph_param_year_month <- gamlss(Phosphate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                   family = JSU(), data = Phosphate,
                                   #mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# create function to get predicted param from model
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- fitted(model, what = "mu")
  data$sigma_hat <- fitted(model, what = "sigma")
  data$nu_hat <- fitted(model, what = "nu")
  data$tau_hat <- fitted(model, what = "tau")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
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

model_pred_param <- best_model_param(poly_ph_param_year_month, Phosphate)
print(model_pred_param$by_year, n=33)
print(model_pred_param$by_month)


## PLOT comparison per year 

# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,1), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "mean")
lines(summary$mean ~ summary$year, col = "azure4", lwd = 2, lty = 2)

# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,0.5), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "variance")
lines(summary$var ~ summary$year, col = "azure4", lwd = 2, lty = 2)

# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$year, type = "l", 
    ylim = c(-0.4, 7), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "skewness")
lines(summary$skew ~ summary$year, col = "azure4", lwd = 2, lty = 2)
abline(h = 0, lty = 3, col = "red")

# plot empirical tau vs predicted tau
plot(model_pred_param$mean_tau_hat ~ model_pred_param$year, type = "l", 
     ylim = c(1, 7), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "kurtosis")
lines(summary$kurt ~ summary$year, col = "azure4", lwd = 2, lty = 2)


## PLOT comparison per month 
# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0,1.2), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "mean")
lines(month_moments$mean ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0,0.5), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "variance")
lines(month_moments$var ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$month, type = "l", 
     ylim = c(-10, 19), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "skew")
lines(month_moments$skew ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)
abline(h = 0, lty = 3, col = "red")

# plot empirical tau vs predicted tau
plot(model_pred_param$mean_tau_hat ~ model_pred_param$month, type = "l", 
     ylim = c(1, 12), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "kurtosis")
lines(month_moments$kurt ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

################################################## Moment bucket 
# intercept model 
moment_bucket(phJSUm1) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean only 
moment_bucket(phJSUm2) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(phJSUm3) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(phJSUpolym1) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

###################################################################


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



############ TRENDS ###############


# create day of the year for nitrate and phytopl
Phosphate <- Phosphate %>%
  mutate(DOY  = yday(Date)) 

# look at phosphate patterns for different years 
years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)

ggplot(Phosphate %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Phosphate, color = factor(year))) +
  geom_line() +
  scale_colour_manual(values = c("Other" = "lightgrey",
                                 "1963" = "chocolate",
                                 "1980" = "cornflowerblue",
                                 "1994" = "darkseagreen4")) + 
  geom_vline(xintercept = 355, linetype = "dashed") + # start of winter
  geom_vline(xintercept = 80, linetype = "dashed") + # start of spring
  geom_vline(xintercept = 173, linetype = "dashed") + # start of summer
  geom_vline(xintercept = 266, linetype = "dashed") + # start of fall
  annotate("text", x = 37, y = 3, label = "Winter", size = 4) +
  annotate("text", x = 130, y = 3, label = "Spring", size = 4) +
  annotate("text", x = 220, y = 3, label = "Summer", size = 4) +
  annotate("text", x = 315, y = 3, label = "Fall", size = 4) +
  labs(x = "Day of Year", y = "Phosphate (µmol/l)", colour = "Year") +
  theme_classic()

####################################################
### Nitrite and fit models of different families ###
####################################################

# Create own df to remove NAs
Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite, year, month)

ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "darkgrey") + 
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
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(Nitrite$Nitrite), linetype = "dashed", color = "azure4") +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### Other family distribution MODELS ###

niiTF2m <- gamlss(Nitrite ~ Date, family = TF2(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niiSN1m <- gamlss(Nitrite ~ Date, family = SN1(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niiNOm <- gamlss(Nitrite ~ Date, family = NO(), data = Nitrite, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))



# TEST
#============================================================================================
# intercept model 

### LINEAR ###

## just mean
#  ~ Date
# ~ year 
# ~ year + month 
# ~ year + month + DOY

## all param basic
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
# ~ year, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date


# intercept 
ni_Intercept <- gamlss(Nitrite ~ 1, family = SST(), data = Nitrite,
                    mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

## just mean
#  ~ Date
ni_mean_date <- gamlss(Nitrite ~ Date, family = SST(), data = Nitrite,
                    mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year 
ni_mean_year <- gamlss(Nitrite ~ year, family = SST(), data = Nitrite,
                    mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month
ni_mean_year_month <- gamlss(Nitrite ~ year + month, family = SST(), data = Nitrite,
                          mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY
ni_mean_year_month_DOY <- gamlss(Nitrite ~ year + month + DOY, family = SST(), data = Nitrite,
                              mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(ni_mean_date)
AIC(ni_mean_year)
AIC(ni_mean_year_month)
AIC(ni_mean_year_month_DOY) # best


## all param basic 
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
ni_param_date <- gamlss(Nitrite ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                        family = SST(), data = Nitrite,
                     #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
ni_param_year <- gamlss(Nitrite ~ year, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                     family = SST(), data = Nitrite,
                     #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
ni_param_year_month <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                           family = SST(), data = Nitrite,
                           mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY
ni_param_year_month_DOY <- gamlss(Nitrite ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY, tau.fo = ~ year, 
                                  family = SST(), data = Nitrite,
                                     mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                     method = mixed(10,200),
                                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(ni_param_date)
AIC(ni_param_year)
AIC(ni_param_year_month) # tau had to be ~ year only ## best
AIC(ni_mean_year_month_DOY) #tau had to be ~ year only

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
ni_param_date_sw_year <- gamlss(Nitrite ~ Date, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year,
                                family = SST(), data = Nitrite,
                             mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
ni_param_year_month_sw_year <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year,
                                   family = SST(), data = Nitrite,
                                   #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
ni_param_year_month_sw_date <- gamlss(Nitrite ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                      family = SST(), data = Nitrite,
                                   #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
ni_param_year_month_DOY_sw_year <- gamlss(Nitrite ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                                          family = SST(), data = Nitrite,
                                       #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
ni_param_year_month_DOY_sw_year_month <- gamlss(Nitrite ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                                family = SST(), data = Nitrite,
                                             #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                             method = mixed(10,200),
                                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date
ni_param_year_month_DOY_sw_date <- gamlss(Nitrite ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                          family = SST(), data = Nitrite,
                                       #mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


AIC(ni_param_date_sw_year) 
AIC(ni_param_year_month_sw_year)
AIC(ni_param_year_month_sw_date)
AIC(ni_param_year_month_DOY_sw_year)
AIC(ni_param_year_month_DOY_sw_year_month) # tau had to be ~ year only ## best
AIC(ni_param_year_month_DOY_sw_date)

summary(ni_param_year_month) # BEST MODEL
summary(ni_param_year_month_DOY_sw_year_month) # SECOND BEST MODEL
summary(ni_param_year_month_sw_year) # THIRD BEST MODEL

# adding polynomial to best model
poly_ni_param_year_month <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                              family = SST(), data = Nitrite,
                              mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
AIC(poly_ni_param_year_month) # BETTER YOHOO
summary(poly_ni_param_year_month)

# compare different family for the best model
JSU_poly_ni_param_year_month <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = JSU(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

SHASHo_poly_ni_param_year_month <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = SHASHo(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

NET_poly_ni_param_year_month <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = NET(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


AIC(JSU_poly_ni_param_year_month)
AIC(SHASHo_poly_ni_param_year_month)
AIC(NET_poly_ni_param_year_month)


#=====================================================================================

#=====================================================================================
# Intercept model
niiSSTm1 <- gamlss(Nitrite ~ 1, family = SST(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niiSSTm1)

# Only mean changing through time 
niiSSTm2 <- gamlss(Nitrite ~ Date, family = SST(), data = Nitrite, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niiSSTm2)

# Mean, sigma and nu changing through time
niiSSTm3 <- gamlss(Nitrite ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Nitrite, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niiSSTm3)

# Mean, sigma nu and tau changing through time
niiSSTm4 <- gamlss(Nitrite ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), data = Nitrite, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niiSSTm4)

# mu sig and nu changing with seasonality
niiSSTm4_sea <- gamlss(Nitrite ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                      family = SST(), data = Nitrite,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(niiSSTm4_sea)
#=====================================================================================
# predicted mu param for seasonaility model 
Nitrite$mu_hat <- fitted(ni_param_year_month, what = "mu")
Nitrite$mu_hat2 <- fitted(poly_ni_param_year_month, what = "mu")

# plot that sucka
ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "grey") + 
  geom_abline(intercept = ni_mean_date$mu.coefficients[1], slope = ni_mean_date$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # mean only changing through time 
  geom_abline(intercept = ni_param_date$mu.coefficients[1], slope = ni_param_date$mu.coefficients[2],
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu changing through time
  geom_line(aes(y = mu_hat2), color = "darkolivegreen", linewidth = 1) +
  #geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, color = "darkolivegreen") +
  geom_vline(xintercept = as.Date("1963-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1980-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1994-04-13"), linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = 0.81659, linetype = "dashed", color = "black") +  # intercept
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()


############################################ PDF

# get the pdf and param for three dates
# m1 (intercept model)
resm1 <- pdf_SST_at_date(niiSSTm1, Nitrite, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm1b <- pdf_SST_at_date(niiSSTm1, Nitrite, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm1c <- pdf_SST_at_date(niiSSTm1, Nitrite, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m2 (mean only)
resm2 <- pdf_SST_at_date(niiSSTm2, Nitrite, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm2b <- pdf_SST_at_date(niiSSTm2, Nitrite, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm2c <- pdf_SST_at_date(niiSSTm2, Nitrite, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m3 (mean sigma nu and tau)
resm3 <- pdf_SST_at_date(niiSSTm4, Nitrite, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm3b <- pdf_SST_at_date(niiSSTm4, Nitrite, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm3c <- pdf_SST_at_date(niiSSTm4, Nitrite, "1994-04-13",
                          time_var = "Date", date_var = "Date")


par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,1.2),
     ylab = "Density", xlab = "Value",
     main = "1963-09-23") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3$x, resm3$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,1.2),
     ylab = "Density", xlab = "Value",
     main = "1980-09-23") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3b$x, resm3b$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")
# plot the models together for date 3: 1994-04-13
plot(resm1c$x, resm1c$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,1.2),
     ylab = "Density", xlab = "Value",
     main = "1994-04-13") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3c$x, resm3c$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")


### parameters for the three dates ###
# ylim: 1.2 for all

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


# # plot all three dates 
# par(mfrow = c(1, 3))
# plot(res1$x, res1$density, type = "l", ylim=c(0,1.2))
# plot(res2$x, res2$density, type = "l", ylim=c(0,1.2)) # 1.2 for all
# plot(res3$x, res3$density, type = "l", ylim=c(0,1.2))
# par(mfrow = c(1, 1)) # reset


#######################################################

###### Look at moment changes every n years ########


summary <- Nitrite %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Nitrite, na.rm = TRUE),
    var = var(Nitrite, na.rm = TRUE),
    skew = skewness(Nitrite, na.rm = TRUE),
    kurt = kurtosis(Nitrite, na.rm = TRUE)
  )

print(summary, n=33)

###### Look at moment changes average months ########

month_moments <- Nitrite %>%
  mutate(
    month_num = month(Date),
    month_lab = month(Date, label = TRUE, abbr = TRUE)
  ) %>%
  group_by(month_num, month_lab) %>%
  summarise(
    n     = sum(!is.na(Nitrite)),
    mean  = mean(Nitrite, na.rm = TRUE),
    var   = var(Nitrite, na.rm = TRUE),
    skew  = skewness(Nitrite, na.rm = TRUE),
    kurt  = kurtosis(Nitrite, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_num)

month_moments

# compare the param with the predicted param of the best model 
poly_ni_param_year_month <- gamlss(Nitrite ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = SST(), data = Nitrite,
                                   mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# create function to get predicted param from model
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- fitted(model, what = "mu")
  data$sigma_hat <- fitted(model, what = "sigma")
  data$nu_hat <- fitted(model, what = "nu")
  data$tau_hat <- fitted(model, what = "tau")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
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

model_pred_param <- best_model_param(poly_ni_param_year_month, Nitrite)
print(model_pred_param, n=33)
print(model_pred_param$by_month)


## PLOT comparison per year 

# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,1.5), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "mean")
lines(summary$mean ~ summary$year, col = "azure4", lwd = 2, lty = 2)

# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,1), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "variance")
lines(summary$var ~ summary$year, col = "azure4", lwd = 2, lty = 2)

# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0, 5), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "skewness")
lines(summary$skew ~ summary$year, col = "azure4", lwd = 2, lty = 2)
abline(h = 0, lty = 3, col = "red")

# plot empirical tau vs predicted tau
#plot(model_pred_param$mean_tau_hat ~ model_pred_param$year, type = "l", 
     #ylim = c(1, 7), col = "goldenrod", lwd = 2,
     #xlab = "year", ylab = "kurtosis")
plot(summary$kurt ~ summary$year, col = "azure4", lwd = 2, type = "l", lty = 2)
abline(h = 3, lty = 3, col = "red") # normal tails


## PLOT comparison per month 
# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0,1.5), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "mean")
lines(month_moments$mean ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0,1), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "variance")
lines(month_moments$var ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0, 7), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "skew")
lines(month_moments$skew ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical tau vs predicted tau
#plot(model_pred_param$mean_tau_hat ~ model_pred_param$month, type = "l", 
     #ylim = c(1, 12), col = "goldenrod", lwd = 1,
     #xlab = "month", ylab = "kurtosis")
plot(month_moments$kurt ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2, type = "l")
abline(h = 3, lty = 3, col = "red") # normal tails

################################################## Moment bucket 
# intercept model 
moment_bucket(niiSSTm1) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean only 
moment_bucket(niiSSTm2) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(niiSSTm3) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

###################################################################


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


############ TRENDS ###############


# create day of the year for nitrate and phytopl
Nitrite <- Nitrite %>%
  mutate(DOY  = yday(Date)) 

# look at phosphate patterns for different years 
years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)

ggplot(Nitrite %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Nitrite, color = factor(year))) +
  geom_line() +
  scale_colour_manual(values = c("Other" = "lightgrey",
                                 "1963" = "chocolate",
                                 "1980" = "cornflowerblue",
                                 "1994" = "darkseagreen4")) + 
  geom_vline(xintercept = 355, linetype = "dashed") + # start of winter
  geom_vline(xintercept = 80, linetype = "dashed") + # start of spring
  geom_vline(xintercept = 173, linetype = "dashed") + # start of summer
  geom_vline(xintercept = 266, linetype = "dashed") + # start of fall
  annotate("text", x = 37, y = 4, label = "Winter", size = 4) +
  annotate("text", x = 130, y = 4, label = "Spring", size = 4) +
  annotate("text", x = 220, y = 4, label = "Summer", size = 4) +
  annotate("text", x = 315, y = 4, label = "Fall", size = 4) +
  labs(x = "Day of Year", y = "Nitrite (µmol/l)", colour = "Year") +
  theme_classic()





################################################
### DIN and fit models of different families ###
################################################

# create it's own df to remove NAs
DIN <- df %>% filter(!is.na(DIN)) %>% 
  select(Date, DIN, year, month)

# plot that sucka
ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "darkgrey") +
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
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(DIN$DIN), linetype = "dashed", color = "azure4") +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()



### Different family distribution MODELS ###

DINTF2m <- gamlss(DIN ~ Date, family = TF2(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

DINSN1m <- gamlss(DIN ~ Date, family = SN1(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

DINNOm <- gamlss(DIN ~ Date, family = NO(), data = DIN, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


#====================================================================================
# Intercept model 
DINSSTm1 <- gamlss(DIN ~ 1, family = SST(), data = DIN,
                  method = mixed(10, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(DINSSTm1)

# Mean only changing through time
DINSSTm2 <- gamlss(DIN ~ Date, family = SST(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(DINSSTm2)

# Mean, sigma and nu changing through time
DINSSTm3 <- gamlss(DIN ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = DIN, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(DINSSTm3)

# Mean, sigma nu and tau changing through time
DINSSTm4 <- gamlss(DIN ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), data = DIN, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(DINSSTm4)

# mu sig and nu changing with seasonality
DINSSTm4_sea <- gamlss(DIN ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                       family = SST(), data = DIN,
                       #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(DINSSTm4_sea)
#====================================================================================

# plot that sucka
ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "grey") +
  geom_abline(intercept = DINSSTm2$mu.coefficients[1], slope = DINSSTm2$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # only mean as function of time 
  geom_abline(intercept = DINSSTm4$mu.coefficients[1], slope = DINSSTm4$mu.coefficients[2], 
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu as function of time
  geom_hline(yintercept = 23.182, linetype = "dashed", color = "black") + # mean intercept
  geom_vline(xintercept = as.Date("1963-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1980-09-23"), linetype = "dashed", color = "darkred") + 
  geom_vline(xintercept = as.Date("1994-04-13"), linetype = "dashed", color = "darkred") + 
  labs(x = "Time", y = "DIN (µmol/l)") +
  theme_minimal()

############################################ PDF

# get the pdf and param for three dates
# m1 (intercept model)
resm1 <- pdf_SST_at_date(DINSSTm1, DIN, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm1b <- pdf_SST_at_date(DINSSTm1, DIN, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm1c <- pdf_SST_at_date(DINSSTm1, DIN, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m2 (mean only)
resm2 <- pdf_SST_at_date(DINSSTm2, DIN, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm2b <- pdf_SST_at_date(DINSSTm2, DIN, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm2c <- pdf_SST_at_date(DINSSTm2, DIN, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m3 (mean sigma nu and tau)
resm3 <- pdf_SST_at_date(DINSSTm4, DIN, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm3b <- pdf_SST_at_date(DINSSTm4, DIN, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm3c <- pdf_SST_at_date(DINSSTm4, DIN, "1994-04-13",
                          time_var = "Date", date_var = "Date")


par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.075),
     ylab = "Density", xlab = "Value",
     main = "1963-09-23") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3$x, resm3$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.075),
     ylab = "Density", xlab = "Value",
     main = "1980-09-23") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3b$x, resm3b$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")
# plot the models together for date 3:1994-04-13
plot(resm1c$x, resm1c$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.075),
     ylab = "Density", xlab = "Value",
     main = "1994-04-13") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3c$x, resm3c$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(1,2,1), bty = "n")
par(mfrow = c(1, 1)) #reset

### parameters for the three dates ###
# ylim: 0.07 for all

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

# # plot all three dates 
# par(mfrow = c(1, 3))
# plot(res1$x, res1$density, type = "l", ylim=c(0,0.07)) 
# plot(res2$x, res2$density, type = "l", ylim=c(0,0.07)) # 0.04 for m1 and m2
# plot(res3$x, res3$density, type = "l", ylim=c(0,0.07)) # 0.07 for m3 and m4
# par(mfrow = c(1, 1)) # reset




#######################################################


################################################## Moment bucket 
# intercept model 
moment_bucket(DINSSTm1) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean only 
moment_bucket(DINSSTm2) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(DINSSTm3) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

###################################################################


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


########## TRENDS ################

# create day of the year for nitrate and phytopl
DIN <- DIN %>%
  mutate(DOY  = yday(Date)) 

# look at phosphate patterns for different years 
years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)

ggplot(DIN %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = DIN, color = factor(year))) +
  geom_line() +
  scale_colour_manual(values = c("Other" = "lightgrey",
                                 "1963" = "chocolate",
                                 "1980" = "cornflowerblue",
                                 "1994" = "darkseagreen4")) + 
  geom_vline(xintercept = 355, linetype = "dashed") + # start of winter
  geom_vline(xintercept = 80, linetype = "dashed") + # start of spring
  geom_vline(xintercept = 173, linetype = "dashed") + # start of summer
  geom_vline(xintercept = 266, linetype = "dashed") + # start of fall
  annotate("text", x = 37, y = 150, label = "Winter", size = 4) +
  annotate("text", x = 130, y = 150, label = "Spring", size = 4) +
  annotate("text", x = 220, y = 150, label = "Summer", size = 4) +
  annotate("text", x = 315, y = 150, label = "Fall", size = 4) +
  labs(x = "Day of Year", y = "DIN (µmol/l)", colour = "Year") +
  theme_classic()




#####################################################
### Silicate and fit models of different families ###
#####################################################

# create own df to remove NAs 
Silicate <- df %>% filter(!is.na(Silicate)) %>% 
  select(Date, Silicate, year, month)

# plot a for a specific year 
ggplot(dplyr::filter(Silicate, year(Date) == 1972),
       aes(Date, Silicate)) +
  geom_line(color = "black")

# rescale 
Silicate$y_rescaled <- Silicate$Silicate + 10

### CHECK DISTRIBUTION ###

# histogram 
hist(Silicate$Silicate)

# look at the density plot 
ggplot(Silicate, aes(y = Silicate)) +
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(Silicate$Silicate), linetype = "dashed", color = "azure4") +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### Different family distribution MODELS ###

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



# intercept model (there are 0s in the data, that SST doesn't like)

# siSSTm1 <- gamlss(y_rescaled ~ 1, family = SST(), data = Silicate, 
#                  mu.start = mean(Silicate$y_rescaled), sigma.start = max(sd(Silicate$y_rescaled), 1),
#                  method = mixed(10, 200),
#                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(siSSTm1)


# TEST
#============================================================================================
# intercept model 

### LINEAR ###

## just mean
#  ~ Date
# ~ year 
# ~ year + month 
# ~ year + month + DOY

## all param basic
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
# ~ year, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date


# intercept 
si_Intercept <- gamlss(Silicate ~ 1, family = JSU(), data = Silicate,
                    mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

## just mean
#  ~ Date
si_mean_date <- gamlss(Silicate ~ Date, family = JSU(), data = Silicate,
                    mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year 
si_mean_year <- gamlss(Silicate ~ year, family = JSU(), data = Silicate,
                    mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month
si_mean_year_month <- gamlss(Silicate ~ year + month, family = JSU(), data = Silicate,
                          mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY
si_mean_year_month_DOY <- gamlss(Silicate ~ year + month + DOY, family = JSU(), data = Silicate,
                              mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(si_mean_date)
AIC(si_mean_year)
AIC(si_mean_year_month) # best
AIC(si_mean_year_month_DOY)


## all param basic 
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
si_param_date <- gamlss(Silicate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                        family = JSU(), data = Silicate,
                     #mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
si_param_year <- gamlss(Silicate ~ year, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                        family = JSU(), data = Silicate,
                     mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
si_param_year_month <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                              family = JSU(), data = Silicate,
                           mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY
si_param_year_month_DOY <- gamlss(Silicate ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY, tau.fo = ~ year + month + DOY, 
                                  family = JSU(), data = Silicate,
                                     mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                     method = mixed(10,200),
                                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(si_param_date)
AIC(si_param_year)
AIC(si_param_year_month) # best
AIC(si_mean_year_month_DOY) 

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
si_param_date_sw_year <- gamlss(Silicate ~ Date, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                                family = JSU(), data = Silicate,
                             mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
si_param_year_month_sw_year <- gamlss(Silicate ~ year + month, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                                      family = JSU(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
si_param_year_month_sw_date <- gamlss(Silicate ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date,  tau.fo = ~ Date, 
                                      family = JSU(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
si_param_year_month_DOY_sw_year <- gamlss(Silicate ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year,  
                                          family = JSU(), data = Silicate,
                                       mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
si_param_year_month_DOY_sw_year_month <- gamlss(Silicate ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                                family = JSU(), data = Silicate,
                                             mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                             method = mixed(10,200),
                                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date
si_param_year_month_DOY_sw_date <- gamlss(Silicate ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                          family = JSU(), data = Silicate,
                                       mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


AIC(si_param_date_sw_year) 
AIC(si_param_year_month_sw_year)
AIC(si_param_year_month_sw_date)
AIC(si_param_year_month_DOY_sw_year)
AIC(si_param_year_month_DOY_sw_year_month)
AIC(si_param_year_month_DOY_sw_date)

summary(si_param_year_month) # BEST MODEL
summary(si_param_year_month_DOY_sw_year_month) # SECOND BEST MODEL
summary(si_param_year_month_sw_date) # THIRD BEST MODEL


# adding poly to best model 
poly_si_param_year_month <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                              family = JSU(), data = Silicate,
                              mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
AIC(poly_si_param_year_month) # BETTER YOOHOO
summary(poly_si_param_year_month)

# compare different families on the best model 
## too many errors
SST_poly_si_param_year_month <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ 1, 
                                   family = SST(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
AIC(SST_poly_si_param_year_month) # does not do well
#==========================================================================================


#==========================================================================================
# Intercept model with SEP3 => mu is mode, not mean
# siSEPm1 <- gamlss(Silicate ~ 1, family = SEP3(), data = Silicate, 
#                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
#                  method = mixed(20, 200),
#                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(siSEPm1)
# mfv(Silicate$Silicate)

# Intercept model
siJSUm1 <- gamlss(Silicate ~ 1, family = JSU(), data = Silicate, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(siJSUm1)

# Mean only 
siJSUm2 <- gamlss(Silicate ~ Date, family = JSU(), data = Silicate, 
                  method = mixed(10, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(siJSUm2)

# Mean sigma and nu changing through time
siJSUm3 <- gamlss(Silicate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = JSU(), data = Silicate, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(siJSUm3)

# mu sig and nu changing with seasonality
siJSUm3_sea <- gamlss(Silicate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, 
                       family = JSU(), data = Silicate,
                       #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(siJSUm3_sea)
#==========================================================================================
# predict mu from seasonailty model
Silicate$mu_hat <- fitted(si_param_year_month, what = "mu")
Silicate$mu_hat2 <- fitted(poly_si_param_year_month, what = "mu")

# plot that sucka 
ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "grey") +
  geom_abline(intercept = si_mean_date$mu.coefficients[1], slope = si_mean_date$mu.coefficients[2],
              color = "steelblue", linewidth = 1) + # mean only
  geom_abline(intercept = si_param_date$mu.coefficients[1], slope = si_param_date$mu.coefficients[2],
              color = "goldenrod", linewidth = 1) + # mean sigma and nu changing through time
  geom_line(aes(y = mu_hat2), color = "darkolivegreen", linewidth = 1) +
  geom_hline(yintercept = 7.378, linetype = "dashed", color = "black") +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()

# plot that sucka, with seasonaility 
Silicate$mu_hat <- fitted(siJSUm3_sea, what = "mu")
ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "azure4", na.rm = TRUE) +
  geom_line(aes(y = mu_hat), color = "darkred", linewidth = 1) +
  theme_minimal()

############################################ PDF

# get the pdf and param for three dates
# m1 (intercept model)
resm1 <- pdf_SST_at_date(siJSUm1, Silicate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm1b <- pdf_SST_at_date(siJSUm1, Silicate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm1c <- pdf_SST_at_date(siJSUm1, Silicate, "1992-09-23",
                          time_var = "Date", date_var = "Date")
# m2 (mean only)
resm2 <- pdf_SST_at_date(siJSUm2, Silicate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm2b <- pdf_SST_at_date(siJSUm2, Silicate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm2c <- pdf_SST_at_date(siJSUm2, Silicate, "1992-09-23",
                          time_var = "Date", date_var = "Date")
# m3 (mean sigma nu and tau)
resm3 <- pdf_SST_at_date(siJSUm3, Silicate, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm3b <- pdf_SST_at_date(siJSUm3, Silicate, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm3c <- pdf_SST_at_date(siJSUm3, Silicate, "1992-09-23",
                          time_var = "Date", date_var = "Date")


par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "goldenrod", type = "l", lwd = 2, ylim=c(0,0.075),
     ylab = "Density", xlab = "Value",
     main = "1963") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 2) # mean only model
lines(resm3$x, resm3$density, col = "grey",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("goldenrod", "steelblue", "grey"), cex = 0.7, lty = c(1,2,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "goldenrod", type = "l", lwd = 2, ylim=c(0,0.075),
     ylab = "Density", xlab = "Value",
     main = "1980") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 2) # mean only model
lines(resm3b$x, resm3b$density, col = "grey",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("goldenrod", "steelblue", "grey"), cex = 0.7, lty = c(1,2,1),  bty = "n")
# plot the models together for date 3: 1992-09-23
plot(resm1c$x, resm1c$density, col = "goldenrod", type = "l", lwd = 2, ylim=c(0,0.075),
     ylab = "Density", xlab = "Value",
     main = "1992") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 2) # mean only model
lines(resm3c$x, resm3c$density, col = "grey",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("goldenrod", "steelblue", "grey"), cex = 0.7, lty = c(1,2,1),  bty = "n")
par(mfrow = c(1, 1)) #reset

### parameters for the three dates ###
# ylim: 0.07 for all

# m1 (intercept model)
resm1$params #1963
resm1b$params #1980
resm1c$params #1992
# m2 (mean only)
resm2$params #1963
resm2b$params #1980
resm2c$params #1992
# m3 (mean sigma and nu)
resm3$params #1963
resm3b$params #1980
resm3c$params #1992


# # plot all three dates 
# par(mfrow = c(1, 3))
# plot(res1$x, res1$density, type = "l", ylim=c(0,0.25)) # 0.13 for m1 
# plot(res2$x, res2$density, type = "l", ylim=c(0,0.25)) # 0.17 for m2
# plot(res3$x, res3$density, type = "l", ylim=c(0,0.25)) # 0.25 for m3
# par(mfrow = c(1, 1)) # reset


#######################################################
###### Look at moment changes every n years ########


summary <- Silicate %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Silicate, na.rm = TRUE),
    var = var(Silicate, na.rm = TRUE),
    skew = skewness(Silicate, na.rm = TRUE),
    kurt = kurtosis(Silicate, na.rm = TRUE)
  )

print(summary, n=33)

###### Look at moment changes average months ########

month_moments <- Silicate %>%
  mutate(
    month_num = month(Date),
    month_lab = month(Date, label = TRUE, abbr = TRUE)
  ) %>%
  group_by(month_num, month_lab) %>%
  summarise(
    n     = sum(!is.na(Silicate)),
    mean  = mean(Silicate, na.rm = TRUE),
    var   = var(Silicate, na.rm = TRUE),
    skew  = skewness(Silicate, na.rm = TRUE),
    kurt  = kurtosis(Silicate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_num)

month_moments

# compare the param with the predicted param of the best model 
poly_si_param_year_month <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                   family = JSU(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

poly_si_param_year_month_meanonly <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1, 
                                   family = JSU(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

poly_si_param_year_month_noskew <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ 1, tau.fo = ~ year + month, 
                                   family = JSU(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

poly_si_param_year_month_nokurt <- gamlss(Silicate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ 1, 
                                   family = JSU(), data = Silicate,
                                   mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


AIC(poly_si_param_year_month)
AIC(poly_si_param_year_month_meanonly)
AIC(poly_si_param_year_month_noskew)
AIC(poly_si_param_year_month_nokurt)


# create function to get predicted param from model
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- fitted(model, what = "mu")
  data$sigma_hat <- fitted(model, what = "sigma")
  data$nu_hat <- fitted(model, what = "nu")
  data$tau_hat <- fitted(model, what = "tau")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
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

model_pred_param <- best_model_param(poly_si_param_year_month, Silicate)
model_pred_param2 <- best_model_param(poly_si_param_year_month_meanonly, Silicate)
model_pred_param3 <- best_model_param(poly_si_param_year_month_noskew, Silicate)
model_pred_param4 <- best_model_param(poly_si_param_year_month_nokurt, Silicate)
print(model_pred_param, n=33)
print(model_pred_param)


## PLOT comparison per year 

# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(2,18), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "mean")
lines(model_pred_param2$mean_mu_hat ~ model_pred_param2$year, col = "steelblue", lwd = 2) # mean only 
lines(model_pred_param3$mean_mu_hat ~ model_pred_param3$year, col = "chocolate", lwd = 2) # constant skewness
lines(model_pred_param4$mean_mu_hat ~ model_pred_param4$year, col = "darkolivegreen", lwd = 2) # constant kurtosis
lines(summary$mean ~ summary$year, col = "azure4", lwd = 2, lty = 2)
legend("topleft", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.6, bty = "n")


# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,68), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "variance")
lines(model_pred_param2$mean_sigma_hat ~ model_pred_param2$year, col = "steelblue", lwd = 2) # mean only 
lines(model_pred_param3$mean_sigma_hat ~ model_pred_param3$year, col = "chocolate", lwd = 2) # constant skewness
lines(model_pred_param4$mean_sigma_hat ~ model_pred_param4$year, col = "darkolivegreen", lwd = 2) # constant kurtosis
lines(summary$var ~ summary$year, col = "azure4", lwd = 2, lty = 2)
legend("topleft", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.6, bty = "n")


# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(-0.3,11), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "skewness")
lines(model_pred_param2$mean_nu_hat ~ model_pred_param2$year, col = "steelblue", lwd = 2) # mean only 
lines(model_pred_param3$mean_nu_hat ~ model_pred_param3$year, col = "chocolate", lwd = 2) # constant skewness
lines(model_pred_param4$mean_nu_hat ~ model_pred_param4$year, col = "darkolivegreen", lwd = 2) # constant kurtosis
lines(summary$skew ~ summary$year, col = "azure4", lwd = 2, lty = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.5, bty = "n")
abline(h = 0, lty = 3, col = "red") # symmetrical for empirical


# plot empirical tau vs predicted tau
plot(model_pred_param$mean_tau_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,9), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "kurtosis")
lines(model_pred_param2$mean_tau_hat ~ model_pred_param2$year, col = "steelblue", lwd = 2) # mean only 
lines(model_pred_param3$mean_tau_hat ~ model_pred_param3$year, col = "chocolate", lwd = 2) # constant skewness
lines(model_pred_param4$mean_tau_hat ~ model_pred_param4$year, col = "darkolivegreen", lwd = 2) # constant kurtosis
lines(summary$kurt ~ summary$year, col = "azure4", lwd = 2, lty = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.5, bty = "n")
abline(h = 3, lty = 3, col = "red") # normal tails for empirical



## PLOT comparison per month 
# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$month, type = "l", 
     ylim = c(2,12), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "mean")
lines(model_pred_param2$mean_tau_hat ~ model_pred_param2$month, col = "steelblue", lwd = 2) # mean only 
lines(model_pred_param3$mean_tau_hat ~ model_pred_param3$month, col = "chocolate", lwd = 2) # constant skewness
lines(model_pred_param4$mean_tau_hat ~ model_pred_param4$month, col = "darkolivegreen", lwd = 2) # constant kurtosis
lines(month_moments$mean ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0,67), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "variance")
lines(month_moments$var ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$month, type = "l", 
     ylim = c(0, 20), col = "goldenrod", lwd = 1,
     xlab = "month", ylab = "skew")
lines(month_moments$skew ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2)

# plot empirical tau vs predicted tau
plot(model_pred_param$mean_tau_hat ~ model_pred_param$month, type = "l",
ylim = c(1, 20), col = "goldenrod", lwd = 1,
xlab = "month", ylab = "kurtosis")
lines(month_moments$kurt ~ month_moments$month_lab, col = "azure4", lwd = 2, lty = 2, type = "l")
abline(h = 3, lty = 3, col = "red") # normal tails


################################################## Moment bucket 
# intercept model 
moment_bucket(siJSUm1) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean only 
moment_bucket(siJSUm2) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(siJSUm3) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

###################################################################


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


########## TRENDS ###############

# create day of the year for nitrate and phytopl
Silicate <- Silicate %>%
  mutate(DOY  = yday(Date)) 

# look at phosphate patterns for different years 
years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)

ggplot(Silicate %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Silicate, color = factor(year))) +
  geom_line() +
  scale_colour_manual(values = c("Other" = "lightgrey",
                                 "1963" = "chocolate",
                                 "1980" = "cornflowerblue",
                                 "1994" = "darkseagreen4")) + 
  geom_vline(xintercept = 355, linetype = "dashed") + # start of winter
  geom_vline(xintercept = 80, linetype = "dashed") + # start of spring
  geom_vline(xintercept = 173, linetype = "dashed") + # start of summer
  geom_vline(xintercept = 266, linetype = "dashed") + # start of fall
  annotate("text", x = 37, y = 40, label = "Winter", size = 4) +
  annotate("text", x = 130, y = 40, label = "Spring", size = 4) +
  annotate("text", x = 220, y = 40, label = "Summer", size = 4) +
  annotate("text", x = 315, y = 40, label = "Fall", size = 4) +
  labs(x = "Day of Year", y = "Silicate (µmol/l)", colour = "Year") +
  theme_classic()


#####################################################
### Ammonium and fit models of different families ###
#####################################################

# create own df to remove NAs 
Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium, year, month)

# plot that sucka
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "darkgrey") +
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
  geom_density(linewidth = 0.8) +
  geom_hline(yintercept = mean(Ammonium$Ammonium), linetype = "dashed", color = "azure4") +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


### Different family distribution MODELS ###

amTF2m <- gamlss(Ammonium ~ Date, family = TF2(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

amSN1m <- gamlss(Ammonium ~ Date, family = SN1(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

amNOm <- gamlss(Ammonium ~ Date, family = NO(), data = Ammonium, 
                method = mixed(5, 200),
                control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


# TEST
#============================================================================================
# intercept model 

### LINEAR ###

## just mean
#  ~ Date
# ~ year 
# ~ year + month 
# ~ year + month + DOY

## all param basic
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
# ~ year, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date


# intercept 
am_Intercept <- gamlss(Ammonium ~ 1, family = SST(), data = Ammonium,
                    method = mixed(5,200),
                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

## just mean
#  ~ Date
am_mean_date <- gamlss(Ammonium ~ Date, family = SST(), data = Ammonium,
                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year 
am_mean_year <- gamlss(Ammonium ~ year, family = SST(), data = Ammonium,
                    mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                    method = mixed(10,200),
                    control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month
am_mean_year_month <- gamlss(Ammonium ~ year + month, family = SST(), data = Ammonium,
                          mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY
am_mean_year_month_DOY <- gamlss(Ammonium ~ year + month + DOY, family = SST(), data = Ammonium,
                              mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(am_mean_date)
AIC(am_mean_year)
AIC(am_mean_year_month) # best
AIC(am_mean_year_month_DOY)


## all param basic 
#  ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date 
am_param_date <- gamlss(Ammonium ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date,
                        family = SST(), data = Ammonium,
                     #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
am_param_year <- gamlss(Ammonium ~ year, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                        family = SST(), data = Ammonium,
                     #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                     method = mixed(10,200),
                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month
am_param_year_month <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                              family = SST(), data = Ammonium,
                           mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                           method = mixed(5,200),
                           control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY
am_param_year_month_DOY <- gamlss(Ammonium ~ year + month + DOY, sigma.fo = ~ year + month + DOY, nu.fo = ~ year + month + DOY, tau.fo = ~ year, 
                                  family = SST(), data = Ammonium,
                                     mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                     method = mixed(10,200),
                                     control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(am_param_date)
AIC(am_param_year)
AIC(am_param_year_month) # worked for JSU; for SST, tau has to be ~ year only
AIC(am_mean_year_month_DOY) # tau has to be ~ year only 

# all param switch-a-roo
#  ~ Date, sigma.fo = ~ year, nu.fo = ~ year 
am_param_date_sw_year <- gamlss(Ammonium ~ Date, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                                family = SST(), data = Ammonium,
                             mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ year, nu.fo = ~ year
am_param_year_month_sw_year <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                                      family = SST(), data = Ammonium,
                                   #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date
am_param_year_month_sw_date <- gamlss(Ammonium ~ year + month, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                      family = SST(), data = Ammonium,
                                   #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(10,200),
                                   control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year 
am_param_year_month_DOY_sw_year <- gamlss(Ammonium ~ year + month + DOY, sigma.fo = ~ year, nu.fo = ~ year, tau.fo = ~ year, 
                                          family = SST(), data = Ammonium,
                                       #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month
am_param_year_month_DOY_sw_year_month <- gamlss(Ammonium ~ year + month + DOY, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                                family = SST(), data = Ammonium,
                                             mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                             method = mixed(5,200),
                                             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date
am_param_year_month_DOY_sw_date <- gamlss(Ammonium ~ year + month + DOY, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, 
                                          family = SST(), data = Ammonium,
                                       #mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                       method = mixed(10,200),
                                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


AIC(am_param_date_sw_year) 
AIC(am_param_year_month_sw_year)
AIC(am_param_year_month_sw_date)
AIC(am_param_year_month_DOY_sw_year)
AIC(am_param_year_month_DOY_sw_year_month) # tau had to be ~ year only
AIC(am_param_year_month_DOY_sw_date)

summary(am_param_year_month) # BEST MODEL
summary(am_param_year_month_DOY_sw_year_month) # SECOND BEST MODEL
summary(am_param_year_month_sw_year) # THIRD BEST MODEL

# adding poly to best model 
poly_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                              family = SST(), data = Ammonium,
                              mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                              method = mixed(5,200),
                              control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
AIC(poly_am_param_year_month) # BETTER YOOHOO
summary(poly_am_param_year_month)

# compare different families on the best model
JSU_poly_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = JSU(), data = Ammonium,
                                   mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

SHASHo_poly_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = SHASHo(), data = Ammonium,
                                   mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

SEP3_poly_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = SEP3(), data = Ammonium,
                                   mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

NET_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = NET(), data = Ammonium,
                                   mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

TF2_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, tau.fo = ~ year, 
                                  family = TF2(), data = Ammonium,
                                  mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                  method = mixed(5,200),
                                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

AIC(JSU_poly_am_param_year_month) # tis is better
AIC(SHASHo_poly_am_param_year_month)
AIC(SEP3_poly_am_param_year_month)
AIC(NET_am_param_year_month)
AIC(TF2_am_param_year_month) # not good



# make sigma or skew or tau constant at a time 
poly_am_param_year_month_meanonly <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1, 
                                          family = JSU(), data = Ammonium,
                                          mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                          method = mixed(5,200),
                                          control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

poly_am_param_year_month_noskew <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ 1, tau.fo = ~ year, 
                                   family = JSU(), data = Ammonium,
                                   mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


poly_am_param_year_month_nokurt <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ 1, 
                                          family = JSU(), data = Ammonium,
                                          mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                          method = mixed(5,200),
                                          control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

AIC(JSU_poly_am_param_year_month)
AIC(poly_am_param_year_month_meanonly) # holy moly its so bad
AIC(poly_am_param_year_month_noskew) # does even worse
AIC(poly_am_param_year_month_nokurt) # does worse 

#=========================================================================================

#=========================================================================================
# Intercept model 
amSSTm1 <- gamlss(Ammonium ~ 1, family = SST(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(amSSTm1)

# mean as a function of time 
amSSTm2 <- gamlss(Ammonium ~ Date, family = SST(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(amSSTm2)

# mean, sigma and nu as a function of time 
amSSTm3 <- gamlss(Ammonium ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = Ammonium, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(amSSTm3)

# mean, sigma nu and tau as a function of time 
amSSTm4 <- gamlss(Ammonium ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), data = Ammonium, 
                  method = mixed(5, 200),
                  control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(amSSTm4)

# mu sig and nu changing with seasonality
amSSTm4_sea <- gamlss(Ammonium ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = SST(), data = Ammonium,
                      #mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
summary(amSSTm4_sea)
#=========================================================================================
# predict mu from seasonality model 
Ammonium$mu_hat <- fitted(am_param_year_month, what = "mu")
Ammonium$mu_hat2 <- fitted(poly_am_param_year_month, what = "mu")

# plot that sucka
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey") +
  geom_abline(intercept = am_mean_date$mu.coefficients[1], slope = am_mean_date$mu.coefficients[2], 
              color = "steelblue", linewidth = 1) + # only mean as function of time 
  geom_abline(intercept = am_param_date$mu.coefficients[1], slope = am_param_date$mu.coefficients[2], 
              color = "goldenrod", linewidth = 1) + # mean, sigma and nu as function of time
  geom_line(aes(y = mu_hat), color = "darkolivegreen", linewidth = 1) +
  geom_vline(xintercept = as.Date("1963-09-23"), linetype = "dashed" ,color = "darkred") +
  geom_vline(xintercept = as.Date("1980-09-23"), linetype = "dashed" ,color = "darkred") +
  geom_vline(xintercept = as.Date("1994-04-13"), linetype = "dashed" ,color = "darkred") +
  geom_hline(yintercept = 5.8704, linetype = "dashed", color = "black") +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()

Ammonium$mu_hat2 <- fitted(JSU_poly_am_param_year_month, what = "mu")
Ammonium$mu_hat3 <- fitted(poly_am_param_year_month_meanonly, what = "mu")
Ammonium$mu_hat4 <- fitted(poly_am_param_year_month_noskew, what = "mu")
Ammonium$mu_hat5 <- fitted(poly_am_param_year_month_nokurt, what = "mu")

# plot all seasonaility model 
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "grey") +
  geom_line(aes(y = mu_hat2), color = "goldenrod", linewidth = 1) + # all param
  geom_line(aes(y = mu_hat3), color = "steelblue", linewidth = 1) + # mean only
  geom_line(aes(y = mu_hat4), color = "chocolate", linewidth = 1) + # constant skew
  geom_line(aes(y = mu_hat5), color = "darkolivegreen", linewidth = 1) + # constant kurtosis
  geom_hline(yintercept = 5.8704, linetype = "dashed", color = "black") +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_classic()

############################################ PDF

# get the pdf and param for three dates
# m1 (intercept model)
resm1 <- pdf_SST_at_date(amSSTm1, Ammonium, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm1b <- pdf_SST_at_date(amSSTm1, Ammonium, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm1c <- pdf_SST_at_date(amSSTm1, Ammonium, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m2 (mean only)
resm2 <- pdf_SST_at_date(amSSTm2, Ammonium, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm2b <- pdf_SST_at_date(amSSTm2, Ammonium, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm2c <- pdf_SST_at_date(amSSTm2, Ammonium, "1994-04-13",
                          time_var = "Date", date_var = "Date")
# m3 (mean sigma nu and tau)
resm3 <- pdf_SST_at_date(amSSTm4, Ammonium, "1963-09-23",
                         time_var = "Date", date_var = "Date")
resm3b <- pdf_SST_at_date(amSSTm4, Ammonium, "1980-09-23",
                          time_var = "Date", date_var = "Date")
resm3c <- pdf_SST_at_date(amSSTm4, Ammonium, "1994-04-13",
                          time_var = "Date", date_var = "Date")


par(mfrow = c(1, 3)) 
# plot the models together for date 1: 1963-09-23
plot(resm1$x, resm1$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.3),
     ylab = "Density", xlab = "Value",
     main = "1963-09-23") # intercept model
lines(resm2$x, resm2$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3$x, resm3$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "grey"), cex = 0.7, lty = c(2,1,1), bty = "n")
# plot the models together for date 2: 1980-09-23
plot(resm1b$x, resm1b$density, col = "black", type = "l", lwd = 2, lty = 2, ylim=c(0,0.3),
     ylab = "Density", xlab = "Value",
     main = "1980-09-23") # intercept model
lines(resm2b$x, resm2b$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3b$x, resm3b$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")
# plot the models together for date 3: 1994-04-13
plot(resm1c$x, resm1c$density, col = "black", type = "l", lwd = 2, lty = 2,  ylim=c(0,0.1), xlim =c(5,15),
     ylab = "Density", xlab = "Value",
     main = "1994-04-13") # intercept model
lines(resm2c$x, resm2c$density, col = "steelblue", lwd = 2, lty = 1) # mean only model
lines(resm3c$x, resm3c$density, col = "goldenrod",lwd = 2, lty = 1) # mean sigma and nu 
legend("topright", legend = c("Intercept model", "mean(time) only model", "all param(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.7, lty = c(2,1,1), bty = "n")
par(mfrow = c(1, 1)) #reset

### parameters for the three dates ###
# ylim: 0.3 for all

# m1 (intercept model)
resm1$params #1963
resm1b$params #1980
resm1c$params #1994
# m2 (mean only)
resm2$params #1963
resm2b$params #1980
resm2c$params #1994
# m3 (mean sigma nu and tau)
resm3$params #1963
resm3b$params #1980
resm3c$params #1994


# # plot all three dates 
# par(mfrow = c(1, 3))
# plot(res1$x, res1$density, type = "l", ylim=c(0,0.25)) # 0.13 for m1 
# plot(res2$x, res2$density, type = "l", ylim=c(0,0.25)) # 0.17 for m2
# plot(res3$x, res3$density, type = "l", ylim=c(0,0.25)) # 0.25 for m3
# par(mfrow = c(1, 1)) # reset


#######################################################

#######################################################
###### Look at moment changes every n years ########


summary <- Ammonium %>% 
  mutate(year = as.numeric(format(Date, "%Y"))) %>% 
  group_by(year = floor(year/1)*1) %>% 
  summarise(
    mean = mean(Ammonium, na.rm = TRUE),
    var = var(Ammonium, na.rm = TRUE),
    skew = skewness(Ammonium, na.rm = TRUE),
    kurt = kurtosis(Ammonium, na.rm = TRUE)
  )

print(summary, n=33)

###### Look at moment changes average months ########

month_moments <- Ammonium %>%
  mutate(
    month_num = month(Date),
    month_lab = month(Date, label = TRUE, abbr = TRUE)
  ) %>%
  group_by(month_num, month_lab) %>%
  summarise(
    n     = sum(!is.na(Ammonium)),
    mean  = mean(Ammonium, na.rm = TRUE),
    var   = var(Ammonium, na.rm = TRUE),
    skew  = skewness(Ammonium, na.rm = TRUE),
    kurt  = kurtosis(Ammonium, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_num)

month_moments

# compare the param with the predicted param of the best model 
poly_am_param_year_month <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year, 
                                   family = SST(), data = Ammonium,
                                   mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                   method = mixed(5,200),
                                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

# create function to get predicted param from model
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- fitted(model, what = "mu")
  data$sigma_hat <- fitted(model, what = "sigma")
  data$nu_hat <- fitted(model, what = "nu")
  data$tau_hat <- fitted(model, what = "tau")
  
  # get the average per year 
  average_by_year <- data %>%
    group_by(year) %>%
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

model_pred_param <- best_model_param(JSU_poly_am_param_year_month, Ammonium)
print(model_pred_param, n=33)
print(model_pred_param)
model_pred_param2 <- best_model_param(poly_am_param_year_month_meanonly, Ammonium)
model_pred_param3 <- best_model_param(poly_am_param_year_month_noskew, Ammonium)
model_pred_param4 <- best_model_param(poly_am_param_year_month_nokurt, Ammonium)


## PLOT comparison per year 

# plot empirical mean vs predicted mean
plot(model_pred_param$mean_mu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,13), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "mean")
lines(model_pred_param2$mean_mu_hat ~ model_pred_param2$year, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_mu_hat ~ model_pred_param3$year, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_mu_hat ~ model_pred_param4$year, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(summary$mean ~ summary$year, col = "azure4", lwd = 2, lty = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
        col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.6, bty = "n")

# plot empirical sigma vs predicted sigma
plot(model_pred_param$mean_sigma_hat ~ model_pred_param$year, type = "l", 
     ylim = c(0,22), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "variance")
lines(model_pred_param2$mean_sigma_hat ~ model_pred_param2$year, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_sigma_hat ~ model_pred_param3$year, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_sigma_hat ~ model_pred_param4$year, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(summary$var ~ summary$year, col = "azure4", lwd = 2, lty = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.6, bty = "n")


# plot empirical skew vs predicted skew
plot(model_pred_param$mean_nu_hat ~ model_pred_param$year, type = "l", 
     ylim = c(-0.1, 6), col = "goldenrod", lwd = 2,
     xlab = "year", ylab = "skewness")
lines(model_pred_param2$mean_nu_hat ~ model_pred_param2$year, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_nu_hat ~ model_pred_param3$year, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_nu_hat ~ model_pred_param4$year, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(summary$skew ~ summary$year, col = "azure4", lwd = 2, lty = 2)
abline(h = 0, lty = 3, col = "red") # symmetrical for empirical
abline(h = 1, lty = 2, col = "darkred") # symmetrical for SST
legend("topleft", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.6, bty = "n")


# plot empirical tau vs predicted tau
plot(model_pred_param$mean_tau_hat ~ model_pred_param$year, type = "l",
     ylim = c(1, 20), col = "goldenrod", lwd = 2,
    xlab = "year", ylab = "kurtosis")
lines(model_pred_param2$mean_tau_hat ~ model_pred_param2$year, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_tau_hat ~ model_pred_param3$year, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_tau_hat ~ model_pred_param4$year, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(summary$kurt ~ summary$year, col = "azure4", lwd = 2, type = "l", lty = 2)
abline(h = 3, lty = 3, col = "red") # normal tails
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "azure4" ),  lty = c(1,1,1,1,2), cex = 0.6, bty = "n")



## PLOT comparison per month 
# plot empirical mean vs predicted mean
plot(month_moments$mean ~ month_moments$month_lab, type = "l", 
     ylim = c(4,8), col = "azure4", lwd = 1,
     xlab = "month", ylab = "mean")
lines(model_pred_param2$mean_mu_hat ~ model_pred_param2$month, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_mu_hat ~ model_pred_param3$month, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_mu_hat ~ model_pred_param4$month, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(model_pred_param$mean_mu_hat ~ model_pred_param$month, col = "goldenrod", lwd = 2)
legend("topleft", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "black" ),  lty = c(1,1,1,1,1), cex = 0.6, bty = "n")

# plot empirical sigma vs predicted sigma
plot(month_moments$var ~ month_moments$month_lab, type = "l", 
     ylim = c(0,16), col = "azure4", lwd = 1,
     xlab = "month", ylab = "variance")
lines(model_pred_param2$mean_sigma_hat ~ model_pred_param2$month, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_sigma_hat ~ model_pred_param3$month, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_sigma_hat ~ model_pred_param4$month, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(model_pred_param$mean_sigma_hat ~ model_pred_param$month, col = "goldenrod", lwd = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "black" ),  lty = c(1,1,1,1,1), cex = 0.4, bty = "n")

# plot empirical skew vs predicted skew
plot(month_moments$skew ~ month_moments$month_lab, type = "l", 
     ylim = c(-0.1,15), col = "azure4", lwd = 1,
     xlab = "month", ylab = "skewness")
lines(model_pred_param2$mean_nu_hat ~ model_pred_param2$month, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_nu_hat ~ model_pred_param3$month, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_nu_hat ~ model_pred_param4$month, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(model_pred_param$mean_nu_hat ~ model_pred_param$month, col = "goldenrod", lwd = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "black" ),  lty = c(1,1,1,1,1), cex = 0.6, bty = "n")
abline(h = 0, lty = 3, col = "red") # symmetrical for empirical
abline(h = 1, lty = 2, col = "darkred") # symmetrical for SST

# plot empirical tau vs predicted tau
plot(month_moments$kurt ~ month_moments$month_lab, type = "l", 
     ylim = c(1,15), col = "azure4", lwd = 1,
     xlab = "month", ylab = "kurtosis")
lines(model_pred_param2$mean_tau_hat ~ model_pred_param2$month, type = "l", col = "steelblue", lwd = 2) # mean only
lines(model_pred_param3$mean_tau_hat ~ model_pred_param3$month, type = "l", col = "chocolate", lwd = 2) # no skew
lines(model_pred_param4$mean_tau_hat ~ model_pred_param4$month, type = "l", col = "darkolivegreen", lwd = 2) # no kurt
lines(model_pred_param$mean_tau_hat ~ model_pred_param$month, col = "goldenrod", lwd = 2)
legend("topright", legend = c("all param(time)", "mean(time)", "constant skew", "constant kurtosis", "empirical estimate"), 
       col = c("goldenrod", "steelblue", "chocolate", "darkolivegreen", "black" ),  lty = c(1,1,1,1,1), cex = 0.6, bty = "n")
abline(h = 3, lty = 3, col = "red") # normal tails

################################################## Worm plot 

# JSU FAMILY
resid_wp(JSU_poly_am_param_year_month) # slight s-shape = kurtosis not fitted properly
resid_wp(poly_am_param_year_month_meanonly) # s-shape = kurtosis not fitted properly
resid_wp(poly_am_param_year_month_noskew) # slight s-shape = kurtosis not fitted properly
resid_wp(poly_am_param_year_month_nokurt) # s-shape = kurtosis not fitted properly

#SST FAMILY 
resid_wp(poly_am_param_year_month) # slight s-shape = kurtosis not fitted properly
resid_wp(poly_am_param_year_month_meanonly) # s-shape = kurtosis not fitted properly
resid_wp(poly_am_param_year_month_noskew) # slight s-shape = kurtosis not fitted properly
resid_wp(poly_am_param_year_month_nokurt) 



################################################## Moment bucket 
library(gamlss.ggplots)

moment_bucket(JSU_poly_am_param_year_month, poly_am_param_year_month_meanonly, 
              poly_am_param_year_month_noskew,poly_am_param_year_month_nokurt, cex_text = 2) +
  theme_bw() + 
  theme(legend.position = "none") 

# JSU vs SST
moment_bucket(JSU_poly_am_param_year_month, poly_am_param_year_month, cex_text = 2) +
  theme_bw() + 
  theme(legend.position = "none") 

# intercept model 
moment_bucket(amSSTm1) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean only 
moment_bucket(poly_am_param_year_month_nokurt) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

# mean, sigma and nu changing through time
moment_bucket(JSU_poly_am_param_year_month) + 
  theme_bw() + 
  ggtitle("(c)") +
  theme(legend.position = "none")

###################################################################



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



########## TRENDS ################

# create day of the year for nitrate and phytopl
Ammonium <- Ammonium %>%
  mutate(DOY  = yday(Date)) 

# look at phosphate patterns for different years 
years_to_plot <- c(1963, 1975, 1980, 1966, 1994, 1962, 1978, 1992)
years_to_plot <- c(1963, 1980, 1994)

ggplot(Ammonium %>% filter(year %in% years_to_plot),
       aes(x = DOY, y = Ammonium, color = factor(year))) +
  geom_line() +
  scale_colour_manual(values = c("Other" = "lightgrey",
                                 "1963" = "chocolate",
                                 "1980" = "cornflowerblue",
                                 "1994" = "darkseagreen4")) + 
  geom_vline(xintercept = 355, linetype = "dashed") + # start of winter
  geom_vline(xintercept = 80, linetype = "dashed") + # start of spring
  geom_vline(xintercept = 173, linetype = "dashed") + # start of summer
  geom_vline(xintercept = 266, linetype = "dashed") + # start of fall
  annotate("text", x = 37, y = 22, label = "Winter", size = 4) +
  annotate("text", x = 130, y = 22, label = "Spring", size = 4) +
  annotate("text", x = 220, y = 22, label = "Summer", size = 4) +
  annotate("text", x = 315, y = 22, label = "Fall", size = 4) +
  labs(x = "Day of Year", y = "Ammonium (µmol/l)", colour = "Year") +
  theme_classic()
