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



# subset all nutrients 
Nitrate <- df %>% filter(!is.na(Nitrate)) %>% 
  select(Date, Nitrate, year, month)

Phosphate <- df %>% filter(!is.na(Phosphate)) %>% 
  select(Date, Phosphate, year, month)

Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite, year, month)

Silicate <- df %>% filter(!is.na(Silicate)) %>% 
  select(Date, Silicate, year, month)

DIN <- df %>% filter(!is.na(DIN)) %>% 
  select(Date, DIN, year, month)

Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium, year, month)



# instead of an intercept model, now we add time-varying parameters, using the family distribution 
# that does the best forthe intercept model foe each nutrients 


##########################
### TIME-VARYING MODEL ###
##########################

#### NITRATE #### ========================================================================================

## STEP 1: ----------
# compare models of the same family, with different time-varying or constant parameters

ni_meanonly_SHASHo <- gamlss(Nitrate ~ year + month, family = SHASHo(), data = Nitrate,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


ni_mean_and_sigma_SHASHo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, family = SHASHo(), data = Nitrate,
                                 method = mixed(10,200),
                                 control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


ni_noskew_SHASHo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, tau.fo = ~ year + month, family = SHASHo(), data = Nitrate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


ni_nokurt_SHASHo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, family = SHASHo(), data = Nitrate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


ni_all_SHASHo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = SHASHo(), data = Nitrate,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

ni_poly_all_SHASHo <- gamlss(Nitrate ~ poly(year,2) + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = SHASHo(), data = Nitrate,
                        method = mixed(5,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


ni_all_SHASHotr <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = SHASHotr(), data = Nitrate,
                        mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),nu.start = 2, tau.start = 1,
                        method = mixed(10,100),
                        control = gamlss.control(n.cyc = 100, c.crit = 0.01, trace = TRUE))


# AIC checks 
AIC(ni_meanonly_SHASHo)
AIC(ni_mean_and_sigma_SHASHo)
AIC(ni_noskew_SHASHo)
AIC(ni_nokurt_SHASHo)
AIC(ni_all_SHASHo) # tis the best
AIC(ni_poly_all_SHASHo) # does not do well when added poly?



Nitrate$mu_hat_best <- predict(ni_all_SHASHo, what = "mu", type = "response")
Nitrate$mu_hat_meanonly <- predict(ni_meanonly_SHASHo, what = "mu", type = "response")
Nitrate$mu_hat_BCT <- predict(ni_all_BCTo, what = "mu", type = "response")


# plot at see real quick 
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "grey") +
  #geom_line(aes(y = mu_hat_best), color = "darkolivegreen", linewidth = 1) + 
  #geom_line(aes(y = mu_hat_meanonly), color = "steelblue", linewidth = 1) + 
  geom_line(aes(y = mu_hat_BCT), color = "goldenrod", linewidth = 1) + 
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()




ni_all_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                        family = BCTo(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

AIC(ni_all_BCTo)

summary(ni_all_BCTo)


# the model that does the best is the one where all the parameters are varying though time with 
# seasonality (year + month)


# get predictive value for BCT and SHASHo model 
Nitrate$SHASHo_mu_hat <- predict(ni_all_SHASHo, what = "mu", type = "response")
Nitrate$BCTo_mu_hat <- predict(ni_all_BCTo, what = "mu", type = "response")



# pit histogram for SHASHo
# Get fitted parameters for each observation:
mu_hat_SHASHo    <- predict(ni_all_SHASHo, "mu", type = "response")
sigma_hat_SHASHo <- predict(ni_all_SHASHo, "sigma", type = "response")
nu_hat_SHASHo    <- predict(ni_all_SHASHo, "nu", type = "response")
tau_hat_SHASHo   <- predict(ni_all_SHASHo, "tau", type = "response")

# Compute PIT values (using pSEP4 for the SEP4 family, replace with appropriate p* function for your family):
pit2 <- pSHASHo(Nitrate$Nitrate, mu = mu_hat_SHASHo, sigma = sigma_hat_SHASHo, nu = nu_hat_SHASHo, tau = tau_hat_SHASHo)

# Plot PIT histogram
hist(pit2, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


# pit histogram for SHASHo mean only
# Get fitted parameters for each observation:
mu_hat_meanonly_SHASHo    <- predict(ni_meanonly_SHASHo, "mu", type = "response")
sigma_hat_meanonly_SHASHo <- predict(ni_meanonly_SHASHo, "sigma", type = "response")
nu_hat_meanonly_SHASHo    <- predict(ni_meanonly_SHASHo, "nu", type = "response")
tau_hat_meanonly_SHASHo   <- predict(ni_meanonly_SHASHo, "tau", type = "response")

# Compute PIT values (using pSEP4 for the SEP4 family, replace with appropriate p* function for your family):
pit3 <- pSHASHo(Nitrate$Nitrate, mu = mu_hat_meanonly_SHASHo, sigma = sigma_hat_meanonly_SHASHo, 
                nu = nu_hat_meanonly_SHASHo, tau = tau_hat_meanonly_SHASHo)

# Plot PIT histogram
hist(pit3, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(pit3)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity


# pit histogram for BCTo
# Get fitted parameters for each observation:
mu_hat_BCTo    <- predict(ni_all_BCTo, "mu", type = "response")
sigma_hat_BCTo <- predict(ni_all_BCTo, "sigma", type = "response")
nu_hat_BCTo    <- predict(ni_all_BCTo, "nu", type = "response")
tau_hat_BCTo   <- predict(ni_all_BCTo, "tau", type = "response")

# Compute PIT values (using pSEP4 for the SEP4 family, replace with appropriate p* function for your family):
pit3 <- pBCTo(Nitrate$Nitrate, mu = mu_hat_BCTo, sigma = sigma_hat_BCTo, nu = nu_hat_BCTo, tau = tau_hat_BCTo)

# Plot PIT histogram
hist(pit3, breaks = 20, main = "PIT Histogram", xlab = "PIT Value", ylab = "Frequency", col = "skyblue")
abline(h = length(pit)/20, col = "red", lwd = 2, lty = 2) # Reference line for uniformity



wp(ni_all_SHASHo, ylim.all = 2.5)
wp(ni_meanonly_SHASHo, ylim.all = 2.5)


resid_wp(ni_all_SHASHo)
resid_wp(ni_meanonly_SHASHo)
resid_wp(ni_nokurt_SHASHo)
resid_wp(ni_noskew_SHASHo)

## STEP 2: ----------
# get the pdf from three random years of the data, as well as the models, and see how well the models
# capture the changes in moments 

# get pdf of three random years --------------
Nitrate_1963 <- Nitrate %>%
  filter(year(Date) == 1963)

Nitrate_1975 <- Nitrate %>%
  filter(year(Date) == 1975)

Nitrate_1994 <- Nitrate %>%
  filter(year(Date) == 1994)


####################################################################
# pdf per year function
pdf_SHASHo_at_year <- function(model, data, year, x = NULL, n = 200,
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
    row_idx <- which(as.Date(data[[date_var]]) == as.Date(date_str))
    if (length(row_idx) == 0) return(NULL)
    newdat <- data[row_idx[1], , drop = FALSE]  # Use the first matching row
    mu    <- predict(model, "mu",    newdata = newdat, type = "response")
    sigma <- predict(model, "sigma", newdata = newdat, type = "response")
    nu    <- predict(model, "nu",    newdata = newdat, type = "response")
    tau   <- predict(model, "tau",   newdata = newdat, type = "response")
    dens  <- dSHASHo(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
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


# mean only
ni_meanonly_1963 <- pdf_SHASHo_at_year(ni_meanonly_SHASHo, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni_meanonly_1975 <- pdf_SHASHo_at_year(ni_meanonly_SHASHo, Nitrate, "1975", time_var = "Date", date_var = "Date")
ni_meanonly_1994 <- pdf_SHASHo_at_year(ni_meanonly_SHASHo, Nitrate, "1994", time_var = "Date", date_var = "Date")
# mean and sigma
ni_mean_and_sigma_1963 <- pdf_SHASHo_at_year(ni_mean_and_sigma_SHASHo, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni_mean_and_sigma_1975 <- pdf_SHASHo_at_year(ni_mean_and_sigma_SHASHo, Nitrate, "1975", time_var = "Date", date_var = "Date")
ni_mean_and_sigma_1994 <- pdf_SHASHo_at_year(ni_mean_and_sigma_SHASHo, Nitrate, "1994", time_var = "Date", date_var = "Date")
# no skew
ni_noskew_1963 <- pdf_SHASHo_at_year(ni_noskew_SHASHo, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni_noskew_1975 <- pdf_SHASHo_at_year(ni_noskew_SHASHo, Nitrate, "1975", time_var = "Date", date_var = "Date")
ni_noskew_1994 <- pdf_SHASHo_at_year(ni_noskew_SHASHo, Nitrate, "1994", time_var = "Date", date_var = "Date")
# no kurt 
ni_nokurt_1963 <- pdf_SHASHo_at_year(ni_nokurt_SHASHo, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni_nokurt_1975 <- pdf_SHASHo_at_year(ni_nokurt_SHASHo, Nitrate, "1975", time_var = "Date", date_var = "Date")
ni_nokurt_1994 <- pdf_SHASHo_at_year(ni_nokurt_SHASHo, Nitrate, "1994", time_var = "Date", date_var = "Date")
# all
ni_all_1963 <- pdf_SHASHo_at_year(ni_all_SHASHo, Nitrate, "1963", time_var = "Date", date_var = "Date")
ni_all_1975 <- pdf_SHASHo_at_year(ni_all_SHASHo, Nitrate, "1975", time_var = "Date", date_var = "Date")
ni_all_1994 <- pdf_SHASHo_at_year(ni_all_SHASHo, Nitrate, "1994", time_var = "Date", date_var = "Date")


par(mfrow = c(1, 3)) 
# plot the models together for year 1963
ni_1963_density <- density(Nitrate_1963$Nitrate)
plot(ni_1963_density, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0,0.18),xlim = c(0,18),
     ylab = "Density", xlab = "Value",
     main = "1963")
lines(ni_meanonly_1963$x, ni_meanonly_1963$avg_density, col = "red", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_1963$x, ni_mean_and_sigma_1963$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_1963$x, ni_noskew_1963$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
lines(ni_nokurt_1963$x, ni_nokurt_1963$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_all_1963$x, ni_all_1963$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")


# plot the models together for year 1975
ni_1975_density <- density(Nitrate_1975$Nitrate)
plot(ni_1975_density, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0,0.08),xlim = c(0,80),
     ylab = "Density", xlab = "Value",
     main = "1975")
lines(ni_meanonly_1975$x, ni_meanonly_1975$avg_density, col = "red", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_1975$x, ni_mean_and_sigma_1975$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_1975$x, ni_noskew_1975$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
lines(ni_nokurt_1975$x, ni_nokurt_1975$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_all_1975$x, ni_all_1975$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")

# plot the models together for year 1994
ni_1994_density <- density(Nitrate_1994$Nitrate)
plot(ni_1994_density, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0,0.04),xlim = c(0,152),
     ylab = "Density", xlab = "Value",
     main = "1994")
lines(ni_meanonly_1994$x, ni_meanonly_1994$avg_density, col = "red", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_1994$x, ni_mean_and_sigma_1994$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_1994$x, ni_noskew_1994$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
lines(ni_nokurt_1994$x, ni_nokurt_1994$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_all_1994$x, ni_all_1994$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
legend("topright", legend = c("Intercept model", "mean(time) only model", "mean(time), sigma(time) & nu(time) model"),
       col = c("black", "steelblue", "goldenrod"), cex = 0.5, lty = c(2,1,1), bty = "n")
par(mfrow = c(1,1)) # reset



############################################ MUTIPLE YEARS
# get pdf multiple years 
# get df of three random years ----------------
Nitrate_62to66 <- Nitrate[Nitrate$year %in% c(1962, 1963, 1964, 1965, 1966), ]
# density plots
ni_62to66_pdf <- density(Nitrate_62to66$Nitrate)
plot(ni_62to66_pdf)

# 1967 till 1970
Nitrate_67to70 <- Nitrate[Nitrate$year %in% c(1967, 1968, 1969, 1970), ]
# density plots
ni_67to70_pdf <- density(Nitrate_67to70$Nitrate)
plot(ni_67to70_pdf)

# 1971 till 1974
Nitrate_71to74 <- Nitrate[Nitrate$year %in% c(1971, 1972, 1973, 1974), ]
# density plots
ni_71to74_pdf <- density(Nitrate_71to74$Nitrate)
plot(ni_71to74_pdf)

# 1975 till 1976
Nitrate_75to78 <- Nitrate[Nitrate$year %in% c(1975, 1976, 1977, 1978), ]
# density plots
ni_75to78_pdf <- density(Nitrate_75to78$Nitrate)
plot(ni_75to78_pdf)

# 1977 till 1982
Nitrate_79to82 <- Nitrate[Nitrate$year %in% c(1979, 1980, 1981, 1982), ]
# density plots
ni_79to82_pdf <- density(Nitrate_79to82$Nitrate)
plot(ni_79to82_pdf)

# 1983 till 1986
Nitrate_83to86 <- Nitrate[Nitrate$year %in% c(1983, 1984, 1985, 1986), ]
# density plots
ni_83to86_pdf <- density(Nitrate_83to86$Nitrate)
plot(ni_83to86_pdf)

# 1987 till 1990
Nitrate_87to90 <- Nitrate[Nitrate$year %in% c(1987, 1988, 1989, 1990), ]
# density plots
ni_87to90_pdf <- density(Nitrate_87to90$Nitrate)
plot(ni_87to90_pdf)

# 1991 till 1994
Nitrate_91to94 <- Nitrate[Nitrate$year %in% c(1991, 1992, 1993, 1994), ]
# density plots
ni_91to94_pdf <- density(Nitrate_91to94$Nitrate)
plot(ni_91to94_pdf)


####################################################################
# pdf multiple years function for for each year as well
pdf_multiple_years <- function(model, data, years, x = NULL, n = 200,
                             time_var = "Date", date_var = "Date") {
  lapply(years, function(y) {
    res <- tryCatch({
      pdf_SHASHo_at_year(model, data, y, x = x, n = n, time_var = time_var, date_var = date_var)
    }, error = function(e) {
      warning(paste("Year", y, "failed:", e$message))
      NULL
    })
    names(res) <- y
    return(res)
  })
}

pdf_multiple_years <- function(model, data, years, x = NULL, n = 200,
                             time_var = "Date", date_var = "Date") {
  results <- lapply(years, function(y) {
    tryCatch({
      pdf_SHASHo_at_year(model, data, y, x = x, n = n, time_var = time_var, date_var = date_var)
    }, error = function(e) {
      warning(paste("Year", y, "failed:", e$message))
      NULL
    })
  })
  names(results) <- as.character(years)
  return(results)
}

######################################################################
# pdf for multiple years combined
pdf_multiple_years_combine <- function(model, data, years, x = NULL, n = 200,
                              time_var = "Date", date_var = "Date") {
  # 1. Identify all unique dates in the desired years
  years <- as.character(years)
  all_dates <- unique(data[[date_var]][format(as.Date(data[[date_var]]), "%Y") %in% years])
  if (length(all_dates) == 0) stop("No dates found for specified years.")
  
  # 2. Create common x grid
  if (is.null(x)) {
    x <- seq(min(model$y, na.rm = TRUE),
             max(model$y, na.rm = TRUE),
             length.out = n)
  }
  
  # 3. For each date, calculate params and PDF, always using full data for prediction
  results <- lapply(all_dates, function(date_str) {
    row_idx <- which(as.Date(data[[date_var]]) == as.Date(date_str))
    if (length(row_idx) == 0) return(NULL)
    newdat <- data[row_idx[1], , drop = FALSE]
    mu    <- predict(model, "mu",    newdata = newdat, type = "response")
    sigma <- predict(model, "sigma", newdata = newdat, type = "response")
    nu    <- predict(model, "nu",    newdata = newdat, type = "response")
    tau   <- predict(model, "tau",   newdata = newdat, type = "response")
    dens  <- dSHASHo(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
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

# get the years we want 
year1 <- c(1962, 1963, 1964, 1965, 1966)
year2 <- c(1967, 1968, 1969, 1970)
year3 <- c(1971, 1972, 1973, 1974)
year4 <- c(1975, 1976, 1977, 1978)
year5 <- c(1979, 1980, 1981, 1982)
year6 <- c(1983, 1984, 1985, 1986)
year7 <- c(1987, 1988, 1989, 1990)
year8 <- c(1991, 1992, 1993, 1994)


# mean only
ni_meanonly_62to66 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_meanonly_67to70 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_meanonly_71to74 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_meanonly_75to78 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_meanonly_79to82 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_meanonly_83to86 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_meanonly_87to90 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_meanonly_91to94 <- pdf_multiple_years_combine(ni_meanonly_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")
# mean and sigma 
ni_mean_and_sigma_62to66 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_67to70 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_71to74 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_75to78 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_79to82 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_83to86 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_87to90 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_91to94 <- pdf_multiple_years_combine(ni_mean_and_sigma_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")
# no skew
ni_noskew_62to66 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_67to70 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_71to74 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_75to78 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_79to82 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_83to86 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_87to90 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_91to94 <- pdf_multiple_years_combine(ni_noskew_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date") # 50 warnings
# no kurt
ni_nokurt_62to66 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_nokurt_67to70 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_nokurt_71to74 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_nokurt_75to78 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_nokurt_79to82 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_nokurt_83to86 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_nokurt_87to90 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_nokurt_91to94 <- pdf_multiple_years_combine(ni_nokurt_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")
# all
ni_all_62to66 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_all_67to70 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_all_71to74 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_all_75to78 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_all_79to82 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_all_83to86 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_all_87to90 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_all_91to94 <- pdf_multiple_years_combine(ni_all_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")



## plot -------------------------------------------------
# plot the models together for year 1963 till 1966
plot(ni_62to66_pdf, col = "black", type = "l", lwd = 2, lty = 2, ylim = c(0, 0.14),
     ylab = "Density", xlab = "Value",
     main = "1962 till 1966")
lines(ni_meanonly_62to66$x, ni_meanonly_62to66$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_62to66$x, ni_mean_and_sigma_62to66$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_62to66$x, ni_noskew_62to66$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_62to66$x, ni_nokurt_62to66$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_62to66$x, ni_all_62to66$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1967 till 1970
plot(ni_67to70_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.11),
     ylab = "Density", xlab = "Value",
     main = "1967 till 1970")
lines(ni_meanonly_67to70$x, ni_meanonly_67to70$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_67to70$x, ni_mean_and_sigma_67to70$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_67to70$x, ni_noskew_67to70$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_67to70$x, ni_nokurt_67to70$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_67to70$x, ni_all_67to70$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1971 till 1974
plot(ni_71to74_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.11),
     ylab = "Density", xlab = "Value",
     main = "1971 till 1974")
lines(ni_meanonly_71to74$x, ni_meanonly_71to74$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_71to74$x, ni_mean_and_sigma_71to74$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_71to74$x, ni_noskew_71to74$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_71to74$x, ni_nokurt_71to74$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_71to74$x, ni_all_71to74$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1975 till 1978
plot(ni_75to78_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1975 till 1978")
lines(ni_meanonly_75to78$x, ni_meanonly_75to78$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_75to78$x, ni_mean_and_sigma_75to78$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_75to78$x, ni_noskew_75to78$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_75to78$x, ni_nokurt_75to78$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_75to78$x, ni_all_75to78$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1979 till 1982
plot(ni_79to82_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1975 till 1978")
lines(ni_meanonly_79to82$x, ni_meanonly_79to82$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_79to82$x, ni_mean_and_sigma_79to82$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_79to82$x, ni_noskew_79to82$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_79to82$x, ni_nokurt_79to82$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_79to82$x, ni_all_79to82$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1983 till 1986
plot(ni_83to86_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1983 till 1986")
lines(ni_meanonly_83to86$x, ni_meanonly_83to86$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_83to86$x, ni_mean_and_sigma_83to86$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_83to86$x, ni_noskew_83to86$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_83to86$x, ni_nokurt_83to86$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_83to86$x, ni_all_83to86$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1987 till 1990
plot(ni_87to90_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1987 till 1990")
lines(ni_meanonly_87to90$x, ni_meanonly_87to90$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_87to90$x, ni_mean_and_sigma_87to90$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_87to90$x, ni_noskew_87to90$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_87to90$x, ni_nokurt_87to90$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_87to90$x, ni_all_87to90$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1991 till 1994
plot(ni_91to94_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1991 till 1994")
lines(ni_meanonly_91to94$x, ni_meanonly_91to94$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_91to94$x, ni_mean_and_sigma_91to94$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_91to94$x, ni_noskew_91to94$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_91to94$x, ni_nokurt_91to94$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_91to94$x, ni_all_91to94$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
lines(ni_all_BCT_91to94$x, ni_all_BCT_91to94$avg_density, col = "pink",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")







#### empirical estimate of how  moments changing through time ####

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

# for each year, get the estimated parametres from the model

allyears <- c(1962:1994)

# mean only
ni_meanonly_allyear_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, allyears, time_var = "Date", date_var = "Date")
ni_meanonly_67to70_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_meanonly_71to74_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_meanonly_75to78_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_meanonly_79to82_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_meanonly_83to86_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_meanonly_87to90_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_meanonly_91to94_param <- pdf_multiple_years(ni_meanonly_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")
# mean and sigma 
ni_mean_and_sigma_62to66_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_67to70_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_71to74_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_75to78_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_79to82_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_83to86_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_87to90_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_91to94_param <- pdf_multiple_years(ni_mean_and_sigma_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")
# no skew
ni_noskew_62to66_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_67to70_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_71to74_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_75to78_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_79to82_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_83to86_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_87to90_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date") # 50 warnings
ni_noskew_91to94_param <- pdf_multiple_years(ni_noskew_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date") # 50 warnings
# no kurt
ni_nokurt_62to66_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_nokurt_67to70_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_nokurt_71to74_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_nokurt_75to78_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_nokurt_79to82_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_nokurt_83to86_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_nokurt_87to90_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_nokurt_91to94_param <- pdf_multiple_years(ni_nokurt_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")
# all
ni_all_62to66_param <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_all_67to70 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_all_71to74 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_all_75to78 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_all_79to82 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_all_83to86 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_all_87to90 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_all_91to94 <- pdf_multiple_years(ni_all_SHASHo, Nitrate, year8, time_var = "Date", date_var = "Date")


# create function to get predicted param from model
best_model_param <- function(model, data){
  
  # get the param_hat 
  data$mu_hat <- predict(model, what = "mu", type = "response")
  data$sigma_hat <- predict(model, what = "sigma", type = "response")
  data$nu_hat <- predict(model, what = "nu", type = "response")
  data$tau_hat <- predict(model, what = "tau", type = "response")
  
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
  return(average_by_year)
  #return(average_by_month)
  
}


# get the predicted params for all models
ni_meanonly_param <- best_model_param(ni_meanonly_SHASHo, Nitrate)
ni_mean_and_sigma_param <- best_model_param(ni_mean_and_sigma_SHASHo, Nitrate)
ni_noskew_param <- best_model_param(ni_noskew_SHASHo, Nitrate)
ni_nokurt_param <- best_model_param(ni_nokurt_SHASHo, Nitrate)
ni_all_param <- best_model_param(ni_all_SHASHo, Nitrate)

ni_all_BCTo_param <- best_model_param(ni_all_BCTo, Nitrate)


# predicted mean
plot(ni_moments_summary$mean ~ ni_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(-4,40),)
lines(ni_meanonly_param$mean_mu_hat ~ ni_meanonly_param$year, col = "darkred", lwd = 2)
lines(ni_mean_and_sigma_param$mean_mu_hat ~ ni_mean_and_sigma_param$year, col = "steelblue", lwd = 2)
lines(ni_noskew_param$mean_mu_hat ~ ni_noskew_param$year, col = "chocolate", lwd = 2)
lines(ni_nokurt_param$mean_mu_hat ~ ni_nokurt_param$year, col = "darkolivegreen", lwd = 2)
lines(ni_all_param$mean_mu_hat ~ ni_all_param$year, col = "goldenrod", lwd = 2)
lines(ni_all_BCTo_param$mean_mu_hat ~ ni_all_BCTo_param$year, col = "pink", lwd = 2)
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("azure4", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# predicted variance
plot(ni_moments_summary$var ~ ni_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(0,1000),)
lines(ni_meanonly_param$mean_sigma_hat ~ ni_meanonly_param$year, col = "darkred", lwd = 2)
lines(ni_mean_and_sigma_param$mean_sigma_hat ~ ni_mean_and_sigma_param$year, col = "steelblue", lwd = 2)
lines(ni_noskew_param$mean_sigma_hat ~ ni_noskew_param$year, col = "chocolate", lwd = 2)
lines(ni_nokurt_param$mean_sigma_hat ~ ni_nokurt_param$year, col = "darkolivegreen", lwd = 2)
lines(ni_all_param$mean_sigma_hat ~ ni_all_param$year, col = "goldenrod", lwd = 2)
lines(ni_all_BCTo_param$mean_sigma_hat ~ ni_all_BCTo_param$year, col = "pink", lwd = 2)
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("azure4", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# predicted skew
plot(ni_moments_summary$skew ~ ni_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(0,5),)
lines(ni_meanonly_param$mean_nu_hat ~ ni_meanonly_param$year, col = "darkred", lwd = 2)
lines(ni_mean_and_sigma_param$mean_nu_hat ~ ni_mean_and_sigma_param$year, col = "steelblue", lwd = 2)
lines(ni_noskew_param$mean_nu_hat ~ ni_noskew_param$year, col = "chocolate", lwd = 2)
lines(ni_nokurt_param$mean_nu_hat ~ ni_nokurt_param$year, col = "darkolivegreen", lwd = 2)
lines(ni_all_param$mean_nu_hat ~ ni_all_param$year, col = "goldenrod", lwd = 2)
lines(ni_all_BCTo_param$mean_nu_hat ~ ni_all_BCTo_param$year, col = "pink", lwd = 2)
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("azure4", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# predicted kurt
plot(ni_moments_summary$kurt ~ ni_moments_summary$year, col = "azure4", type = "l", lwd = 2, lty = 2, ylim = c(-1,20),)
lines(ni_meanonly_param$mean_tau_hat ~ ni_meanonly_param$year, col = "darkred", lwd = 2)
lines(ni_mean_and_sigma_param$mean_tau_hat ~ ni_mean_and_sigma_param$year, col = "steelblue", lwd = 2)
lines(ni_noskew_param$mean_tau_hat ~ ni_noskew_param$year, col = "chocolate", lwd = 2)
lines(ni_nokurt_param$mean_tau_hat ~ ni_nokurt_param$year, col = "darkolivegreen", lwd = 2)
lines(ni_all_param$mean_tau_hat ~ ni_all_param$year, col = "goldenrod", lwd = 2)
lines(ni_all_BCTo_param$mean_tau_hat ~ ni_all_BCTo_param$year, col = "pink", lwd = 2)
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("azure4", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")



# get params for each year
ni_meanonly_62to66[["1965"]]$avg_params
  
  
library(plotly)

plot_ly(
  data = ni_all_param,
  x = ~ year,
  y = ~ mean_sigma_hat,
  z = ~ mean_nu_hat,
  type = 'scatter3d',
  mode = 'markers+lines',
  marker = list(size = 5, color = ~ mean_sigma_hat, colorscale = 'Virdis')
) %>% 
  layout(
    scene = list(
      xaxis = list(title = "Year"),
      yaxis = list(title = "Variance"),
      zaxis = list(title = "Skewness")
    )
  )
  















