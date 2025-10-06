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



Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium, year, month)




am_all_JSU <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, 
                                          nu.fo = ~ year + month, tau.fo = ~ year + month, 
                                          family = JSU(), data = Ammonium,
                                          mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                                          method = mixed(5,200),
                                          control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


am_all_SEP2tr <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, 
                     nu.fo = ~ year + month, tau.fo = ~ year + month, 
                     family = SEP2tr(), data = Ammonium,
                     mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                     method = mixed(5,200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

AIC(am_all_JSU)
AIC(am_all_SEP2)



# bounded family distribution #


am_all_BCTo <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, 
                        nu.fo = ~ year + month, tau.fo = ~ year + month, 
                        family = BCTo(), data = Ammonium,
                        mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium), nu.start = 2, tau.start = 1,
                        method = mixed(5,200),
                        control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

am_all_BCPE <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, 
                      nu.fo = ~ year + month, tau.fo = ~ year + month, 
                      family = BCPE(), data = Ammonium,
                      mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                      method = mixed(5,200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

am_all_BCPEo <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, 
                      nu.fo = ~ year + month, tau.fo = ~ year + month, 
                      family = BCPEo(), data = Ammonium,
                      mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                      method = mixed(5,200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

am_all_GB2 <- gamlss(Ammonium ~ poly(year,2) + month, sigma.fo = ~ year + month, 
                      nu.fo = ~ year + month, tau.fo = ~ year + month, 
                      family = GB2(), data = Ammonium,
                      mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                      method = mixed(5,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))






##### NITRATE ###### BCT
pdf_multiple_years_combine2 <- function(model, data, years, x = NULL, n = 200,
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
    dens  <- dBCTo(x, mu = mu, sigma = sigma, nu = nu, tau = tau)
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


# USE BCTo
ni_meanonly_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1,
                      family = BCTo(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

ni_mean_and_sigma_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ 1, tau.fo = ~ 1,
                      family = BCTo(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

ni_noskew_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ 1, tau.fo = ~ year + month,
                      family = BCTo(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

ni_nokurt_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ 1,
                      family = BCTo(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

ni_all_BCTo <- gamlss(Nitrate ~ year + month, sigma.fo = ~ year + month, nu.fo = ~ year + month, tau.fo = ~ year + month,
                      family = BCTo(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# check AIC
AIC(ni_meanonly_BCTo)
AIC(ni_mean_and_sigma_BCTo)
AIC(ni_noskew_BCTo)
AIC(ni_nokurt_BCTo)
AIC(ni_all_BCTo)



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
ni_meanonly_BCT_62to66 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_67to70 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_71to74 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_75to78 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_79to82 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_83to86 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_87to90 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_meanonly_BCT_91to94 <- pdf_multiple_years_combine2(ni_meanonly_BCTo, Nitrate, year8, time_var = "Date", date_var = "Date")
# mean and sigma 
ni_mean_and_sigma_BCT_62to66 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_67to70 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_71to74 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_75to78 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_79to82 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_83to86 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_87to90 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_mean_and_sigma_BCT_91to94 <- pdf_multiple_years_combine2(ni_mean_and_sigma_BCTo, Nitrate, year8, time_var = "Date", date_var = "Date")
# no skew
ni_noskew_BCT_62to66 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year1, time_var = "Date", date_var = "Date") 
ni_noskew_BCT_67to70 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year2, time_var = "Date", date_var = "Date") 
ni_noskew_BCT_71to74 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year3, time_var = "Date", date_var = "Date") 
ni_noskew_BCT_75to78 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year4, time_var = "Date", date_var = "Date") 
ni_noskew_BCT_79to82 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_noskew_BCT_83to86 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year6, time_var = "Date", date_var = "Date") 
ni_noskew_BCT_87to90 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year7, time_var = "Date", date_var = "Date") 
ni_noskew_BCT_91to94 <- pdf_multiple_years_combine2(ni_noskew_BCTo, Nitrate, year8, time_var = "Date", date_var = "Date") 
# no kurt
ni_nokurt_BCT_62to66 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_67to70 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_71to74 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_75to78 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_79to82 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_83to86 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_87to90 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_nokurt_BCT_91to94 <- pdf_multiple_years_combine2(ni_nokurt_BCTo, Nitrate, year8, time_var = "Date", date_var = "Date")
# all
ni_all_BCT_62to66 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year1, time_var = "Date", date_var = "Date")
ni_all_BCT_67to70 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year2, time_var = "Date", date_var = "Date")
ni_all_BCT_71to74 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year3, time_var = "Date", date_var = "Date")
ni_all_BCT_75to78 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year4, time_var = "Date", date_var = "Date")
ni_all_BCT_79to82 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year5, time_var = "Date", date_var = "Date")
ni_all_BCT_83to86 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year6, time_var = "Date", date_var = "Date")
ni_all_BCT_87to90 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year7, time_var = "Date", date_var = "Date")
ni_all_BCT_91to94 <- pdf_multiple_years_combine2(ni_all_BCTo, Nitrate, year8, time_var = "Date", date_var = "Date")



# plot the models together for year 1963 till 1966
plot(ni_62to66_pdf, col = "black", type = "l", lwd = 2, lty = 2, ylim = c(0, 0.14),
     ylab = "Density", xlab = "Value",
     main = "1962 till 1966")
lines(ni_meanonly_BCT_62to66$x, ni_meanonly_BCT_62to66$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_62to66$x, ni_mean_and_sigma_BCT_62to66$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_62to66$x, ni_noskew_BCT_62to66$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_62to66$x, ni_nokurt_BCT_62to66$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_62to66$x, ni_all_BCT_62to66$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1967 till 1970
plot(ni_67to70_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.11),
     ylab = "Density", xlab = "Value",
     main = "1967 till 1970")
lines(ni_meanonly_BCT_67to70$x, ni_meanonly_BCT_67to70$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_67to70$x, ni_mean_and_sigma_BCT_67to70$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_67to70$x, ni_noskew_BCT_67to70$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_67to70$x, ni_nokurt_BCT_67to70$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_67to70$x, ni_all_BCT_67to70$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1971 till 1974
plot(ni_71to74_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.11),
     ylab = "Density", xlab = "Value",
     main = "1971 till 1974")
lines(ni_meanonly_BCT_71to74$x, ni_meanonly_BCT_71to74$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_71to74$x, ni_mean_and_sigma_BCT_71to74$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_71to74$x, ni_noskew_BCT_71to74$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_71to74$x, ni_nokurt_BCT_71to74$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_71to74$x, ni_all_BCT_71to74$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1975 till 1978
plot(ni_75to78_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1975 till 1978")
lines(ni_meanonly_BCT_75to78$x, ni_meanonly_BCT_75to78$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_75to78$x, ni_mean_and_sigma_BCT_75to78$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_75to78$x, ni_noskew_BCT_75to78$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_75to78$x, ni_nokurt_BCT_75to78$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_75to78$x, ni_all_BCT_75to78$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1979 till 1982
plot(ni_79to82_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1975 till 1978")
lines(ni_meanonly_BCT_79to82$x, ni_meanonly_BCT_79to82$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_79to82$x, ni_mean_and_sigma_BCT_79to82$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_79to82$x, ni_noskew_BCT_79to82$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_79to82$x, ni_nokurt_BCT_79to82$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_79to82$x, ni_all_BCT_79to82$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1983 till 1986
plot(ni_83to86_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1983 till 1986")
lines(ni_meanonly_BCT_83to86$x, ni_meanonly_BCT_83to86$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_83to86$x, ni_mean_and_sigma_BCT_83to86$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_83to86$x, ni_noskew_BCT_83to86$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_83to86$x, ni_nokurt_BCT_83to86$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_83to86$x, ni_all_BCT_83to86$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1987 till 1990
plot(ni_87to90_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1987 till 1990")
lines(ni_meanonly_BCT_87to90$x, ni_meanonly_BCT_87to90$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_87to90$x, ni_mean_and_sigma_BCT_87to90$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_87to90$x, ni_noskew_BCT_87to90$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_87to90$x, ni_nokurt_BCT_87to90$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_87to90$x, ni_all_BCT_87to90$avg_density, col = "goldenrod",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")

# plot the models together for year 1991 till 1994
plot(ni_91to94_pdf, col = "black", type = "l", lwd = 2, lty = 2,ylim = c(0, 0.08),
     ylab = "Density", xlab = "Value",
     main = "1991 till 1994")
lines(ni_meanonly_BCT_91to94$x, ni_meanonly_BCT_91to94$avg_density, col = "darkred", lwd = 2, lty = 1) 
lines(ni_mean_and_sigma_BCT_91to94$x, ni_mean_and_sigma_BCT_91to94$avg_density, col = "steelblue", lwd = 2, lty = 1) # mean and sigma
lines(ni_noskew_BCT_91to94$x, ni_noskew_BCT_91to94$avg_density, col = "chocolate",lwd = 2, lty = 1) 
lines(ni_nokurt_BCT_91to94$x, ni_nokurt_BCT_91to94$avg_density, col = "darkolivegreen",lwd = 2, lty = 1) 
lines(ni_all_BCT_91to94$x, ni_all_BCT_91to94$avg_density, col = "pink",lwd = 2, lty = 1) 
legend("topright", legend = c("empirical", "mean(time) only model", "mean(time), sigma(time)", 
                              "constant skew", "constant kurtosis", "all(time)"),
       col = c("black", "darkred", "steelblue",
               "chocolate", "darkolivegreen", "goldenrod"), cex = 0.7, lty = c(2,1,1,1,1,1), bty = "n")





# compare with other family distirbution 
# two-parameter
# NO
ni_all_NOtr <- gamlss(Nitrate ~ year + month,
                      family = NOtr(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# LO
ni_all_LOtr <- gamlss(Nitrate ~ year + month,
                      family = LOtr(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# RG
ni_all_RGtr <- gamlss(Nitrate ~ year + month,
                      family = RGtr(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# exGAUS
ni_all_exGAUStr <- gamlss(Nitrate ~ year + month, # doesnt work 
                      family = exGAUStr(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE
ni_all_PEtr <- gamlss(Nitrate ~ year + month, 
                          family = PEtr(), data = Nitrate,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
ni_all_SN1tr <- gamlss(Nitrate ~ year + month, 
                       family = SN1tr(), data = Nitrate,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN2
ni_all_SN2tr <- gamlss(Nitrate ~ year + month, 
                      family = SN2tr(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
ni_all_TFtr <- gamlss(Nitrate ~ year + month, 
                       family = TFtr(), data = Nitrate,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF2
ni_all_TF2tr <- gamlss(Nitrate ~ year + month, 
                       family = TF2tr(), data = Nitrate,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
ni_all_GTtr <- gamlss(Nitrate ~ year + month,
                       family = GTtr(), data = Nitrate,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSU
ni_all_JSUtr <- gamlss(Nitrate ~ year + month,
                      family = JSUtr(), data = Nitrate,
                      method = mixed(10,200),
                      control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
ni_all_JSUotr <- gamlss(Nitrate ~ year + month,
                       family = JSUotr(), data = Nitrate,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# NET
ni_all_JSUotr <- gamlss(Nitrate ~ year + month,
                        family = JSUotr(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))








