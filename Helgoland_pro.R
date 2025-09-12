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


# plot that sucka
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "black", na.rm = TRUE) +
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()

# plot a specific year 
ggplot(dplyr::filter(Nitrate, year(Date) == 1994),
       aes(Date, Nitrate)) +
  geom_line(color = "black")






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


### MODELS ###

niSSTm <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niTF2m <- gamlss(Nitrate ~ Date, family = TF2(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niSN1m <- gamlss(Nitrate ~ Date, family = SN1(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

niNOm <- gamlss(Nitrate ~ Date, family = NO(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


summary(niSSTm)
summary(niTF2m)
AIC(niTF2m)
AIC(niTF2mtest)


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

omit <- na.omit(df)

##########=====================================================
# niSSTm_1 <- gamlss(Nitrate ~ Date, family = SST(), data = na.omit(df), 
#                    method = mixed(5, 200),
#                    control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_1)
# with(Nitrate, plot(Nitrate ~ Date))
# curve(cbind(1,x)%*%coef(niSSTm_1), add =T, col = "red", lwd=2)
# 
# 
# niSSTm_sig <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, family = SST(), data = na.omit(df), 
#                      method = mixed(5, 200),
#                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_sig)
# 
# niSSTm_nu <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = na.omit(df), 
#                     method = mixed(5, 200),
#                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_nu)
# 
# niSSTm_tau <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), 
#                      data = na.omit(df), 
#                      method = mixed(5, 200),
#                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_tau)
##########======================================================


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


#############################################
############################################# LOG

niSSTm_1b <- gamlss(log(Nitrate) ~ Date, family = SST(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_1b)
with(Nitrate, plot(Nitrate ~ Date))
curve(exp(cbind(1,x)%*%coef(niSSTm_1b)), add =T, col = "red", lwd=2)


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


##########=========================================================
# niSSTm_1b <- gamlss(log(Nitrate) ~ Date, family = SST(), data = na.omit(df), 
#                     method = mixed(5, 200),
#                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_1b)
# with(Nitrate, plot(Nitrate ~ Date))
# curve(exp(cbind(1,x)%*%coef(niSSTm_1b)), add =T, col = "red", lwd=2)
# 
# 
# niSSTm_sigb <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, family = SST(), data = na.omit(df), 
#                       method = mixed(5, 200),
#                       control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_sigb)
# 
# niSSTm_nub <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = na.omit(df), 
#                      method = mixed(5, 200),
#                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_nub)
# 
# niSSTm_taub <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), 
#                       data = na.omit(df), 
#                       method = mixed(5, 200),
#                       control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
# summary(niSSTm_taub)
############==========================================================

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



#############################################
############################################# BOX-COX

niSSTm_1c <- gamlss(Nitrate ~ Date, family = BCT(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_1)
with(Nitrate, plot(Nitrate ~ Date))
curve(cbind(1,x)%*%coef(niSSTm_1), add =T, col = "red", lwd=2)


niSSTm_sigc <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, family = BCT(), data = Nitrate, 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_sig)

niSSTm_nuc <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = BCT(), data = Nitrate, 
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
               exp(coefAll(niSSTm_1c)$nu),"NA", (exp(coefAll(niSSTm_1c)$tau) + 2), "NA", AIC(niSSTm_1c)),
  niSSTm_sigc = c(coefAll(niSSTm_sigc)$mu[1],coefAll(niSSTm_sigc)$mu[2] ,exp(coefAll(niSSTm_sigc)$sigma[1]),
                 exp(coefAll(niSSTm_sigc)$sigma[2]), exp(coefAll(niSSTm_sigc)$nu),"NA", 
                 (exp(coefAll(niSSTm_sigc)$tau) + 2), "NA", AIC(niSSTm_sigc)),
  niSSTm_nuc = c(coefAll(niSSTm_nuc)$mu[1],coefAll(niSSTm_nuc)$mu[2], exp(coefAll(niSSTm_nuc)$sigma[1]),
                exp(coefAll(niSSTm_nuc)$sigma[2]), exp(coefAll(niSSTm_nuc)$nu[1]),exp(coefAll(niSSTm_nuc)$nu[2]), 
                (exp(coefAll(niSSTm_nuc)$tau) + 2), "NA", AIC(niSSTm_nuc)),
  niSSTm_tauc = c(coefAll(niSSTm_tauc)$mu[1],coefAll(niSSTm_tauc)$mu[2], exp(coefAll(niSSTm_tauc)$sigma[1]),
                 exp(coefAll(niSSTm_tauc)$sigma[2]), exp(coefAll(niSSTm_tauc)$nu[1]),exp(coefAll(niSSTm_tauc)$nu[2]), 
                 (exp(coefAll(niSSTm_tauc)$tau[1]) + 2), (exp(coefAll(niSSTm_tauc)$tau[2]) + 2), AIC(niSSTm_tauc))
)
print(param_summary)




#############################################
############################################# Skew exponential 

niSEPm_1 <- gamlss(Nitrate ~ Date, family = SEP3(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSEPm_1)


niSEPm_sig <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, family = SEP3(), data = Nitrate, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))


niSEPm_nu <- gamlss(Nitrate ~ Date, sigma.fo = ~ Date, family = SEP3(), data = Nitrate, 
                      method = mixed(5, 200),
                      control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

# summary 
param_summary <- data.frame(
  param = c("mu","mu(time)" , "SD","SD(time)","nu", "nu(time)", "tau", "tau(time)", "AIC"), 
  niSEPm_1 = c(coefAll(niSEPm_1)$mu[1],coefAll(niSEPm_1)$mu[2] ,exp(coefAll(niSEPm_1)$sigma),"NA", 
                exp(coefAll(niSEPm_1)$nu),"NA", (exp(coefAll(niSEPm_1)$tau) + 2), "NA", AIC(niSEPm_1)),
  niSEPm_sig = c(coefAll(niSEPm_sig)$mu[1],coefAll(niSEPm_sig)$mu[2] ,exp(coefAll(niSEPm_sig)$sigma[1]),
                  exp(coefAll(niSEPm_sig)$sigma[2]), exp(coefAll(niSEPm_sig)$nu),"NA", 
                  (exp(coefAll(niSEPm_sig)$tau) + 2), "NA", AIC(niSEPm_sig)),
  niSEPm_nu = c(coefAll(niSEPm_nu)$mu[1],coefAll(niSEPm_nu)$mu[2], exp(coefAll(niSEPm_nu)$sigma[1]),
                 exp(coefAll(niSEPm_nu)$sigma[2]), exp(coefAll(niSEPm_nu)$nu[1]),exp(coefAll(niSEPm_nu)$nu[2]), 
                 (exp(coefAll(niSEPm_nu)$tau) + 2), "NA", AIC(niSEPm_nu))
)
print(param_summary)


#############################################
############################################# SST left-truncated

library(gamlss.tr)
gen.trun(0,"SST",type="left")

niSSTm_1 <- gamlss(Nitrate ~ Date, family = SSTtr(), data = Nitrate, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 900, c.crit = 0.01, trace = FALSE))








######################################################
### Phosphate and fit models of different families ###
######################################################

# create own df to remove NAs
bluh <- df %>% filter(!is.na(Phosphate)) %>% 
  select(Date, Phosphate)

# plot that sucka
ggplot(bluh, aes(x = Date, y = Phosphate)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()


### CHECK DISTRIBUTION ###

# histogram 
hist(bluh$Phosphate)

# density plot 
ggplot(bluh, aes(y = Phosphate)) +
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

phSSTm <- gamlss(Phosphate ~ poly(Date,2), family = SST(), data = na.omit(df), 
                 sigma.start = max(sd(df$Phosphate, na.rm = TRUE), 1e-6),
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


####################################################
### Nitrite and fit models of different families ###
####################################################

# Create own df to remove NAs
Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite)

# plot that sucka
ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "black") + 
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()


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


################################################
### DIN and fit models of different families ###
################################################

# create it's own df to remove NAs
DIN <- df %>% filter(!is.na(DIN)) %>% 
  select(Date, DIN)

# plot that sucka
ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "DIN (µmol/l)") +
  theme_minimal()


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


DINSSTm_nu <- gamlss(log(DIN) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = DIN, 
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


#####################################################
### Silicate and fit models of different families ###
#####################################################

# create own df to remove NAs 
Silicate <- df %>% filter(!is.na(Silicate)) %>% 
  select(Date, Silicate)

# plot that sucka 
ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()


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


#####################################################
### Ammonium and fit models of different families ###
#####################################################

# create own df to remove NAs 
Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium)

# plot that sucka
ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()


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

amSSTm <- gamlss(Ammonium ~ Date, family = SST(), data = Ammonium, 
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

