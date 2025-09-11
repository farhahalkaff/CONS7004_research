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
  select(Date, Nitrate)

# plot that sucka
ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "black", na.rm = TRUE) +
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()



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


niSSTm_1 <- gamlss(log(Nitrate) ~ Date, family = SST(), data = na.omit(df), 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_1)
with(Nitrate, plot(Nitrate ~ Date))
curve(exp(cbind(1,x)%*%coef(niSSTm_1)), add =T, col = "red", lwd=2)


niSSTm_sig <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, family = SST(), data = na.omit(df), 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_sig)

niSSTm_nu <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, family = SST(), data = na.omit(df), 
                     method = mixed(5, 200),
                     control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))
summary(niSSTm_nu)

niSSTm_tau <- gamlss(log(Nitrate) ~ Date, sigma.fo = ~ Date, nu.fo = ~ Date, tau.fo = ~ Date, family = SST(), 
                     data = na.omit(df), 
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
                 (exp(coefAll(niSSTm_nu)$tau) + 2), "NA", AIC(niSSTm_nu)),
  niSSTm_tau = c(coefAll(niSSTm_tau)$mu[1],coefAll(niSSTm_tau)$mu[2], exp(coefAll(niSSTm_tau)$sigma[1]),
                exp(coefAll(niSSTm_tau)$sigma[2]), exp(coefAll(niSSTm_tau)$nu[1]),exp(coefAll(niSSTm_tau)$nu[2]), 
                (exp(coefAll(niSSTm_tau)$tau[1]) + 2), (exp(coefAll(niSSTm_tau)$tau[2]) + 2), AIC(niSSTm_tau))
)
print(param_summary)


# 
# n <- 6000
# t_index  <- seq_len(n)
# t_scaled <- (t_index - 1)/(n - 1)
# mu0 <- 0
# mu_t <- pmax(0, mu0 + 0.005 * n)
# sigma0 <- 4
# sigma_t <- pmax(0, sigma0 + 1 * n)
# nu0 <- 20
# nu_t <- pmax(0, nu0 + 1 * n)
# tau0 <- 2.1
# 
# y <- rSST(n, mu = mu_t, sigma = sigma_t, nu = nu0, tau = tau0)
# test <- data.frame(
#   t = t_index,
#   t_scaled = n,
#   y = y,
#   sd_true = sigma_t,
#   mu_true = mu_t,
#   nu_true = nu0,
#   tau_true = tau0
# )
#   
# ggplot(test, aes(x = t)) +
#   geom_line(
#     aes(y = y),
#     color = "darkred",
#     linewidth = 0.8
#   ) +
#   labs(
#     x = "Time",
#     y = "Value"
#   ) +
#   theme_minimal(base_size = 14)

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

