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

Nitrate <- df %>% filter(!is.na(Nitrate)) %>% 
  select(Date, Nitrate)

ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
  geom_line(color = "black", na.rm = TRUE) +
  labs(x = "Time", y = "Nitrate (µmol/l)") +
  theme_minimal()

niSSTm <- gamlss(Nitrate ~ Date, family = SST(), data = na.omit(df), 
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

summary(m1)
summary(m2)

AIC(niSSTm)
AIC(niTF2m)
AIC(niSN1m)
AIC(niNOm)


######################################################
### Phosphate and fit models of different families ###
######################################################

bluh <- df %>% filter(!is.na(Phosphate)) %>% 
  select(Date, Phosphate)


ggplot(bluh, aes(x = Date, y = Phosphate)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Phosphate (µmol/l)") +
  theme_minimal()

phSSTm <- gamlss(Phosphate ~ poly(Date,2), family = SST(), data = bluh, 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phTF2m <- gamlss(Phosphate ~ poly(Date,2), family = TF2(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phSN1m <- gamlss(Phosphate ~ poly(Date,2), family = SN1(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))

phNOm <- gamlss(Phosphate ~ poly(Date,1), family = NO(), data = na.omit(df), 
             method = mixed(5, 200),
             control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE))



AIC(phTF2m)
AIC(phSN1m)
AIC(phNOm)


####################################################
### Nitrite and fit models of different families ###
####################################################

Nitrite <- df %>% filter(!is.na(Nitrite)) %>% 
  select(Date, Nitrite)


ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Nitrite (µmol/l)") +
  theme_minimal()


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

DIN <- df %>% filter(!is.na(DIN)) %>% 
  select(Date, DIN)


ggplot(DIN, aes(x = Date, y = DIN)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "DIN (µmol/l)") +
  theme_minimal()


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

Silicate <- df %>% filter(!is.na(Silicate)) %>% 
  select(Date, Silicate)


ggplot(Silicate, aes(x = Date, y = Silicate)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Silicate (µmol/l)") +
  theme_minimal()


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


AIC(siSSTm)
AIC(siTF2m)
AIC(siSN1m)
AIC(siNOm)



#####################################################
### Ammonium and fit models of different families ###
#####################################################

Ammonium <- df %>% filter(!is.na(Ammonium)) %>% 
  select(Date, Ammonium)


ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
  geom_line(color = "black") +
  labs(x = "Time", y = "Ammonium (µmol/l)") +
  theme_minimal()


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

