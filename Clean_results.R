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



########################################
##### HISTOGRAM WITH DENSITY PLOT ######
########################################

p1 <- ggplot(Nitrate, aes(x = Nitrate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 10) + 
  xlab("Nitrate (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

p2 <- ggplot(Nitrite, aes(x = Nitrite)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  xlab("Nitrite (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,6)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

p3 <- ggplot(Ammonium, aes(x = Ammonium)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("Ammonium (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,35)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

p4 <- ggplot(DIN, aes(x = DIN)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 10) + 
  xlab("DIN (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

p5 <- ggplot(Silicate, aes(x = Silicate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 5) + 
  xlab("Silicate (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

p6 <- ggplot(Phosphate, aes(x = Phosphate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.2) + 
  xlab("Phosphate (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,2.5)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

# stack em together
 p1 + p2 + p3 + p4 + p5 + p6

 
 
 ###########################################
 ###### RAW PLOT WITH INTERCEPT MODEL ######
 ###########################################
 
 # Nitrate
 pp1 <- ggplot(Nitrate, aes(x = Date, y = Nitrate)) +
   geom_line(color = "azure4", na.rm = TRUE) +
   geom_hline(yintercept = mean(Nitrate$Nitrate), linetype = "dashed", color = "black") + # intercept mean
   labs(x = "Time", y = "Nitrate (µmol/l)") +
   theme_bw()

 # Nitrite
 pp2 <- ggplot(Nitrite, aes(x = Date, y = Nitrite)) +
   geom_line(color = "azure4", na.rm = TRUE) +
   geom_hline(yintercept = mean(Nitrite$Nitrite), linetype = "dashed", color = "black") + # intercept mean
   labs(x = "Time", y = "Nitrite (µmol/l)") +
   theme_bw()
 
 #Ammonium 
 pp3 <- ggplot(Ammonium, aes(x = Date, y = Ammonium)) +
   geom_line(color = "azure4", na.rm = TRUE) +
   geom_hline(yintercept = mean(Ammonium$Ammonium), linetype = "dashed", color = "black") + # intercept mean
   labs(x = "Time", y = "Ammonium (µmol/l)") +
   theme_bw()
 
 # DIN
 pp4 <- ggplot(DIN, aes(x = Date, y = DIN)) +
   geom_line(color = "azure4", na.rm = TRUE) +
   geom_hline(yintercept = mean(DIN$DIN), linetype = "dashed", color = "black") + # intercept mean
   labs(x = "Time", y = "DIN (µmol/l)") +
   theme_bw()
 
 # Silicate
 pp5 <- ggplot(Silicate, aes(x = Date, y = Silicate)) +
   geom_line(color = "azure4", na.rm = TRUE) +
   geom_hline(yintercept = mean(Silicate$Silicate), linetype = "dashed", color = "black") + # intercept mean
   labs(x = "Time", y = "Silicate (µmol/l)") +
   theme_bw()
 
 # Phosphate
 pp6 <- ggplot(Phosphate, aes(x = Date, y = Phosphate)) +
   geom_line(color = "azure4", na.rm = TRUE) +
   geom_hline(yintercept = mean(Phosphate$Phosphate), linetype = "dashed", color = "black") + # intercept mean
   labs(x = "Time", y = "Phosphate (µmol/l)") +
   theme_bw()

# stack em 
layout_matrix <- matrix(c(1, 2,
                          3, 4,
                          4, 6), ncol = 2, byrow = TRUE)
 
rawplot <- plot_grid(pp1, pp2, pp3, pp4, pp5, pp6, 
          ncol = 2, nrow = 3)
ggsave("rawplots.pdf", plot = rawplot, width = 8, height = 6, units = "in", dpi = 300)



##################################################
#### intercept models with different families ####
##################################################


#### NITRATE #### ========================================================================================

# NO
ni_static_NO <- gamlss(Nitrate ~ 1, family = NO(), data = Nitrate) # two-parameter
# NO2
ni_static_NO2 <- gamlss(Nitrate ~ 1, family = NO2(), data = Nitrate)
# GU
ni_static_GU <- gamlss(Nitrate ~ 1, family = GU(), data = Nitrate)
# LO
ni_static_LO <- gamlss(Nitrate ~ 1, family = LO(), data = Nitrate)
# RG
ni_static_RG <- gamlss(Nitrate ~ 1, family = RG(), data = Nitrate)
# exGAUS
ni_static_exGAUS <- gamlss(Nitrate ~ 1, family = exGAUS(), data = Nitrate) # three-parameter
# NOF
ni_static_NOF <- gamlss(Nitrate ~ 1, family = NOF(), data = Nitrate)
# PE
ni_static_PE <- gamlss(Nitrate ~ 1, family = PE(), data = Nitrate)
# PE2
ni_static_PE2 <- gamlss(Nitrate ~ 1, family = PE2(), data = Nitrate, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
ni_static_SN1 <- gamlss(Nitrate ~ 1, family = SN1(), data = Nitrate)
# TF
ni_static_TF <- gamlss(Nitrate ~ 1, family = TF(), data = Nitrate)
# TF2
ni_static_TF2 <- gamlss(Nitrate ~ 1, family = TF2(), data = Nitrate)
# GT
ni_static_GT <- gamlss(Nitrate ~ 1, family = GT(), data = Nitrate, 
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
ni_static_JSU <- gamlss(Nitrate ~ 1, family = JSU(), data = Nitrate, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
ni_static_JSUo <- gamlss(Nitrate ~ 1, family = JSUo(), data = Nitrate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# NET
ni_static_NET <- gamlss(Nitrate ~ 1, family = NET(), data = Nitrate)
# SHASH
ni_static_SHASH <- gamlss(Nitrate ~ 1, family = SHASH(), data = Nitrate,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
ni_static_SHASHo <- gamlss(Nitrate ~ 1, family = SHASHo(), data = Nitrate,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
ni_static_SHASHo2 <- gamlss(Nitrate ~ 1, family = SHASHo2(), data = Nitrate,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
ni_static_SEP2 <- gamlss(Nitrate ~ 1, family = SEP2(), data = Nitrate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
ni_static_SEP3 <- gamlss(Nitrate ~ 1, family = SEP3(), data = Nitrate,
                         mu.start = mean(Nitrate$Nitrate), sigma.start = sd(Nitrate$Nitrate),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue #
# SEP4
ni_static_SEP4 <- gamlss(Nitrate ~ 1, family = SEP4(), data = Nitrate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue
# SST
ni_static_SST <- gamlss(Nitrate ~ 1, family = SST(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST1
ni_static_ST1 <- gamlss(Nitrate ~ 1, family = ST1(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
ni_static_ST2 <- gamlss(Nitrate ~ 1, family = ST2(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
ni_static_ST3 <- gamlss(Nitrate ~ 1, family = ST3(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
ni_static_ST4 <- gamlss(Nitrate ~ 1, family = ST4(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
ni_static_ST5 <- gamlss(Nitrate ~ 1, family = ST5(), data = Nitrate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
ni_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2", "SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(ni_static_NO, k = 2),GAIC(ni_static_NO2, k = 2), GAIC(ni_static_GU, k = 2), GAIC(ni_static_LO, k = 2),
          GAIC(ni_static_RG, k = 2),GAIC(ni_static_exGAUS, k = 2), GAIC(ni_static_NOF, k = 2), GAIC(ni_static_PE, k = 2),
          GAIC(ni_static_PE2, k = 2),GAIC(ni_static_SN1, k = 2), GAIC(ni_static_TF, k = 2), GAIC(ni_static_TF2, k = 2),
          GAIC(ni_static_GT, k = 2),GAIC(ni_static_JSU, k = 2), GAIC(ni_static_JSUo, k = 2), GAIC(ni_static_NET, k = 2),
          GAIC(ni_static_SHASH, k = 2),GAIC(ni_static_SHASHo, k = 2), GAIC(ni_static_SHASHo2, k = 2), GAIC(ni_static_SEP2, k = 2),
          GAIC(ni_static_SST, k = 2), GAIC(ni_static_ST1, k = 2),
          GAIC(ni_static_ST2, k = 2),GAIC(ni_static_ST3, k = 2), GAIC(ni_static_ST4, k = 2), GAIC(ni_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4")
  )


# Plot

custom_colors <- c("2" = "gray88", "3" = "grey48", "4" = "black")

pAIC_1 <- ggplot(ni_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors, name = "param") +
  labs(x = "AIC", y = "Intercept models (Nitrate)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



#### NITRITE #### ========================================================================================

# NO
nii_static_NO <- gamlss(Nitrite ~ 1, family = NO(), data = Nitrite) # two-parameter
# NO2
nii_static_NO2 <- gamlss(Nitrite ~ 1, family = NO2(), data = Nitrite)
# GU
nii_static_GU <- gamlss(Nitrite ~ 1, family = GU(), data = Nitrite)
# LO
nii_static_LO <- gamlss(Nitrite ~ 1, family = LO(), data = Nitrite)
# RG
nii_static_RG <- gamlss(Nitrite ~ 1, family = RG(), data = Nitrite)
# exGAUS
nii_static_exGAUS <- gamlss(Nitrite ~ 1, family = exGAUS(), data = Nitrite) # three-parameter
# NOF
nii_static_NOF <- gamlss(Nitrite ~ 1, family = NOF(), data = Nitrite)
# PE
nii_static_PE <- gamlss(Nitrite ~ 1, family = PE(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
nii_static_PE2 <- gamlss(Nitrite ~ 1, family = PE2(), data = Nitrite, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
nii_static_SN1 <- gamlss(Nitrite ~ 1, family = SN1(), data = Nitrite)
# TF
nii_static_TF <- gamlss(Nitrite ~ 1, family = TF(), data = Nitrite)
# TF2
nii_static_TF2 <- gamlss(Nitrite ~ 1, family = TF2(), data = Nitrite,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
nii_static_GT <- gamlss(Nitrite ~ 1, family = GT(), data = Nitrite, 
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
nii_static_JSU <- gamlss(Nitrite ~ 1, family = JSU(), data = Nitrite, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
nii_static_JSUo <- gamlss(Nitrite ~ 1, family = JSUo(), data = Nitrite,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# NET
nii_static_NET <- gamlss(Nitrite ~ 1, family = NET(), data = Nitrite)
# SHASH
nii_static_SHASH <- gamlss(Nitrite ~ 1, family = SHASH(), data = Nitrite,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
nii_static_SHASHo <- gamlss(Nitrite ~ 1, family = SHASHo(), data = Nitrite,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
nii_static_SHASHo2 <- gamlss(Nitrite ~ 1, family = SHASHo2(), data = Nitrite,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
nii_static_SEP2 <- gamlss(Nitrite ~ 1, family = SEP2(), data = Nitrite,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
nii_static_SEP3 <- gamlss(Nitrite ~ 1, family = SEP3(), data = Nitrite,
                         mu.start = mean(Nitrite$Nitrite), sigma.start = sd(Nitrite$Nitrite),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue #
# SEP4
nii_static_SEP4 <- gamlss(Nitrite ~ 1, family = SEP4(), data = Nitrite,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue
# SST
nii_static_SST <- gamlss(Nitrite ~ 1, family = SST(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST1
nii_static_ST1 <- gamlss(Nitrite ~ 1, family = ST1(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
nii_static_ST2 <- gamlss(Nitrite ~ 1, family = ST2(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
nii_static_ST3 <- gamlss(Nitrite ~ 1, family = ST3(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
nii_static_ST4 <- gamlss(Nitrite ~ 1, family = ST4(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
nii_static_ST5 <- gamlss(Nitrite ~ 1, family = ST5(), data = Nitrite,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
nii_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "SEP4" ,"SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(nii_static_NO, k = 2),GAIC(nii_static_NO2, k = 2), GAIC(nii_static_GU, k = 2), GAIC(nii_static_LO, k = 2), GAIC(nii_static_RG, k = 2),
          GAIC(nii_static_exGAUS, k = 2), GAIC(nii_static_NOF, k = 2), GAIC(nii_static_PE, k = 2), GAIC(nii_static_PE2, k = 2),GAIC(nii_static_SN1, k = 2), GAIC(nii_static_TF, k = 2), GAIC(nii_static_TF2, k = 2),
          GAIC(nii_static_GT, k = 2),GAIC(nii_static_JSU, k = 2), GAIC(nii_static_JSUo, k = 2), GAIC(nii_static_NET, k = 2),
          GAIC(nii_static_SHASH, k = 2),GAIC(nii_static_SHASHo, k = 2), GAIC(nii_static_SHASHo2, k = 2), GAIC(nii_static_SEP2, k = 2), GAIC(nii_static_SEP3, k = 2),
          GAIC(nii_static_SEP4, k = 2), GAIC(nii_static_SST, k = 2), GAIC(nii_static_ST1, k = 2),
          GAIC(nii_static_ST2, k = 2),GAIC(nii_static_ST3, k = 2), GAIC(nii_static_ST4, k = 2), GAIC(nii_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4")
)



# Plot
pAIC_2 <- ggplot(nii_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors, name = "param") +
  labs(x = "AIC", y = "Intercept models (Nitrite)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


#### AMMONIUM #### ========================================================================================

# NO
am_static_NO <- gamlss(Ammonium ~ 1, family = NO(), data = Ammonium) # two-parameter
# NO2
am_static_NO2 <- gamlss(Ammonium ~ 1, family = NO2(), data = Ammonium)
# GU
am_static_GU <- gamlss(Ammonium ~ 1, family = GU(), data = Ammonium)
# LO
am_static_LO <- gamlss(Ammonium ~ 1, family = LO(), data = Ammonium)
# RG
am_static_RG <- gamlss(Ammonium ~ 1, family = RG(), data = Ammonium)
# exGAUS
am_static_exGAUS <- gamlss(Ammonium ~ 1, family = exGAUS(), data = Ammonium,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
am_static_NOF <- gamlss(Ammonium ~ 1, family = NOF(), data = Ammonium)
# PE
am_static_PE <- gamlss(Ammonium ~ 1, family = PE(), data = Ammonium,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
am_static_PE2 <- gamlss(Ammonium ~ 1, family = PE2(), data = Ammonium, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
am_static_SN1 <- gamlss(Ammonium ~ 1, family = SN1(), data = Ammonium)
# TF
am_static_TF <- gamlss(Ammonium ~ 1, family = TF(), data = Ammonium)
# TF2
am_static_TF2 <- gamlss(Ammonium ~ 1, family = TF2(), data = Ammonium,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
am_static_GT <- gamlss(Ammonium ~ 1, family = GT(), data = Ammonium, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
am_static_JSU <- gamlss(Ammonium ~ 1, family = JSU(), data = Ammonium, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
am_static_JSUo <- gamlss(Ammonium ~ 1, family = JSUo(), data = Ammonium,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# NET
am_static_NET <- gamlss(Ammonium ~ 1, family = NET(), data = Ammonium)
# SHASH
am_static_SHASH <- gamlss(Ammonium ~ 1, family = SHASH(), data = Ammonium,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
am_static_SHASHo <- gamlss(Ammonium ~ 1, family = SHASHo(), data = Ammonium,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
am_static_SHASHo2 <- gamlss(Ammonium ~ 1, family = SHASHo2(), data = Ammonium,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
am_static_SEP2 <- gamlss(Ammonium ~ 1, family = SEP2(), data = Ammonium,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
am_static_SEP3 <- gamlss(Ammonium ~ 1, family = SEP3(), data = Ammonium,
                          mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
am_static_SEP4 <- gamlss(Ammonium ~ 1, family = SEP4(), data = Ammonium,
                         mu.start = mean(Ammonium$Ammonium), sigma.start = sd(Ammonium$Ammonium),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue
# ST1
am_static_ST1 <- gamlss(Ammonium ~ 1, family = ST1(), data = Ammonium,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
am_static_ST2 <- gamlss(Ammonium ~ 1, family = ST2(), data = Ammonium,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
am_static_ST3 <- gamlss(Ammonium ~ 1, family = ST3(), data = Ammonium,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
am_static_ST4 <- gamlss(Ammonium ~ 1, family = ST4(), data = Ammonium,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
am_static_ST5 <- gamlss(Ammonium ~ 1, family = ST5(), data = Ammonium,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
am_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(am_static_NO, k = 2),GAIC(am_static_NO2, k = 2), GAIC(am_static_GU, k = 2), GAIC(am_static_LO, k = 2), 
          GAIC(am_static_RG, k = 2), GAIC(am_static_exGAUS, k = 2), GAIC(am_static_NOF, k = 2), GAIC(am_static_PE, k = 2), 
          GAIC(am_static_PE2, k = 2),GAIC(am_static_SN1, k = 2), GAIC(am_static_TF, k = 2), GAIC(am_static_TF2, k = 2),
          GAIC(am_static_GT, k = 2),GAIC(am_static_JSU, k = 2), GAIC(am_static_JSUo, k = 2), GAIC(am_static_NET, k = 2),
          GAIC(am_static_SHASH, k = 2),GAIC(am_static_SHASHo, k = 2), GAIC(am_static_SHASHo2, k = 2), GAIC(am_static_SEP2, k = 2), 
          GAIC(am_static_SEP3, k = 2), GAIC(am_static_ST1, k = 2),
          GAIC(am_static_ST2, k = 2),GAIC(am_static_ST3, k = 2), GAIC(am_static_ST4, k = 2), GAIC(am_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4")
)


# Plot
pAIC_3 <- ggplot(am_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors, name = "param") +
  labs(x = "AIC", y = "Intercept models (Ammonium)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



#### DIN #### ========================================================================================

# NO
DIN_static_NO <- gamlss(DIN ~ 1, family = NO(), data = DIN) # two-parameter
# NO2
DIN_static_NO2 <- gamlss(DIN ~ 1, family = NO2(), data = DIN)
# GU
DIN_static_GU <- gamlss(DIN ~ 1, family = GU(), data = DIN)
# LO
DIN_static_LO <- gamlss(DIN ~ 1, family = LO(), data = DIN)
# RG
DIN_static_RG <- gamlss(DIN ~ 1, family = RG(), data = DIN)
# exGAUS
DIN_static_exGAUS <- gamlss(DIN ~ 1, family = exGAUS(), data = DIN,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
DIN_static_NOF <- gamlss(DIN ~ 1, family = NOF(), data = DIN)
# PE
DIN_static_PE <- gamlss(DIN ~ 1, family = PE(), data = DIN,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
DIN_static_PE2 <- gamlss(DIN ~ 1, family = PE2(), data = DIN, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
DIN_static_SN1 <- gamlss(DIN ~ 1, family = SN1(), data = DIN)
# TF
DIN_static_TF <- gamlss(DIN ~ 1, family = TF(), data = DIN)
# TF2
DIN_static_TF2 <- gamlss(DIN ~ 1, family = TF2(), data = DIN,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
DIN_static_GT <- gamlss(DIN ~ 1, family = GT(), data = DIN, 
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
DIN_static_JSU <- gamlss(DIN ~ 1, family = JSU(), data = DIN, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
DIN_static_JSUo <- gamlss(DIN ~ 1, family = JSUo(), data = DIN,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# NET
DIN_static_NET <- gamlss(DIN ~ 1, family = NET(), data = DIN)
# SHASH
DIN_static_SHASH <- gamlss(DIN ~ 1, family = SHASH(), data = DIN,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
DIN_static_SHASHo <- gamlss(DIN ~ 1, family = SHASHo(), data = DIN,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
DIN_static_SHASHo2 <- gamlss(DIN ~ 1, family = SHASHo2(), data = DIN,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
DIN_static_SEP2 <- gamlss(DIN ~ 1, family = SEP2(), data = DIN,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
DIN_static_SEP3 <- gamlss(DIN ~ 1, family = SEP3(), data = DIN,
                         mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
DIN_static_SEP4 <- gamlss(DIN ~ 1, family = SEP4(), data = DIN,
                         mu.start = mean(DIN$DIN), sigma.start = sd(DIN$DIN),
                         method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
DIN_static_SST <- gamlss(DIN ~ 1, family = SST(), data = DIN,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST1
DIN_static_ST1 <- gamlss(DIN ~ 1, family = ST1(), data = DIN,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
DIN_static_ST2 <- gamlss(DIN ~ 1, family = ST2(), data = DIN,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
DIN_static_ST3 <- gamlss(DIN ~ 1, family = ST3(), data = DIN,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
DIN_static_ST4 <- gamlss(DIN ~ 1, family = ST4(), data = DIN,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
DIN_static_ST5 <- gamlss(DIN ~ 1, family = ST5(), data = DIN,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
DIN_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4" ,"SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(DIN_static_NO, k = 2),GAIC(DIN_static_NO2, k = 2), GAIC(DIN_static_GU, k = 2), GAIC(DIN_static_LO, k = 2), 
          GAIC(DIN_static_RG, k = 2), GAIC(DIN_static_exGAUS, k = 2), GAIC(DIN_static_NOF, k = 2), GAIC(DIN_static_PE, k = 2), 
          GAIC(DIN_static_PE2, k = 2),GAIC(DIN_static_SN1, k = 2), GAIC(DIN_static_TF, k = 2), GAIC(DIN_static_TF2, k = 2),
          GAIC(DIN_static_GT, k = 2),GAIC(DIN_static_JSU, k = 2), GAIC(DIN_static_JSUo, k = 2), GAIC(DIN_static_NET, k = 2),
          GAIC(DIN_static_SHASH, k = 2),GAIC(DIN_static_SHASHo, k = 2), GAIC(DIN_static_SHASHo2, k = 2), GAIC(DIN_static_SEP2, k = 2), 
          GAIC(DIN_static_SEP3, k = 2), GAIC(DIN_static_SEP4, k = 2), GAIC(DIN_static_SST, k = 2), GAIC(DIN_static_ST1, k = 2),
          GAIC(DIN_static_ST2, k = 2),GAIC(DIN_static_ST3, k = 2), GAIC(DIN_static_ST4, k = 2), GAIC(DIN_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4")
)


# Plot
pAIC_4 <- ggplot(DIN_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors, name = "param") +
  labs(x = "AIC", y = "Intercept models (DIN)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



#### SILICATE #### ========================================================================================

# NO
si_static_NO <- gamlss(Silicate ~ 1, family = NO(), data = Silicate) # two-parameter
# NO2
si_static_NO2 <- gamlss(Silicate ~ 1, family = NO2(), data = Silicate)
# GU
si_static_GU <- gamlss(Silicate ~ 1, family = GU(), data = Silicate)
# LO
si_static_LO <- gamlss(Silicate ~ 1, family = LO(), data = Silicate)
# RG
si_static_RG <- gamlss(Silicate ~ 1, family = RG(), data = Silicate)
# exGAUS
si_static_exGAUS <- gamlss(Silicate ~ 1, family = exGAUS(), data = Silicate,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
si_static_NOF <- gamlss(Silicate ~ 1, family = NOF(), data = Silicate)
# PE
si_static_PE <- gamlss(Silicate ~ 1, family = PE(), data = Silicate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
si_static_PE2 <- gamlss(Silicate ~ 1, family = PE2(), data = Silicate, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
si_static_SN1 <- gamlss(Silicate ~ 1, family = SN1(), data = Silicate)
# TF
si_static_TF <- gamlss(Silicate ~ 1, family = TF(), data = Silicate)
# TF2
si_static_TF2 <- gamlss(Silicate ~ 1, family = TF2(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
si_static_GT <- gamlss(Silicate ~ 1, family = GT(), data = Silicate, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
si_static_JSU <- gamlss(Silicate ~ 1, family = JSU(), data = Silicate, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
si_static_JSUo <- gamlss(Silicate ~ 1, family = JSUo(), data = Silicate,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# NET
si_static_NET <- gamlss(Silicate ~ 1, family = NET(), data = Silicate)
# SHASH
si_static_SHASH <- gamlss(Silicate ~ 1, family = SHASH(), data = Silicate,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
si_static_SHASHo <- gamlss(Silicate ~ 1, family = SHASHo(), data = Silicate,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
si_static_SHASHo2 <- gamlss(Silicate ~ 1, family = SHASHo2(), data = Silicate,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
si_static_SEP2 <- gamlss(Silicate ~ 1, family = SEP2(), data = Silicate,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
si_static_SEP3 <- gamlss(Silicate ~ 1, family = SEP3(), data = Silicate,
                          mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate), # convergence issue
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
si_static_SEP4 <- gamlss(Silicate ~ 1, family = SEP4(), data = Silicate,
                          mu.start = mean(Silicate$Silicate), sigma.start = sd(Silicate$Silicate),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
si_static_SST <- gamlss(Silicate ~ 1, family = SST(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # not working 
# ST1
si_static_ST1 <- gamlss(Silicate ~ 1, family = ST1(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
si_static_ST2 <- gamlss(Silicate ~ 1, family = ST2(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
si_static_ST3 <- gamlss(Silicate ~ 1, family = ST3(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
si_static_ST4 <- gamlss(Silicate ~ 1, family = ST4(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
si_static_ST5 <- gamlss(Silicate ~ 1, family = ST5(), data = Silicate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


# scale the AIC 
si_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP4", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(si_static_NO, k = 2),GAIC(si_static_NO2, k = 2), GAIC(si_static_GU, k = 2), GAIC(si_static_LO, k = 2), 
          GAIC(si_static_RG, k = 2), GAIC(si_static_exGAUS, k = 2), GAIC(si_static_NOF, k = 2), GAIC(si_static_PE, k = 2), 
          GAIC(si_static_PE2, k = 2),GAIC(si_static_SN1, k = 2), GAIC(si_static_TF, k = 2), GAIC(si_static_TF2, k = 2),
          GAIC(si_static_GT, k = 2),GAIC(si_static_JSU, k = 2), GAIC(si_static_JSUo, k = 2), GAIC(si_static_NET, k = 2),
          GAIC(si_static_SHASH, k = 2),GAIC(si_static_SHASHo, k = 2), GAIC(si_static_SHASHo2, k = 2), GAIC(si_static_SEP2, k = 2), 
          GAIC(si_static_SEP4, k = 2), GAIC(si_static_ST1, k = 2),
          GAIC(si_static_ST2, k = 2),GAIC(si_static_ST3, k = 2), GAIC(si_static_ST4, k = 2), GAIC(si_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4")
)


# Plot
pAIC_5 <- ggplot(si_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors, name = "param") +
  labs(x = "AIC", y = "Intercept models (Silicate)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



#### PHOSPHATE #### ========================================================================================

# NO
ph_static_NO <- gamlss(Phosphate ~ 1, family = NO(), data = Phosphate) # two-parameter
# NO2
ph_static_NO2 <- gamlss(Phosphate ~ 1, family = NO2(), data = Phosphate)
# GU
ph_static_GU <- gamlss(Phosphate ~ 1, family = GU(), data = Phosphate)
# LO
ph_static_LO <- gamlss(Phosphate ~ 1, family = LO(), data = Phosphate)
# RG
ph_static_RG <- gamlss(Phosphate ~ 1, family = RG(), data = Phosphate)
# exGAUS
ph_static_exGAUS <- gamlss(Phosphate ~ 1, family = exGAUS(), data = Phosphate,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
ph_static_NOF <- gamlss(Phosphate ~ 1, family = NOF(), data = Phosphate)
# PE
ph_static_PE <- gamlss(Phosphate ~ 1, family = PE(), data = Phosphate,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
ph_static_PE2 <- gamlss(Phosphate ~ 1, family = PE2(), data = Phosphate, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
ph_static_SN1 <- gamlss(Phosphate ~ 1, family = SN1(), data = Phosphate)
# TF
ph_static_TF <- gamlss(Phosphate ~ 1, family = TF(), data = Phosphate)
# TF2
ph_static_TF2 <- gamlss(Phosphate ~ 1, family = TF2(), data = Phosphate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
ph_static_GT <- gamlss(Phosphate ~ 1, family = GT(), data = Phosphate, 
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
ph_static_JSU <- gamlss(Phosphate ~ 1, family = JSU(), data = Phosphate, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
ph_static_JSUo <- gamlss(Phosphate ~ 1, family = JSUo(), data = Phosphate,
                         mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue
# NET
ph_static_NET <- gamlss(Phosphate ~ 1, family = NET(), data = Phosphate)
# SHASH
ph_static_SHASH <- gamlss(Phosphate ~ 1, family = SHASH(), data = Phosphate,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
ph_static_SHASHo <- gamlss(Phosphate ~ 1, family = SHASHo(), data = Phosphate,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
ph_static_SHASHo2 <- gamlss(Phosphate ~ 1, family = SHASHo2(), data = Phosphate,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
ph_static_SEP2 <- gamlss(Phosphate ~ 1, family = SEP2(), data = Phosphate,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
ph_static_SEP3 <- gamlss(Phosphate ~ 1, family = SEP3(), data = Phosphate,
                         mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate), 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
ph_static_SEP4 <- gamlss(Phosphate ~ 1, family = SEP4(), data = Phosphate,
                         mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
ph_static_SST <- gamlss(Phosphate ~ 1, family = SST(), data = Phosphate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # not working 
# ST1
ph_static_ST1 <- gamlss(Phosphate ~ 1, family = ST1(), data = Phosphate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
ph_static_ST2 <- gamlss(Phosphate ~ 1, family = ST2(), data = Phosphate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
ph_static_ST3 <- gamlss(Phosphate ~ 1, family = ST3(), data = Phosphate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
ph_static_ST4 <- gamlss(Phosphate ~ 1, family = ST4(), data = Phosphate,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
ph_static_ST5 <- gamlss(Phosphate ~ 1, family = ST5(), data = Phosphate,
                        mu.start = mean(Phosphate$Phosphate), sigma.start = sd(Phosphate$Phosphate),
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # convergence issue


# scale the AIC 
ph_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS", "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "ST1", "ST2", "ST3", "ST4"),
  AIC = c(GAIC(ph_static_NO, k = 2),GAIC(ph_static_NO2, k = 2), GAIC(ph_static_GU, k = 2), GAIC(ph_static_LO, k = 2), 
          GAIC(ph_static_RG, k = 2), GAIC(ph_static_exGAUS, k = 2), GAIC(ph_static_NOF, k = 2), GAIC(ph_static_PE, k = 2), 
          GAIC(ph_static_PE2, k = 2),GAIC(ph_static_SN1, k = 2), GAIC(ph_static_TF, k = 2), GAIC(ph_static_TF2, k = 2),
          GAIC(ph_static_GT, k = 2),GAIC(ph_static_JSU, k = 2), GAIC(ph_static_NET, k = 2),
          GAIC(ph_static_SHASH, k = 2),GAIC(ph_static_SHASHo, k = 2), GAIC(ph_static_SHASHo2, k = 2), GAIC(ph_static_SEP2, k = 2), 
          GAIC(ph_static_SEP3, k = 2), GAIC(ph_static_SEP4, k = 2), GAIC(ph_static_ST1, k = 2),
          GAIC(ph_static_ST2, k = 2),GAIC(ph_static_ST3, k = 2), GAIC(ph_static_ST4, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4")
)


# Plot
pAIC_6 <- ggplot(ph_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = custom_colors, name = "param") +
  labs(x = "AIC", y = "Intercept models (Phosphate)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# all AIC plots together
pAIC_1 + pAIC_2 + pAIC_3 + pAIC_4 + pAIC_5 + pAIC_6 


#############################################################################################

# Combine the plots the three plots together for each nutrient
# Nitrate
(p <- (pp1 | (p1 / pAIC_1)) + 
    plot_layout(widths = c(3, 1)))   


#############################################################################################

#####################################################
######## OVERLAP MODEL PDF WITH DATA HIST ###########
#####################################################

# Function to simulate from a fitted gamlss model
sim_from_model <- function(model, n = 1000) {
  fam <- family(model)[[1]]   # distribution family object
  rfun <- get(paste0("r", fam)) 
  mu    <- predict(model, "mu", type = "response")[1]
  sigma <- fitted(model, "sigma", type = "response")[1]
  nu    <- tryCatch(predict(model, "nu", type = "response")[1], error = function(e) NULL)
  tau   <- tryCatch(predict(model, "tau", type = "response")[1], error = function(e) NULL)
  
  # match args by what the family uses
  args <- list(n = n, mu = mu, sigma = sigma)
  if (!is.null(nu))  args$nu  <- nu
  if (!is.null(tau)) args$tau <- tau
  
  tibble(value = do.call(rfun, args))
}


# Nitrate best models for each n-parameter distribution
models_ni <- list(SEP2 = ni_static_SEP2, exGAUS = ni_static_exGAUS, RG = ni_static_RG)
# Simulate and bind results
set.seed(123)
sims <- bind_rows(
  lapply(names(models_ni), function(name) {
    sim_from_model(models_ni[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)
# Plot overlapping densities
# Nitrate
ggplot(Nitrate, aes(x = Nitrate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 10) + 
  xlab("Nitrate (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(data = sims, aes(x = value, color = model)) +
  geom_density(colour = "black", alpha = 0.5) +
  theme_bw() 


# Nitrite best models for each n-parameter distribution
models_nii <- list(JSU = nii_static_JSU, exGAUS = nii_static_exGAUS, RG = nii_static_RG)
# Simulate and bind results
set.seed(123)
sims <- bind_rows(
  lapply(names(models_nii), function(name) {
    sim_from_model(models_nii[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)
# Plot overlapping densities
ggplot(Nitrite, aes(x = Nitrite)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  xlab("Nitrite (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,5)) +
  geom_density(data = sims, aes(x = value, color = model)) +
  geom_density(colour = "black", alpha = 0.5) +
  theme_bw() 


# Ammonium best models for each n-parameter distribution
models_am <- list(SEP2 = am_static_SEP2, exGAUS = am_static_exGAUS, RG = am_static_RG)
# Simulate and bind results
set.seed(123)
sims <- bind_rows(
  lapply(names(models_am), function(name) {
    sim_from_model(models_am[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)
# Plot overlapping densities
ggplot(Ammonium, aes(x = Ammonium)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("Ammonium (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,35)) +
  geom_density(data = sims, aes(x = value, color = model)) +
  geom_density(colour = "black", alpha = 0.5) +
  theme_bw() 


# DIN best models for each n-parameter distribution
models_DIN <- list(JSU = DIN_static_JSU, exGAUS = DIN_static_exGAUS, RG = DIN_static_RG)
# Simulate and bind results
set.seed(123)
sims <- bind_rows(
  lapply(names(models_DIN), function(name) {
    sim_from_model(models_DIN[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)
# Plot overlapping densities
ggplot(DIN, aes(x = DIN)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 10) + 
  xlab("DIN (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(data = sims, aes(x = value, color = model)) +
  geom_density(colour = "black", alpha = 0.5) +
  theme_bw() 


# Silicate best models for each n-parameter distribution
models_si <- list(SEP2 = si_static_SEP2, exGAUS = si_static_exGAUS, RG = si_static_RG)
# Simulate and bind results
set.seed(123)
sims <- bind_rows(
  lapply(names(models_si), function(name) {
    sim_from_model(models_si[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)
# Plot overlapping densities
ggplot(Silicate, aes(x = Silicate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 5) + 
  xlab("Silicate (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(data = sims, aes(x = value, color = model)) +
  geom_density(colour = "black", alpha = 0.5) +
  theme_bw() 


# Phosphate best models for each n-parameter distribution
models_ph <- list(SEP4 = ph_static_SEP4, PE = ph_static_PE, NO = ph_static_NO)
# Simulate and bind results
set.seed(123)
sims <- bind_rows(
  lapply(names(models_ph), function(name) {
    sim_from_model(models_ph[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)
# Plot overlapping densities
ggplot(Phosphate, aes(x = Phosphate)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.2) + 
  xlab("Phosphate (µmol/l)") + 
  ylab("Density") +
  xlim(x = c(0,2.5)) +
  geom_density(data = sims, aes(x = value, color = model)) +
  geom_density(colour = "black", alpha = 0.5) +
  theme_bw() 





