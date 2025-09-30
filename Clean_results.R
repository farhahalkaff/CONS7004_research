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
  theme_minimal() +
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
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



