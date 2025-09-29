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

