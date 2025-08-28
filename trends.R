H62 <- read.csv("dataset/Helgoland_1962.csv")
H63 <- read.csv("dataset/Helgoland_1963.csv")
H64 <- read.csv("dataset/Helgoland_1964.csv")
H65 <- read.csv("dataset/Helgoland_1965.csv")
H66 <- read.csv("dataset/Helgoland_1966.csv")
H67 <- read.csv("dataset/Helgoland_1967.csv")
H68 <- read.csv("dataset/Helgoland_1968.csv")
H69 <- read.csv("dataset/Helgoland_1969.csv")
H70 <- read.csv("dataset/Helgoland_1970.csv")
H71 <- read.csv("dataset/Helgoland_1971.csv")
H72 <- read.csv("dataset/Helgoland_1972.csv")
H73 <- read.csv("dataset/Helgoland_1973.csv")
H74 <- read.csv("dataset/Helgoland_1974.csv")
H75 <- read.csv("dataset/Helgoland_1975.csv")
H76 <- read.csv("dataset/Helgoland_1976.csv")
H77 <- read.csv("dataset/Helgoland_1977.csv")
H78 <- read.csv("dataset/Helgoland_1978.csv")
H79 <- read.csv("dataset/Helgoland_1979.csv")
H80 <- read.csv("dataset/Helgoland_1980.csv")
H81 <- read.csv("dataset/Helgoland_1981.csv")
H82 <- read.csv("dataset/Helgoland_1982.csv")
H83 <- read.csv("dataset/Helgoland_1983.csv")
H84 <- read.csv("dataset/Helgoland_1984.csv")
H85 <- read.csv("dataset/Helgoland_1985.csv")
H86 <- read.csv("dataset/Helgoland_1986.csv")
H87 <- read.csv("dataset/Helgoland_1987.csv")
H88 <- read.csv("dataset/Helgoland_1988.csv")
H89 <- read.csv("dataset/Helgoland_1989.csv")
H90 <- read.csv("dataset/Helgoland_1990.csv")
H91 <- read.csv("dataset/Helgoland_1991.csv")
H92 <- read.csv("dataset/Helgoland_1992.csv")
H93 <- read.csv("dataset/Helgoland_1993.csv")
H94 <- read.csv("dataset/Helgoland_1994.csv")
H98 <- read.csv("dataset/Helgoland_1998.csv")
H2002 <- read.csv("dataset/Helgoland_2002.csv")

library(dplyr)
library(ggplot2)
library(patchwork)

# Combine all datasets into one
df <- bind_rows(H62, H63, H64, H65, H66, H67, H68, H69, H70, H71, H72, H73, H74, H75, H76,
                H77, H78, H79, H80, H81, H82, H83, H84, H85, H86, H87, H88, H89, H90, H91,
                H92, H93, H94)

library(lubridate)

# add year as a column
df$year <- year(ymd_hm(df$Date))
df$Date <- as.Date(df$Date)


##############
### TRENDS ###
##############

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ---- Setup ----
# Replace `df` with your data frame name
# Edit this list to match your nutrient column names:
nutrients <- c("Nitrate", "Nitrite", "Silicate", "DIN") # add more if you have them
nitrate <- "Nitrate"
DIN <- "DIN"
phosphate <- "Phosphate"
phyto_col <- "Phytopl"
phyto_biom <- "Phytopl_C"


# If you already have a 'year' column instead of 'date', replace the mutate() with:
df_year <- df %>%
  group_by(year) %>%
  summarise(across(all_of(c(phosphate, phyto_biom)), ~ mean(., na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-year, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup()

ggplot(df_year, aes(x = year, y = z, color = variable)) +
 # geom_line(linewidth = 0.7, alpha = 0.9) +
  geom_smooth(se = FALSE, method = "loess", span = 0.3, linewidth = 0.6) +
  labs(x = "Year", y = "Standardized (z-score)",
       color = "Variable",
       title = "Yearly trends (standardized) of nutrients and phytoplankton") +
  theme_minimal()



### Looking at trends per year ###

target_year <- 1968

# If you have a Date column named `date`:
df_one <- df %>%
  mutate(year = year(Date), month = month(Date, label = TRUE)) %>%
  filter(year == target_year) %>%
  group_by(month) %>%
  summarise(across(all_of(c(phosphate, phyto_col)), ~mean(., na.rm = TRUE)),
            .groups = "drop")

# (If you already have a `year` column, drop mutate(year=...) and filter(year == target_year);
# add/create a `month` column if you have dates.)

df_long <- df_one %>%
  pivot_longer(-month, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  mutate(z = as.numeric(scale(value))) %>%
  ungroup()

ggplot(df_long, aes(month, z, color = variable, group = variable)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  labs(title = paste("Monthly trends (standardized) in", target_year),
       x = "Month", y = "Z-score", color = "Variable") +
  theme_minimal()


