## Simple simulation to prove first hypothesis; 
## 1. H1: heavy-tailed distribution is common in ecological time-series
## now however, I'm adding platy and negative skew 


## SIMULATION WITH DIFFERENCE SCENARIOS ##
#--- sim 10: -ve skewness to +ve skewness,  normal tails (mesokurtic)
#--- sim 11: -ve skewness to +ve skewness,  heavy tails (leptokurtic)
#--- sim 12: +ve skewness, light to heavy tails 
#--- sim 13: -ve skewness to +ve skewness, light to heavy tails 


# now I'm adding time-varying factor, where the tails are not consistently heavy, but can be light and 
# normal at times. similar with skewness, it can be negative, postive and symetrical at times 


library(knitr)
library(gamlss)
library(gamlss.dist)
library(ggplot2)
library(dplyr)
library(ggpubr)



##### SIM 10: complex, -ve skewness to +ve skewness,  normal tails --> SHASHo() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# changing through time
beta_nu <- 0.009
nu_t <- -0.5   + beta_nu * t_index2

# observation randomly taken from SST
y10 <- rSHASHo(t2, mu = mu0, sigma = sigma0, nu = nu_t, tau = 1)

# create dataset
dat10 <- data.frame(
  t = t_index2,
  y = y10,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = nu_t,
  tau_true = 1,
  dist = "SHASHo"
)

# look at density plot
ggplot(dat10, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat10$y), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim10_static_NO <- gamlss(y ~ 1, family = NO(), data = dat10) # two-parameter
# NO2
sim10_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat10)
# GU
sim10_static_GU <- gamlss(y ~ 1, family = GU(), data = dat10)
# LO
sim10_static_LO <- gamlss(y ~ 1, family = LO(), data = dat10)
# RG
sim10_static_RG <- gamlss(y ~ 1, family = RG(), data = dat10)
# exGAUS
sim10_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat10, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim10_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat10,
                          method = mixed(),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim10_static_PE <- gamlss(y ~ 1, family = PE(), data = dat10,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim10_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat10, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim10_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat10) # got a warning NANs produced
# SN2
sim10_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat10,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim10_static_TF <- gamlss(y ~ 1, family = TF(), data = dat10)
# TF2
sim10_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat10,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim10_static_GT <- gamlss(y ~ 1, family = GT(), data = dat10, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim10_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat10, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim10_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat10,
                           mu.start = mean(dat10$y), sigma.start = sd(dat10$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim10_static_NET <- gamlss(y ~ 1, family = NET(), data = dat10)
# SHASH
sim10_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat10,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim10_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat10,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim10_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat10,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim10_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat10,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim10_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat10,
                           mu.start = mean(dat10$y), sigma.start = sd(dat10$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim10_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat10,
                           mu.start = mean(dat10$y), sigma.start = sd(dat10$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim10_static_SST <- gamlss(y ~ 1, family = SST(), data = dat10, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim10_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat10,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim10_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat10,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim10_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat10,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim10_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat10,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim10_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat10,
                          mu.start = mean(dat10$y), sigma.start = sd(dat10$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim10_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "SST","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim10_static_NO, k = 2),GAIC(sim10_static_NO2, k = 2), GAIC(sim10_static_GU, k = 2), GAIC(sim10_static_LO, k = 2), 
          GAIC(sim10_static_RG, k = 2), GAIC(sim10_static_exGAUS, k = 2), GAIC(sim10_static_NOF, k = 2), GAIC(sim10_static_PE, k = 2), 
          GAIC(sim10_static_PE2, k = 2),GAIC(sim10_static_SN1, k = 2), GAIC(sim10_static_SN2, k = 2),GAIC(sim10_static_TF, k = 2), GAIC(sim10_static_TF2, k = 2),
          GAIC(sim10_static_GT, k = 2),GAIC(sim10_static_JSU, k = 2), GAIC(sim10_static_JSUo, k = 2), GAIC(sim10_static_NET, k = 2),
          GAIC(sim10_static_SHASH, k = 2),GAIC(sim10_static_SHASHo, k = 2), GAIC(sim10_static_SHASHo2, k = 2), GAIC(sim10_static_SEP2, k = 2), 
          GAIC(sim10_static_SEP3, k = 2), GAIC(sim10_static_SEP4, k = 2), GAIC(sim10_static_SST, k = 2), GAIC(sim10_static_ST1, k = 2),
          GAIC(sim10_static_ST2, k = 2),GAIC(sim10_static_ST3, k = 2), GAIC(sim10_static_ST4, k = 2), GAIC(sim10_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both", "both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto", "lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim10_pAIC_1a <- ggplot(sim10_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(x = "AIC", y = "models") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


# param recovery ------------------------------------------------------------------
param_table10 <- data.frame(
  models = c("true","SHASHo", "SN2"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10",round(sim10_static_SHASHo$mu.coefficients,3), round(sim10_static_SN2$mu.coefficients,3)),
  sigma = c("2",round(exp(sim10_static_SHASHo$sigma.coefficients),3), round(exp(sim10_static_SN2$sigma.coefficients),3)),
  nu = c("-0.5",round(sim10_static_SHASHo$nu.coefficients[1],3), round(exp(sim10_static_SN2$nu.coefficients[1]),3)),
  tau = c("1",round(exp(sim10_static_SHASHo$tau.coefficients),3), "-")
)


# compare the two models PDF -----------------------------------------------------
model_compare10 <- list(SHASHo = sim10_static_SHASHo, SN2 = sim10_static_SN2)

# Function to simulate from a fitted gamlss model
sim_from_model <- function(model, n = 1000) {
  fam <- family(model)[[1]]   # distribution family object
  rfun <- get(paste0("r", fam)) 
  mu    <- predict(model, "mu", type = "response")[1]
  sigma <- predict(model, "sigma", type = "response")[1]
  nu    <- tryCatch(predict(model, "nu", type = "response")[1], error = function(e) NULL)
  tau   <- tryCatch(predict(model, "tau", type = "response")[1], error = function(e) NULL)
  
  # match args by what the family uses
  args <- list(n = n, mu = mu, sigma = sigma)
  if (!is.null(nu))  args$nu  <- nu
  if (!is.null(tau)) args$tau <- tau
  
  tibble(value = do.call(rfun, args))
}

# Simulate and bind results
set.seed(123)
sim10_sim <- bind_rows(
  lapply(names(model_compare10), function(name) {
    sim_from_model(model_compare10[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim10_p1 <- ggplot(dat10, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 3) + 
  labs(title = "Sim 5: -ve to +ve skew, normal tails", x = "Value", y = "Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim10_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table10 <- ggtexttable(param_table10, rows = NULL, theme = ttheme("light", base_size = 8.5))
table10 <- table_cell_bg(table10, row = 2, column = 1:5, linewidth = 5,
                        fill="grey", color = "grey")
# combeanation
sim10_plots <- (sim10_pAIC_1a  +  sim10_p1 / table10)




##### SIM 11: complex, -ve skewness to +ve skewness,  heavy tails --> SHASHo() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# changing through time
beta_nu <- 0.009
nu_t <- -0.5   + beta_nu * t_index2

# observation randomly taken from SST
y11 <- rSHASHo(t2, mu = mu0, sigma = sigma0, nu = nu_t, tau = 0.6)

# create dataset
dat11 <- data.frame(
  t = t_index2,
  y = y11,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = nu_t,
  tau_true = 0.2,
  dist = "SHASHo"
)

ggplot(dat11, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 

# look at density plot
ggplot(dat11, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat11$y), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim11_static_NO <- gamlss(y ~ 1, family = NO(), data = dat11) # two-parameter
# NO2
sim11_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat11)
# GU
sim11_static_GU <- gamlss(y ~ 1, family = GU(), data = dat11)
# LO
sim11_static_LO <- gamlss(y ~ 1, family = LO(), data = dat11)
# RG
sim11_static_RG <- gamlss(y ~ 1, family = RG(), data = dat11)
# exGAUS
sim11_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat11, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim11_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat11,
                           method = mixed(),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim11_static_PE <- gamlss(y ~ 1, family = PE(), data = dat11,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim11_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat11, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim11_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat11) # got a warning NANs produced
# SN2
sim11_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat11,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim11_static_TF <- gamlss(y ~ 1, family = TF(), data = dat11)
# TF2
sim11_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat11,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim11_static_GT <- gamlss(y ~ 1, family = GT(), data = dat11, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim11_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat11, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim11_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat11,
                            mu.start = mean(dat11$y), sigma.start = sd(dat11$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim11_static_NET <- gamlss(y ~ 1, family = NET(), data = dat11)
# SHASH
sim11_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat11,
                            method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim11_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat11,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim11_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat11,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim11_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat11,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim11_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat11,
                            mu.start = mean(dat11$y), sigma.start = sd(dat11$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim11_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat11, # not working
                            mu.start = mean(dat11$y), sigma.start = sd(dat11$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim11_static_SST <- gamlss(y ~ 1, family = SST(), data = dat11, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim11_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat11,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim11_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat11,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim11_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat11,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim11_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat11,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim11_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat11,
                           mu.start = mean(dat11$y), sigma.start = sd(dat11$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim11_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "SST","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim11_static_NO, k = 2),GAIC(sim11_static_NO2, k = 2), GAIC(sim11_static_GU, k = 2), GAIC(sim11_static_LO, k = 2), 
          GAIC(sim11_static_RG, k = 2), GAIC(sim11_static_exGAUS, k = 2), GAIC(sim11_static_NOF, k = 2), GAIC(sim11_static_PE, k = 2), 
          GAIC(sim11_static_PE2, k = 2),GAIC(sim11_static_SN1, k = 2), GAIC(sim11_static_SN2, k = 2),GAIC(sim11_static_TF, k = 2), GAIC(sim11_static_TF2, k = 2),
          GAIC(sim11_static_GT, k = 2),GAIC(sim11_static_JSU, k = 2), GAIC(sim11_static_JSUo, k = 2), GAIC(sim11_static_NET, k = 2),
          GAIC(sim11_static_SHASH, k = 2),GAIC(sim11_static_SHASHo, k = 2), GAIC(sim11_static_SHASHo2, k = 2), GAIC(sim11_static_SEP2, k = 2), 
          GAIC(sim11_static_SEP3, k = 2),  GAIC(sim11_static_SST, k = 2), GAIC(sim11_static_ST1, k = 2),
          GAIC(sim11_static_ST2, k = 2),GAIC(sim11_static_ST3, k = 2), GAIC(sim11_static_ST4, k = 2), GAIC(sim11_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both", "both", "both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","lepto" ,"lepto", "lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim11_pAIC_1a <- ggplot(sim11_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(x = "AIC", y = "models") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


# param recovery ------------------------------------------------------------------
param_table11 <- data.frame(
  models = c("true","JSUo", "SHASHo2"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10", round(sim11_static_JSUo$mu.coefficients,3), round(sim11_static_SHASHo2$mu.coefficients,3)),
  sigma = c("2", round(exp(sim11_static_JSUo$sigma.coefficients),3), round(exp(sim11_static_SHASHo2$sigma.coefficients),3)),
  nu = c("-0.5", round(sim11_static_JSUo$nu.coefficients[1],3),round(sim11_static_SHASHo2$nu.coefficients[1],3)),
  tau = c("0.2",round(exp(sim11_static_JSUo$tau.coefficients),3), round(exp(sim11_static_SHASHo2$tau.coefficients),3))
)


# compare the two models PDF -----------------------------------------------------
model_compare11 <- list(JSUo = sim11_static_JSUo, SHASHo2 = sim11_static_SHASHo2)

# Simulate and bind results
set.seed(123)
sim11_sim <- bind_rows(
  lapply(names(model_compare11), function(name) {
    sim_from_model(model_compare11[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim11_p1 <- ggplot(dat11, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 10) + 
  labs(title = "Sim 6: -ve to +ve skew, heavy tails", x = "Value", y = "Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim11_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table11 <- ggtexttable(param_table11, rows = NULL, theme = ttheme("light", base_size = 8.5))
table11 <- table_cell_bg(table11, row = 2, column = 1:5, linewidth = 5,
                         fill="grey", color = "grey")
# combeanation
sim11_plots <- (sim11_pAIC_1a  +  sim11_p1 / table11)




##### SIM 12: complex, +ve skewness,  light to heavy tails --> SHASHo() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# changing through time
beta_tau <- -0.0025
tau_t <- 1.5   + beta_tau * t_index2

# observation randomly taken from SST
y12 <- rSHASHo(t2, mu = mu0, sigma = sigma0, nu = 2, tau = tau_t)

# create dataset
dat12 <- data.frame(
  t = t_index2,
  y = y12,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 2,
  tau_true = tau_t,
  dist = "SHASHo"
)

ggplot(dat12, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 

# look at density plot
ggplot(dat12, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat12$y), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim12_static_NO <- gamlss(y ~ 1, family = NO(), data = dat12) # two-parameter
# NO2
sim12_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat12)
# GU
sim12_static_GU <- gamlss(y ~ 1, family = GU(), data = dat12)
# LO
sim12_static_LO <- gamlss(y ~ 1, family = LO(), data = dat12)
# RG
sim12_static_RG <- gamlss(y ~ 1, family = RG(), data = dat12)
# exGAUS
sim12_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat12, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim12_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat12,
                           method = mixed(),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim12_static_PE <- gamlss(y ~ 1, family = PE(), data = dat12,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim12_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat12, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim12_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat12) 
# SN2
sim12_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat12,  # convergence issue
                           mu.start = mean(dat12$y), sigma.start = sd(dat12$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim12_static_TF <- gamlss(y ~ 1, family = TF(), data = dat12)
# TF2
sim12_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat12,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim12_static_GT <- gamlss(y ~ 1, family = GT(), data = dat12, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim12_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat12, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim12_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat12,
                            mu.start = mean(dat12$y), sigma.start = sd(dat12$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim12_static_NET <- gamlss(y ~ 1, family = NET(), data = dat12)
# SHASH
sim12_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat12,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim12_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat12,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim12_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat12,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim12_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat12,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim12_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat12,
                            mu.start = mean(dat12$y), sigma.start = sd(dat12$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim12_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat12, # not working
                            mu.start = mean(dat12$y), sigma.start = sd(dat12$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim12_static_SST <- gamlss(y ~ 1, family = SST(), data = dat12,  # not working
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim12_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat12,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim12_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat12,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim12_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat12,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim12_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat12,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim12_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat12,
                           mu.start = mean(dat12$y), sigma.start = sd(dat12$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim12_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim12_static_NO, k = 2),GAIC(sim12_static_NO2, k = 2), GAIC(sim12_static_GU, k = 2), GAIC(sim12_static_LO, k = 2), 
          GAIC(sim12_static_RG, k = 2), GAIC(sim12_static_exGAUS, k = 2), GAIC(sim12_static_NOF, k = 2), GAIC(sim12_static_PE, k = 2), 
          GAIC(sim12_static_PE2, k = 2),GAIC(sim12_static_SN1, k = 2), GAIC(sim12_static_TF, k = 2), GAIC(sim12_static_TF2, k = 2),
          GAIC(sim12_static_GT, k = 2),GAIC(sim12_static_JSU, k = 2), GAIC(sim12_static_JSUo, k = 2), GAIC(sim12_static_NET, k = 2),
          GAIC(sim12_static_SHASH, k = 2),GAIC(sim12_static_SHASHo, k = 2), GAIC(sim12_static_SHASHo2, k = 2), GAIC(sim12_static_SEP2, k = 2), 
          GAIC(sim12_static_SEP3, k = 2), GAIC(sim12_static_ST1, k = 2),
          GAIC(sim12_static_ST2, k = 2),GAIC(sim12_static_ST3, k = 2), GAIC(sim12_static_ST4, k = 2), GAIC(sim12_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both", "both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both" ,"lepto", "lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim12_pAIC_1a <- ggplot(sim12_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(x = "AIC", y = "models") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


# param recovery ------------------------------------------------------------------
param_table12 <- data.frame(
  models = c("true","JSUo", "SEP3"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10", round(sim12_static_JSUo$mu.coefficients,3), round(sim12_static_SEP3$mu.coefficients,3)),
  sigma = c("2", round(exp(sim12_static_JSUo$sigma.coefficients),3), round(exp(sim12_static_SEP3$sigma.coefficients),3)),
  nu = c("2", round(sim12_static_JSUo$nu.coefficients[1],3),round(exp(sim12_static_SEP3$nu.coefficients[1]),3)),
  tau = c("1.5",round(exp(sim12_static_JSUo$tau.coefficients),3), round(exp(sim12_static_SEP3$tau.coefficients),3))
)


# compare the two models PDF -----------------------------------------------------
model_compare12 <- list(JSUo = sim12_static_JSUo, SEP3 = sim12_static_SEP3)

# Simulate and bind results
set.seed(123)
sim12_sim <- bind_rows(
  lapply(names(model_compare12), function(name) {
    sim_from_model(model_compare12[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim12_p1 <- ggplot(dat12, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 5) + 
  labs(title = "Sim 7: +ve skew, light yo heavy tails", x = "Value", y = "Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim12_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table12 <- ggtexttable(param_table12, rows = NULL, theme = ttheme("light", base_size = 8.5))
table12 <- table_cell_bg(table12, row = 2, column = 1:5, linewidth = 5,
                         fill="grey", color = "grey")
# combeanation
sim12_plots <- (sim12_pAIC_1a  +  sim12_p1 / table12)




##### SIM 13: complex, -ve to +ve skewness,  light to heavy tails --> SHASHo() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# changing through time
# tails
beta_tau <- -0.0025
tau_t <- 1.5   + beta_tau * t_index2
# skew
beta_nu <- 0.009
nu_t <- -0.5   + beta_nu * t_index2

# observation randomly taken from SST
y13 <- rSHASHo(t2, mu = mu0, sigma = sigma0, nu = nu_t, tau = tau_t)

# create dataset
dat13 <- data.frame(
  t = t_index2,
  y = y13,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = nu_t,
  tau_true = tau_t,
  dist = "SHASHo"
)

ggplot(dat13, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 

# look at density plot
ggplot(dat13, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat13$y), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim13_static_NO <- gamlss(y ~ 1, family = NO(), data = dat13) # two-parameter
# NO2
sim13_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat13)
# GU
sim13_static_GU <- gamlss(y ~ 1, family = GU(), data = dat13)
# LO
sim13_static_LO <- gamlss(y ~ 1, family = LO(), data = dat13)
# RG
sim13_static_RG <- gamlss(y ~ 1, family = RG(), data = dat13)
# exGAUS
sim13_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat13, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim13_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat13,
                           method = mixed(),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim13_static_PE <- gamlss(y ~ 1, family = PE(), data = dat13,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim13_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat13, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim13_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat13) 
# SN2
sim13_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat13,  
                           mu.start = mean(dat13$y), sigma.start = sd(dat13$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim13_static_TF <- gamlss(y ~ 1, family = TF(), data = dat13)
# TF2
sim13_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat13,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim13_static_GT <- gamlss(y ~ 1, family = GT(), data = dat13, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim13_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat13, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim13_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat13,
                            mu.start = mean(dat13$y), sigma.start = sd(dat13$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim13_static_NET <- gamlss(y ~ 1, family = NET(), data = dat13)
# SHASH
sim13_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat13,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim13_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat13,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim13_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat13,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim13_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat13,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim13_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat13,
                            mu.start = mean(dat13$y), sigma.start = sd(dat13$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim13_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat13,  # convergnce issue
                            mu.start = mean(dat13$y), sigma.start = sd(dat13$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim13_static_SST <- gamlss(y ~ 1, family = SST(), data = dat13,  
                           mu.start = mean(dat13$y), sigma.start = sd(dat13$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim13_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat13,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim13_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat13,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim13_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat13,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim13_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat13,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim13_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat13,
                           mu.start = mean(dat13$y), sigma.start = sd(dat13$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim13_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "SST","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim13_static_NO, k = 2),GAIC(sim13_static_NO2, k = 2), GAIC(sim13_static_GU, k = 2), GAIC(sim13_static_LO, k = 2), 
          GAIC(sim13_static_RG, k = 2), GAIC(sim13_static_exGAUS, k = 2), GAIC(sim13_static_NOF, k = 2), GAIC(sim13_static_PE, k = 2), 
          GAIC(sim13_static_PE2, k = 2),GAIC(sim13_static_SN1, k = 2), GAIC(sim13_static_SN2, k = 2),GAIC(sim13_static_TF, k = 2), GAIC(sim13_static_TF2, k = 2),
          GAIC(sim13_static_GT, k = 2),GAIC(sim13_static_JSU, k = 2), GAIC(sim13_static_JSUo, k = 2), GAIC(sim13_static_NET, k = 2),
          GAIC(sim13_static_SHASH, k = 2),GAIC(sim13_static_SHASHo, k = 2), GAIC(sim13_static_SHASHo2, k = 2), GAIC(sim13_static_SEP2, k = 2), 
          GAIC(sim13_static_SEP3, k = 2), GAIC(sim13_static_SST, k = 2), GAIC(sim13_static_ST1, k = 2),
          GAIC(sim13_static_ST2, k = 2),GAIC(sim13_static_ST3, k = 2), GAIC(sim13_static_ST4, k = 2), GAIC(sim13_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both", "both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","lepto" ,"lepto", "lepto", "lepto", "lepto", "lepto")
)



# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim13_pAIC_1a <- ggplot(sim13_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(x = "AIC", y = "models") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


# param recovery ------------------------------------------------------------------
param_table13 <- data.frame(
  models = c("true","SHASH"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10", round(sim13_static_SHASH$mu.coefficients,3)),
  sigma = c("2", round(exp(sim13_static_SHASH$sigma.coefficients),3)),
  nu = c("-0.5", round(sim13_static_SHASH$nu.coefficients[1],3)),
  tau = c("1.5",round(exp(sim12_static_SHASH$tau.coefficients),3))
)


# compare the two models PDF -----------------------------------------------------
model_compare13 <- list(SHASH = sim13_static_SHASH)

# Simulate and bind results
set.seed(123)
sim13_sim <- bind_rows(
  lapply(names(model_compare13), function(name) {
    sim_from_model(model_compare13[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim13_p1 <- ggplot(dat13, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  labs(title = "Sim 8: -ve to +ve skew, light to heavy tails", x = "Value", y = "Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim12_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table13 <- ggtexttable(param_table13, rows = NULL, theme = ttheme("light", base_size = 8.5))
table13 <- table_cell_bg(table13, row = 2, column = 1:5, linewidth = 5,
                         fill="grey", color = "grey")
# combeanation
sim13_plots <- (sim13_pAIC_1a  +  sim13_p1 / table13)

