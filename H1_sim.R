## Simple simulation to prove first hypothesis; 
## 1. H1: heavy-tailed distribution is common in ecological time-series

## SIMULATION WITH DIFFERENCE SCENARIOS ##
#--- sim 6: symmetrical,  normal tails (mesokurtic)
#--- sim 7: +ve skewness, normal tails (mesokurtic)
#--- sim 8: symmetrical,  heavy tails (leptokurtic)
#--- sim 9: +ve skewness, heavy tails (leptokurtic)


## we want to prove that flexible family distribution does not win all the time. 
## data that are symmetrical and normal tails, we would assume normal distribution will do best,
## data that are +ve skewed and normal tails, we would assume family distribution that assumes normal 
## tails will do best.
## If flexible distribution (can model both +ve skew and -ve skew, and lepto and platy) does win, that
## the results of our first hypothesis on empirical data would be inconclusive. We won't know if
## the data are actually heavy-tailed, because the flexible distribution has the lowest AIC, and yet
## even if we simulate data to me symmetrical and normal tails, flexible distribution also wins, therefore
## we can't know for sure that our data is in fact heavy tailed 

library(knitr)
library(gamlss)
library(gamlss.dist)
library(ggplot2)
library(dplyr)
library(ggpubr)



##### SIM 6: simple, sym, meso --> SST() ==================================================================


t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# observation randomly taken from SST
set.seed(123)
y6 <- rSST(t2, mu = mu0, sigma = sigma0, nu = 1, tau = 1000)

# create dataset
dat6 <- data.frame(
  t = t_index2,
  y = y6,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 1,
  tau_true = 1000,
  dist = "SST"
)

# look at density plot
ggplot(dat6, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat6$y), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim6_static_NO <- gamlss(y ~ 1, family = NO(), data = dat6) # two-parameter
# NO2
sim6_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat6)
# GU
sim6_static_GU <- gamlss(y ~ 1, family = GU(), data = dat6)
# LO
sim6_static_LO <- gamlss(y ~ 1, family = LO(), data = dat6)
# RG
sim6_static_RG <- gamlss(y ~ 1, family = RG(), data = dat6)
# exGAUS
sim6_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat6, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim6_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat6,
                          method = mixed(),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim6_static_PE <- gamlss(y ~ 1, family = PE(), data = dat6,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim6_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat6, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim6_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat6)
# SN2
sim6_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat6,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim6_static_TF <- gamlss(y ~ 1, family = TF(), data = dat6)
# TF2
sim6_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat6,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim6_static_GT <- gamlss(y ~ 1, family = GT(), data = dat6, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim6_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat6, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim6_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat6,
                           mu.start = mean(dat6$y), sigma.start = sd(dat6$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim6_static_NET <- gamlss(y ~ 1, family = NET(), data = dat6)
# SHASH
sim6_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat6,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim6_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat6,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim6_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat6,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim6_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat6,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim6_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat6,
                           mu.start = mean(dat6$y), sigma.start = sd(dat6$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim6_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat6,
                           mu.start = mean(dat6$y), sigma.start = sd(dat6$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim6_static_SST <- gamlss(y ~ 1, family = SST(), data = dat6, # not working
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim6_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat6,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim6_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat6,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim6_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat6,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim6_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat6,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim6_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat6,
                          mu.start = mean(dat6$y), sigma.start = sd(dat6$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim6_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim6_static_NO, k = 2),GAIC(sim6_static_NO2, k = 2), GAIC(sim6_static_GU, k = 2), GAIC(sim6_static_LO, k = 2), 
          GAIC(sim6_static_RG, k = 2), GAIC(sim6_static_exGAUS, k = 2), GAIC(sim6_static_NOF, k = 2), GAIC(sim6_static_PE, k = 2), 
          GAIC(sim6_static_PE2, k = 2),GAIC(sim6_static_SN1, k = 2), GAIC(sim6_static_SN2, k = 2),GAIC(sim6_static_TF, k = 2), GAIC(sim6_static_TF2, k = 2),
          GAIC(sim6_static_GT, k = 2),GAIC(sim6_static_JSU, k = 2), GAIC(sim6_static_JSUo, k = 2), GAIC(sim6_static_NET, k = 2),
          GAIC(sim6_static_SHASH, k = 2),GAIC(sim6_static_SHASHo, k = 2), GAIC(sim6_static_SHASHo2, k = 2), GAIC(sim6_static_SEP2, k = 2), 
          GAIC(sim6_static_SEP3, k = 2), GAIC(sim6_static_SEP4, k = 2), GAIC(sim6_static_ST1, k = 2),
          GAIC(sim6_static_ST2, k = 2),GAIC(sim6_static_ST3, k = 2), GAIC(sim6_static_ST4, k = 2), GAIC(sim6_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim6_pAIC_1a <- ggplot(sim6_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
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

sim6_pAIC_1a + sim7_pAIC_1a + sim8_pAIC_1a + sim9_pAIC_1a
# param recovery ------------------------------------------------------------------
param_table6 <- data.frame(
  models = c("true","NO", "SEP3"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10",round(sim6_static_NO$mu.coefficients,3), round(sim6_static_SEP3$mu.coefficients,3)),
  sigma = c("2",round(exp(sim6_static_NO$sigma.coefficients),3), round(exp(sim6_static_SEP3$sigma.coefficients),3)),
  nu = c("1","-", round(exp(sim6_static_SEP3$nu.coefficients),3)),
  tau = c("1000","-", round(exp(sim6_static_SEP3$tau.coefficients),3))
)


# copmpare the two models PDF -----------------------------------------------------
model_compare6 <- list(NO = sim6_static_NO, SEP3 = sim6_static_SEP3)

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
sim6_sim <- bind_rows(
  lapply(names(model_compare6), function(name) {
    sim_from_model(model_compare6[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim6_p1 <- ggplot(dat6, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  labs(title = "Sim 1: symmetrical, normal tails", x = "Value", y = "Density") +
  xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim6_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table6 <- ggtexttable(param_table6, rows = NULL, theme = ttheme("light", base_size = 12))
table6 <- table_cell_bg(table6, row = 2, column = 1:5, linewidth = 5,
                     fill="grey", color = "grey")
# combeanation
sim6_plots <- (sim6_pAIC_1a  +  sim6_p1 / table6)



## AWESOME!  as expected, a model using a normal distirbution does best when the data is simulated to be 
## symmetrical and normal tails. Since, a flexible family disitrbution did not win here, tells us that
## our data cannot be assumed with normal distirbution, and is infact skewed and heavy tailed (for now)

################################################################################

# Best model: Normal distribution

sim6_static_NO <- gamlss(y ~ 1, family = NO(), data = dat6) 

# get model's predicted parameters
dat6$mu_hat <- predict(sim6_static_NO, what = "mu", type = "response")
dat6$sigma_hat <- predict(sim6_static_NO, what = "sigma", type = "response")
dat6$nu_hat <- predict(sim6_static_NO, what = "nu", type = "response")
dat6$tau_hat <- predict(sim6_static_NO, what = "tau", type = "response")

### get sim's sample moments
kurtosis(y6) - 3 # -0.02255661
skewness(y6) # 0.03352438
var(y6) # 3.574835
sd(y6) # 1.890723
mean(y6) # 10.01656

# use model's predicted parameters to simulate 
sim6_mu   <-     10.01656
sim6_sigma <-    1.887569

set.seed(5)
sim6_sim <- rNO(300, mu = sim6_mu, sigma = sim6_sigma)

# get skewness and kurtosis 
kurtosis(sim6_sim) - 3 # 0.1123197
skewness(sim6_sim) # -0.03767664
var(sim6_sim) # 3.43292
sd(sim6_sim) # 1.852814
mean(sim6_sim) # 10.04596

##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 

set.seed(1)
bootstrap_moments_sim_only <- function(y, B = 500) {
  n <- length(y)
  sk <- numeric(B); exk <- numeric(B); mu <- numeric(B); s2 <- numeric(B)
  for(b in 1:B) {
    samp <- sample(y, size = n, replace = TRUE)
    mu[b] <- mean(samp)
    s2[b] <- sd(samp)              # sample variance (n-1)
    sk[b] <- mean((samp - mean(samp))^3)/sd(samp)^3
    exk[b] <- mean((samp - mean(samp))^4)/sd(samp)^4 - 3
  }
  list(
    se = c(mean = sd(mu), var = sd(s2), skew = sd(sk), exkurt = sd(exk)),
    ci = rbind(
      mean = quantile(mu, c(0.025,0.975)),
      var  = quantile(s2, c(0.025,0.975)),
      skew = quantile(sk, c(0.025,0.975)),
      exk  = quantile(exk, c(0.025,0.975))
    ),
    boot_reps = list(mean = mu, var = s2, skew = sk, exkurt = exk)
  )
}

# get them errors gurll 
y6_se <- bootstrap_moments_sim_only(y6) 
sim6_sim_se <- bootstrap_moments_sim_only(sim6_sim) 

y6_se$se
sim6_sim_se$se

sim6_df <- tibble::tibble(
  sims = c("prefit", "postfit"),
  mu = c(10.01656,10.04596),
  mu_se = c(0.10856012, 0.1110101),
  var = c(1.890723,1.852814),
  var_se = c(0.07836133, 0.0758970),
  skew = c(0.03352438,-0.03767664),
  skew_se = c(0.14184126, 0.1486353),
  kurt = c(-0.02255661,0.1123197),
  kurt_se = c(0.28927074, 0.3200140)
)

################################################################################


##### SIM 7: simple, +ve skew, meso --> SST() ==================================================================


t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
set.seed(123)
y7 <- rSST(t2, mu = mu0, sigma = sigma0, nu = 3, tau = 1000)

# create dataset
dat7 <- data.frame(
  t = t_index2,
  y = y7,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 3,
  tau_true = 1000,
  dist = "SST"
)

# plot them to see
ggplot(dat7, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat7, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat7$y), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ### -------------------------------

# NO
sim7_static_NO <- gamlss(y ~ 1, family = NO(), data = dat7) # two-parameter
# NO2
sim7_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat7)
# GU
sim7_static_GU <- gamlss(y ~ 1, family = GU(), data = dat7)
# LO
sim7_static_LO <- gamlss(y ~ 1, family = LO(), data = dat7)
# RG
sim7_static_RG <- gamlss(y ~ 1, family = RG(), data = dat7)
# exGAUS
sim7_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat7, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim7_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat7)
# PE
sim7_static_PE <- gamlss(y ~ 1, family = PE(), data = dat7,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim7_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat7, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim7_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat7)
# SN2
sim7_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat7,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim7_static_TF <- gamlss(y ~ 1, family = TF(), data = dat7)
# TF2
sim7_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat7,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim7_static_GT <- gamlss(y ~ 1, family = GT(), data = dat7, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim7_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat7, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim7_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat7,
                           mu.start = mean(dat7$y), sigma.start = sd(dat7$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim7_static_NET <- gamlss(y ~ 1, family = NET(), data = dat7)
# SHASH
sim7_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat7,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim7_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat7,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim7_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat7,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim7_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat7, # not working
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim7_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat7,
                           mu.start = mean(dat7$y), sigma.start = sd(dat7$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim7_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat7,
                           mu.start = mean(dat7$y), sigma.start = sd(dat7$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim7_static_SST <- gamlss(y ~ 1, family = SST(), data = dat7, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim7_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat7, # not working
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim7_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat7, # not working
                          mu.start = mean(dat7$y), sigma.start = sd(dat7$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim7_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat7,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim7_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat7,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim7_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat7,
                          mu.start = mean(dat7$y), sigma.start = sd(dat7$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# compare models AIC ------------------------------------------------------------
sim7_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP3","SEP4", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim7_static_NO, k = 2),GAIC(sim7_static_NO2, k = 2), GAIC(sim7_static_GU, k = 2), GAIC(sim7_static_LO, k = 2), 
          GAIC(sim7_static_RG, k = 2), GAIC(sim7_static_exGAUS, k = 2), GAIC(sim7_static_NOF, k = 2), GAIC(sim7_static_PE, k = 2), 
          GAIC(sim7_static_PE2, k = 2),GAIC(sim7_static_SN1, k = 2), GAIC(sim7_static_SN2, k = 2),GAIC(sim7_static_TF, k = 2), GAIC(sim7_static_TF2, k = 2),
          GAIC(sim7_static_GT, k = 2),GAIC(sim7_static_JSU, k = 2), GAIC(sim7_static_JSUo, k = 2), GAIC(sim7_static_NET, k = 2),
          GAIC(sim7_static_SHASH, k = 2),GAIC(sim7_static_SHASHo, k = 2), GAIC(sim7_static_SHASHo2, k = 2),
          GAIC(sim7_static_SEP3, k = 2), GAIC(sim7_static_SEP4, k = 2), 
          GAIC(sim7_static_ST3, k = 2), GAIC(sim7_static_ST4, k = 2), GAIC(sim7_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis 
sim7_pAIC_1a <- ggplot(sim7_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
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


# parameter recovery ---------------------------------------------------------------
param_table7 <- data.frame(
  models = c("true","SN2", "SEP3"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10",round(sim7_static_SN2$mu.coefficients,3), round(sim7_static_SEP3$mu.coefficients,3)),
  sigma = c("2",round(exp(sim7_static_SN2$sigma.coefficients),3), round(exp(sim7_static_SEP3$sigma.coefficients),3)),
  nu = c("3",round(exp(sim7_static_SN2$nu.coefficients),3), round(exp(sim7_static_SEP3$nu.coefficients),3)),
  tau = c("1000","-", round(exp(sim7_static_SEP3$tau.coefficients),3))
)



# compare the two models PDF -----------------------------------------------------
model_compare7 <- list(SN2 = sim7_static_SN2, SEP3 = sim7_static_SEP3)


# Simulate and bind results
set.seed(123)
sim7_sim <- bind_rows(
  lapply(names(model_compare7), function(name) {
    sim_from_model(model_compare7[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate and model PDF --------------------
sim7_p1 <- ggplot(dat7, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  labs(title = "Sim 2: +ve skewness, normal tails",x = "Value", y = "Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim7_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table7 <- ggtexttable(param_table7, rows = NULL, theme = ttheme("light", base_size = 12))
table7 <- table_cell_bg(table7, row = 2, column = 1:5, linewidth = 5,
                        fill="grey", color = "grey")
# combeanation
sim7_plots <- (sim7_pAIC_1a  +  sim7_p1 / table7)



################################################################################

# Best model: SN2

sim7_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat7,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))


# get model's predicted parameters
dat7$mu_hat <- predict(sim7_static_SN2, what = "mu", type = "response")
dat7$sigma_hat <- predict(sim7_static_SN2, what = "sigma", type = "response")
dat7$nu_hat <- predict(sim7_static_SN2, what = "nu", type = "response")
dat7$tau_hat <- predict(sim7_static_SN2, what = "tau", type = "response")

### get sim's sample moments
kurtosis(y7) - 3 # 0.1937399
skewness(y7) # 0.8661981
var(y7) # 3.687103
sd(y7) # 1.920183
mean(y7) # 9.974478

# use model's predicted parameters to simulate 
sim7_mu   <-     7.937007
sim7_sigma <-    1.091041
sim7_nu <-    2.727621

set.seed(5)
sim7_sim <- rSN2(300, mu = sim7_mu, sigma = sim7_sigma, nu = sim7_nu)

# get skewness and kurtosis 
kurtosis(sim7_sim) - 3 # 0.1503916
skewness(sim7_sim) # 0.8429814
var(sim7_sim) # 3.750385
sd(sim7_sim) # 1.936591
mean(sim7_sim) # 9.978622

##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 


# get them errors gurll 
set.seed(1)
y7_se <- bootstrap_moments_sim_only(y7) 
sim7_sim_se <- bootstrap_moments_sim_only(sim7_sim) 

y7_se$se
sim7_sim_se$se

sim7_df <- tibble::tibble(
  sims = c("prefit", "postfit"),
  mu = c(9.974478,9.978622),
  mu_se = c(0.11127801, 0.11173760),
  var = c(1.920183,1.936591),
  var_se = c(0.08221278, 0.08226808),
  skew = c(0.8661981,0.8429814),
  skew_se = c(0.10807973, 0.11150256),
  kurt = c(0.1937399,0.1503916),
  kurt_se = c(0.31469208, 0.30930100)
)

################################################################################



##### SIM 8: simple, sym, lepto --> SST() ==================================================================


t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
set.seed(123)
y8 <- rSST(t2, mu = mu0, sigma = sigma0, nu = 1, tau = 4)

# create dataset
dat8 <- data.frame(
  t = t_index2,
  y = y8,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 1,
  tau_true = 4,
  dist = "SST"
)

# plot them to see
ggplot(dat8, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 

# look at density plot
ggplot(dat8, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat8$y), linetype = "dashed", col = "grey") +
  theme_minimal()



### fit models of different family distribution ### --------------------------------

# NO
sim8_static_NO <- gamlss(y ~ 1, family = NO(), data = dat8) # two-parameter
# NO2
sim8_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat8)
# GU
sim8_static_GU <- gamlss(y ~ 1, family = GU(), data = dat8)
# LO
sim8_static_LO <- gamlss(y ~ 1, family = LO(), data = dat8)
# RG
sim8_static_RG <- gamlss(y ~ 1, family = RG(), data = dat8)
# exGAUS
sim8_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat8, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim8_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat8)
# PE
sim8_static_PE <- gamlss(y ~ 1, family = PE(), data = dat8,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim8_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat8, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim8_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat8)
# SN2
sim8_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat8,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim8_static_TF <- gamlss(y ~ 1, family = TF(), data = dat8)
# TF2
sim8_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat8,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim8_static_GT <- gamlss(y ~ 1, family = GT(), data = dat8, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim8_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat8, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim8_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat8,
                           mu.start = mean(dat8$y), sigma.start = sd(dat8$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim8_static_NET <- gamlss(y ~ 1, family = NET(), data = dat8)
# SHASH
sim8_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat8,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim8_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat8,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim8_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat8,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim8_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat8, # not working
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim8_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat8,
                           mu.start = mean(dat8$y), sigma.start = sd(dat8$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim8_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat8,
                           mu.start = mean(dat8$y), sigma.start = sd(dat8$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim8_static_SST <- gamlss(y ~ 1, family = SST(), data = dat8, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim8_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat8, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim8_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat8, 
                          mu.start = mean(dat8$y), sigma.start = sd(dat8$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim8_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat8,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim8_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat8,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim8_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat8,
                          mu.start = mean(dat8$y), sigma.start = sd(dat8$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# compare model AIC -------------------------------------------------------------
sim8_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "SST","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim8_static_NO, k = 2),GAIC(sim8_static_NO2, k = 2), GAIC(sim8_static_GU, k = 2), GAIC(sim8_static_LO, k = 2), 
          GAIC(sim8_static_RG, k = 2), GAIC(sim8_static_exGAUS, k = 2), GAIC(sim8_static_NOF, k = 2), GAIC(sim8_static_PE, k = 2), 
          GAIC(sim8_static_PE2, k = 2),GAIC(sim8_static_SN1, k = 2), GAIC(sim8_static_SN2, k = 2),GAIC(sim8_static_TF, k = 2), GAIC(sim8_static_TF2, k = 2),
          GAIC(sim8_static_GT, k = 2),GAIC(sim8_static_JSU, k = 2), GAIC(sim8_static_JSUo, k = 2), GAIC(sim8_static_NET, k = 2),
          GAIC(sim8_static_SHASH, k = 2),GAIC(sim8_static_SHASHo, k = 2), GAIC(sim8_static_SHASHo2, k = 2), GAIC(sim8_static_SEP2, k = 2), 
          GAIC(sim8_static_SEP3, k = 2), GAIC(sim8_static_SEP4, k = 2), GAIC(sim8_static_SST, k = 2),GAIC(sim8_static_ST1, k = 2),
          GAIC(sim8_static_ST2, k = 2),GAIC(sim8_static_ST3, k = 2), GAIC(sim8_static_ST4, k = 2), GAIC(sim8_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto","lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis --------------------------------
sim8_pAIC_1a <- ggplot(sim8_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
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


# parameter recovery ---------------------------------------------------------------
param_table8 <- data.frame(
  models = c("true","SST", "TF", "SEP2"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10",round(sim8_static_SST$mu.coefficients,3), round(sim8_static_TF$mu.coefficients,3), round(sim8_static_SEP2$mu.coefficients,3)),
  sigma = c("2",round(exp(sim8_static_SST$sigma.coefficients),3), round(exp(sim8_static_TF$sigma.coefficients),3), round(exp(sim8_static_SEP2$sigma.coefficients),3)),
  nu = c("1",round(exp(sim8_static_SST$nu.coefficients),3), "-", round(sim8_static_SEP2$nu.coefficients,3)),
  tau = c("4", round((exp(sim8_static_SST$tau.coefficients)+2),3), round(exp(sim8_static_TF$nu.coefficients),3),round(exp(sim8_static_SEP2$tau.coefficients),3))
)



# compare the two models PDF -----------------------------------------------------
model_compare8 <- list(SST = sim8_static_SST, TF = sim8_static_TF ,SEP2 = sim8_static_SEP2)


# Simulate and bind results
set.seed(123)
sim8_sim <- bind_rows(
  lapply(names(model_compare8), function(name) {
    sim_from_model(model_compare8[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim8_p1 <- ggplot(dat8, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  labs(title = "Sim 3: symmetrical, heavy tails", x = "Value", y = "Density") +
  xlim(x = c(0,20)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim8_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 


# combine density plot, AIC ranks and param recovery table -----------------------
table8 <- ggtexttable(param_table8, rows = NULL, theme = ttheme("light", base_size = 12))
table8 <- table_cell_bg(table8, row = 2, column = 1:5, linewidth = 5,
                        fill="grey", color = "grey")
# combeanation
sim8_plots <- (sim8_pAIC_1a  +  sim8_p1 / table8)


################################################################################

# Best model: SST

sim8_static_SST <- gamlss(y ~ 1, family = SST(), data = dat8, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  

resid_wp(sim8_static_SST)
moment_bucket(sim8_static_SST)

# get model's predicted parameters
dat8$mu_hat <- predict(sim8_static_SST, what = "mu", type = "response")
dat8$sigma_hat <- predict(sim8_static_SST, what = "sigma", type = "response")
dat8$nu_hat <- predict(sim8_static_SST, what = "nu", type = "response")
dat8$tau_hat <- predict(sim8_static_SST, what = "tau", type = "response")

### get sim's sample moments
kurtosis(y8) - 3 # 5.551602
skewness(y8) # -0.6345648
var(y8) # 3.280258
sd(y8)  # 1.811148
mean(y8) # 10.01194

# use model's predicted parameters to simulate 
sim8_mu   <-     10.04387
sim8_sigma <-    1.807381
sim8_nu <-    1.099435
sim8_tau <-    4.561637

set.seed(5)
sim8_sim <- rSST(300, mu = sim8_mu, sigma = sim8_sigma, nu = sim8_nu, tau = sim8_tau)

# get skewness and kurtosis 
kurtosis(sim8_sim) - 3 # 0.9121725
skewness(sim8_sim) # 0.454769
var(sim8_sim) # 3.030753
sd(sim8_sim) # 1.740906
mean(sim8_sim) # 10.03635

##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 


# get them errors gurll 
set.seed(1)
y8_se <- bootstrap_moments_sim_only(y8) 
sim8_sim_se <- bootstrap_moments_sim_only(sim8_sim) 

y8_se$se
sim8_sim_se$se

sim8_df <- tibble::tibble(
  sims = c("prefit", "postfit"),
  mu = c(10.01194,10.03635),
  mu_se = c(0.1036661, 0.09951247),
  var = c(1.811148,1.740906),
  var_se = c(0.1457338, 0.08605180),
  skew = c(-0.6345648,0.454769),
  skew_se = c(0.6645056, 0.17948539),
  kurt = c(5.551602,0.9121725),
  kurt_se = c(3.1506029, 0.39187347)
)

#################################################################################


##### SIM 9: simple, +ve, lepto --> SST() ==================================================================


t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
set.seed(123)
y9 <- rSST(t2, mu = mu0, sigma = sigma0, nu = 3, tau = 4)

# create dataset
dat9 <- data.frame(
  t = t_index2,
  y = y9,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 3,
  tau_true = 4,
  dist = "SST"
)

# plot them to see
ggplot(dat9, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat9, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat9$y), linetype = "dashed", col = "grey") +
  theme_minimal()



### fit models of different family distribution ### ------------------------------

# NO
sim9_static_NO <- gamlss(y ~ 1, family = NO(), data = dat9) # two-parameter
# NO2
sim9_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat9)
# GU
sim9_static_GU <- gamlss(y ~ 1, family = GU(), data = dat9)
# LO
sim9_static_LO <- gamlss(y ~ 1, family = LO(), data = dat9)
# RG
sim9_static_RG <- gamlss(y ~ 1, family = RG(), data = dat9)
# exGAUS
sim9_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat9, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim9_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat9)
# PE
sim9_static_PE <- gamlss(y ~ 1, family = PE(), data = dat9,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim9_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat9, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim9_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat9)
# SN2
sim9_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat9,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim9_static_TF <- gamlss(y ~ 1, family = TF(), data = dat9)
# TF2
sim9_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat9,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim9_static_GT <- gamlss(y ~ 1, family = GT(), data = dat9, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim9_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat9, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim9_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat9,
                           mu.start = mean(dat9$y), sigma.start = sd(dat9$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim9_static_NET <- gamlss(y ~ 1, family = NET(), data = dat9)
# SHASH
sim9_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat9,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim9_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat9,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim9_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat9,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim9_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat9, # not working
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim9_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat9,
                           mu.start = mean(dat9$y), sigma.start = sd(dat9$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim9_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat9,
                           mu.start = mean(dat9$y), sigma.start = sd(dat9$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim9_static_SST <- gamlss(y ~ 1, family = SST(), data = dat9, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim9_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat9, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim9_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat9, 
                          mu.start = mean(dat9$y), sigma.start = sd(dat9$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim9_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat9,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim9_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat9,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim9_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat9,
                          mu.start = mean(dat9$y), sigma.start = sd(dat9$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# compare model AIC ---------------------------------------------------------------
sim9_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "SST","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim9_static_NO, k = 2),GAIC(sim9_static_NO2, k = 2), GAIC(sim9_static_GU, k = 2), GAIC(sim9_static_LO, k = 2), 
          GAIC(sim9_static_RG, k = 2), GAIC(sim9_static_exGAUS, k = 2), GAIC(sim9_static_NOF, k = 2), GAIC(sim9_static_PE, k = 2), 
          GAIC(sim9_static_PE2, k = 2),GAIC(sim9_static_SN1, k = 2), GAIC(sim9_static_SN2, k = 2),GAIC(sim9_static_TF, k = 2), GAIC(sim9_static_TF2, k = 2),
          GAIC(sim9_static_GT, k = 2),GAIC(sim9_static_JSU, k = 2), GAIC(sim9_static_JSUo, k = 2), GAIC(sim9_static_NET, k = 2),
          GAIC(sim9_static_SHASH, k = 2),GAIC(sim9_static_SHASHo, k = 2), GAIC(sim9_static_SHASHo2, k = 2), GAIC(sim9_static_SEP2, k = 2), 
          GAIC(sim9_static_SEP3, k = 2), GAIC(sim9_static_SEP4, k = 2), GAIC(sim9_static_SST, k = 2),GAIC(sim9_static_ST1, k = 2),
          GAIC(sim9_static_ST2, k = 2),GAIC(sim9_static_ST3, k = 2), GAIC(sim9_static_ST4, k = 2), GAIC(sim9_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto","lepto", "lepto", "lepto", "lepto")
)


# plot ----------------------------------------------------------------------------
sim9_pAIC_1a <- ggplot(sim9_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
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



# parameter recovery ---------------------------------------------------------------
param_table9 <- data.frame(
  models = c("true","ST3", "SEP2"),
  #param = c("-", "2", "4"),
  #AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SEP3),2)),
  #deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10",round(sim9_static_ST3$mu.coefficients,3), round(sim9_static_SEP2$mu.coefficients,3)),
  sigma = c("2",round(exp(sim9_static_ST3$sigma.coefficients),3), round(exp(sim9_static_SEP2$sigma.coefficients),3)),
  nu = c("1",round(exp(sim8_static_ST3$nu.coefficients),3), round(sim9_static_SEP2$nu.coefficients,3)),
  tau = c("4", round(exp(sim8_static_ST3$tau.coefficients),3), round(exp(sim9_static_SEP2$tau.coefficients),3))
)



# compare the two models PDF -----------------------------------------------------
model_compare9 <- list(ST3 = sim9_static_ST3, SEP2 = sim9_static_SEP2)


# Simulate and bind results
set.seed(123)
sim9_sim <- bind_rows(
  lapply(names(model_compare9), function(name) {
    sim_from_model(model_compare9[[name]], n = 2000) %>%
      mutate(model = name)
  }),
  .id = "model_id"
)


# histogram and non-parametric density estimate ----------------------------------
sim9_p1 <- ggplot(dat9, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  labs(title = "Sim 4: +ve skewness, heavy tails", x = "Value", y = "Density") +
  xlim(x = c(6,20)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  geom_density(data = sim9_sim, aes(x = value, color = model), linewidth = 1) +
  theme_bw() 



# combine density plot, AIC ranks and param recovery table -----------------------
table9 <- ggtexttable(param_table9, rows = NULL, theme = ttheme("light", base_size = 12))
table9 <- table_cell_bg(table9, row = 2, column = 1:5, linewidth = 5,
                        fill="grey", color = "grey")
# combeanation
sim9_plots <- (sim9_pAIC_1a  +  sim9_p1 / table9)


################################################################################

# Best model: ST3

sim9_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat9,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

moment_bucket(sim9_static_ST3)

# get model's predicted parameters
dat9$mu_hat <- predict(sim9_static_ST3, what = "mu", type = "response")
dat9$sigma_hat <- predict(sim9_static_ST3, what = "sigma", type = "response")
dat9$nu_hat <- predict(sim9_static_ST3, what = "nu", type = "response")
dat9$tau_hat <- predict(sim9_static_ST3, what = "tau", type = "response")

### get sim's sample moments
kurtosis(y9) - 3 # 3.21435
skewness(y9) # 1.612569
var(y9) # 3.036752
sd(y9) # 1.742628
mean(y9) # 9.946443

# use model's predicted parameters to simulate 
sim9_mu   <-     8.287107
sim9_sigma <-    0.6162365
sim9_nu <-    3.054366
sim9_tau <-    3.94429

set.seed(5)
sim9_sim <- rST3(300, mu = sim9_mu, sigma = sim9_sigma, nu = sim9_nu, tau = sim9_tau)

# get skewness and kurtosis 
kurtosis(sim9_sim) - 3 # 4.270903
skewness(sim9_sim) # 1.80133
var(sim9_sim) # 3.302383
sd(sim9_sim) # 1.817246
mean(sim9_sim) # 9.961737

##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 

# get them errors gurll 
set.seed(1)
y9_se <- bootstrap_moments_sim_only(y9) 
sim9_sim_se <- bootstrap_moments_sim_only(sim9_sim) 

y9_se$se
sim9_sim_se$se

sim9_df <- tibble::tibble(
  sims = c("prefit", "postfit"),
  mu = c(9.946443,9.961737),
  mu_se = c(0.1008930, 0.1052851),
  var = c(1.742628,1.817246),
  var_se = c(0.1144540, 0.1284349),
  skew = c(1.612569,1.80133),
  skew_se = c(0.2129338, 0.2193494),
  kurt = c(3.21435,4.270903),
  kurt_se = c(1.1787498, 1.2083367)
)

#################################################################################


##### SIM 15:  sym,  light --> SHASHo() ==================================================================

t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# changing through time
beta_tau   <- -0.005    
tau_t   <- 2.4   + beta_tau * t_index2

# observation randomly taken from SST
set.seed(123)
y15 <- rSHASHo(t2, mu = mu0, sigma = sigma0, nu = 0, tau = 2)

# create dataset
dat15 <- data.frame(
  t = t_index2,
  y = y15,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 0,
  tau_true = 2,
  dist = "SHASHo"
)

ggplot(dat15, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 

# look at density plot
ggplot(dat14, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat14), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim15_static_NO <- gamlss(y ~ 1, family = NO(), data = dat15) # two-parameter
# NO2
sim15_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat15)
# GU
sim15_static_GU <- gamlss(y ~ 1, family = GU(), data = dat15)
# LO
sim15_static_LO <- gamlss(y ~ 1, family = LO(), data = dat15)
# RG
sim15_static_RG <- gamlss(y ~ 1, family = RG(), data = dat15)
# exGAUS
sim15_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat15, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim15_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat15,
                           method = mixed(),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim15_static_PE <- gamlss(y ~ 1, family = PE(), data = dat15,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim15_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat15, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim15_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat15) 
# SN2
sim15_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat15, 
                           mu.start = mean(dat15$y), sigma.start = sd(dat15$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim15_static_TF <- gamlss(y ~ 1, family = TF(), data = dat15)
# TF2
sim15_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat15,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim15_static_GT <- gamlss(y ~ 1, family = GT(), data = dat15, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim15_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat15, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim15_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat15,
                            mu.start = mean(dat15$y), sigma.start = sd(dat15$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim15_static_NET <- gamlss(y ~ 1, family = NET(), data = dat15)
# SHASH
sim15_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat15,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim15_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat15,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim15_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat15,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim15_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat15,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim15_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat15,
                            mu.start = mean(dat15$y), sigma.start = sd(dat15$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim15_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat15, 
                            mu.start = mean(dat15$y), sigma.start = sd(dat15$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim15_static_SST <- gamlss(y ~ 1, family = SST(), data = dat15,  # error
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim15_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat15,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim15_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat15,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim15_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat15,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim15_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat15,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim15_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat15,
                           mu.start = mean(dat15$y), sigma.start = sd(dat15$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim15_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "SEP4","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim15_static_NO, k = 2),GAIC(sim15_static_NO2, k = 2), GAIC(sim15_static_GU, k = 2), GAIC(sim15_static_LO, k = 2), 
          GAIC(sim15_static_RG, k = 2), GAIC(sim15_static_exGAUS, k = 2), GAIC(sim15_static_NOF, k = 2), GAIC(sim15_static_PE, k = 2), 
          GAIC(sim15_static_PE2, k = 2),GAIC(sim15_static_SN1, k = 2), GAIC(sim15_static_SN2, k = 2),GAIC(sim15_static_TF, k = 2), GAIC(sim15_static_TF2, k = 2),
          GAIC(sim15_static_GT, k = 2),GAIC(sim15_static_JSU, k = 2), GAIC(sim15_static_JSUo, k = 2), GAIC(sim15_static_NET, k = 2),
          GAIC(sim15_static_SHASH, k = 2),GAIC(sim15_static_SHASHo, k = 2), GAIC(sim15_static_SHASHo2, k = 2), GAIC(sim15_static_SEP2, k = 2), 
          GAIC(sim15_static_SEP3, k = 2),  GAIC(sim15_static_SEP4, k = 2), GAIC(sim15_static_ST1, k = 2),
          GAIC(sim15_static_ST2, k = 2),GAIC(sim15_static_ST3, k = 2), GAIC(sim15_static_ST4, k = 2), GAIC(sim15_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both", "both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim15_pAIC_1a <- ggplot(sim15_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "+ve skew, light tails",x = "AIC", y = "models") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



################################################################################

# Best model: PE

sim15_static_PE <- gamlss(y ~ 1, family = PE(), data = dat15,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# get model's predicted parameters
dat15$mu_hat <- predict(sim15_static_PE, what = "mu", type = "response")
dat15$sigma_hat <- predict(sim15_static_PE, what = "sigma", type = "response")
dat15$tau_hat <- predict(sim15_static_PE, what = "nu", type = "response")

### get sim's sample moments
kurtosis(y15) - 3 # -0.6839482
skewness(y15) # 0.06492624
var(y15) # 0.6499832
sd(y15) # 0.8062154
mean(y15) # 10.00357

# use model's predicted parameters to simulate 
sim15_mu   <-     10.02
sim15_sigma <-    0.8050112
sim15_tau <-    3.374241

set.seed(5)
sim15_sim <- rPE(300, mu = sim15_mu, sigma = sim15_sigma, nu = sim15_tau)

# get skewness and kurtosis 
kurtosis(sim15_sim) - 3 # -0.8885583
skewness(sim15_sim) # 0.07785684
var(sim15_sim) # 0.672405
sd(sim15_sim) # 0.8200031
mean(sim15_sim) # 10.00699


##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 

# get them errors gurll 
set.seed(1)
y15_se <- bootstrap_moments_sim_only(y15) 
sim15_sim_se <- bootstrap_moments_sim_only(sim15_sim) 

y15_se$se
sim15_sim_se$se

sim15_df <- tibble::tibble(
  sims = c("prefit", "postfit"),
  mu = c(10.00357,10.00699),
  mu_se = c(0.04638783, 0.04668439),
  var = c(0.8062154,0.8200031),
  var_se = c(0.02713392, 0.02542587),
  skew = c(0.06492624,0.07785684),
  skew_se = c(0.09458095, 0.08656070),
  kurt = c(-0.6839482,-0.8885583),
  kurt_se = c(0.13143055, 0.09570030)
)

#################################################################################


##### SIM 14:, +ve skewness,  light  --> SHASHo() ==================================================================


t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# changing through time
beta_tau   <- -0.005    
tau_t   <- 2.4   + beta_tau * t_index2

# observation randomly taken from SST
set.seed(123)
y14 <- rSHASHo(t2, mu = mu0, sigma = sigma0, nu = 2, tau = 2)

# create dataset
dat14 <- data.frame(
  t = t_index2,
  y = y14,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = 2,
  tau_true = 2,
  dist = "SHASHo"
)

ggplot(dat14, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 

# look at density plot
ggplot(dat14, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat14), linetype = "dashed", col = "grey") +
  theme_minimal()


### fit models of different family distribution ###-----------------------------

# NO
sim14_static_NO <- gamlss(y ~ 1, family = NO(), data = dat14) # two-parameter
# NO2
sim14_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat14)
# GU
sim14_static_GU <- gamlss(y ~ 1, family = GU(), data = dat14)
# LO
sim14_static_LO <- gamlss(y ~ 1, family = LO(), data = dat14)
# RG
sim14_static_RG <- gamlss(y ~ 1, family = RG(), data = dat14)
# exGAUS
sim14_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat14, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim14_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat14,
                           method = mixed(),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.0001, trace = FALSE))
# PE
sim14_static_PE <- gamlss(y ~ 1, family = PE(), data = dat14,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim14_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat14, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim14_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat14) 
# SN2
sim14_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat14, 
                           mu.start = mean(dat14$y), sigma.start = sd(dat14$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim14_static_TF <- gamlss(y ~ 1, family = TF(), data = dat14)
# TF2
sim14_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat14,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim14_static_GT <- gamlss(y ~ 1, family = GT(), data = dat14, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim14_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat14, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim14_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat14,
                            mu.start = mean(dat14$y), sigma.start = sd(dat14$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim14_static_NET <- gamlss(y ~ 1, family = NET(), data = dat14)
# SHASH
sim14_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat14,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim14_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat14,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim14_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat14,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim14_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat14,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim14_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat14,
                            mu.start = mean(dat14$y), sigma.start = sd(dat14$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim14_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat14, 
                            mu.start = mean(dat14$y), sigma.start = sd(dat14$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim14_static_SST <- gamlss(y ~ 1, family = SST(), data = dat14,  # error
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim14_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat14,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim14_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat14,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim14_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat14,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim14_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat14,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim14_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat14,
                           mu.start = mean(dat14$y), sigma.start = sd(dat14$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# get the AIC  -------------------------------------------------------------------
sim14_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "SN2","TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3", "SEP4","ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim14_static_NO, k = 2),GAIC(sim14_static_NO2, k = 2), GAIC(sim14_static_GU, k = 2), GAIC(sim14_static_LO, k = 2), 
          GAIC(sim14_static_RG, k = 2), GAIC(sim14_static_exGAUS, k = 2), GAIC(sim14_static_NOF, k = 2), GAIC(sim14_static_PE, k = 2), 
          GAIC(sim14_static_PE2, k = 2),GAIC(sim14_static_SN1, k = 2), GAIC(sim14_static_SN2, k = 2),GAIC(sim14_static_TF, k = 2), GAIC(sim14_static_TF2, k = 2),
          GAIC(sim14_static_GT, k = 2),GAIC(sim14_static_JSU, k = 2), GAIC(sim14_static_JSUo, k = 2), GAIC(sim14_static_NET, k = 2),
          GAIC(sim14_static_SHASH, k = 2),GAIC(sim14_static_SHASHo, k = 2), GAIC(sim14_static_SHASHo2, k = 2), GAIC(sim14_static_SEP2, k = 2), 
          GAIC(sim14_static_SEP3, k = 2),  GAIC(sim14_static_SEP4, k = 2), GAIC(sim14_static_ST1, k = 2),
          GAIC(sim14_static_ST2, k = 2),GAIC(sim14_static_ST3, k = 2), GAIC(sim14_static_ST4, k = 2), GAIC(sim14_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3k","3k","3s","3s","3k","3k",
            "4k","4","4","2","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both","both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both", "both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "meso","meso", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto", "lepto", "lepto", "lepto")
)


# plot colored by skewness and shape by kurtosis----------------------------------
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")
# plot them
sim14_pAIC_1a <- ggplot(sim14_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "+ve skew, light tails",x = "AIC", y = "models") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


################################################################################

# Best model: SHASH

sim14_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat14,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))

# get model's predicted parameters
dat14$mu_hat <- predict(sim14_static_SHASH, what = "mu", type = "response")
dat14$sigma_hat <- predict(sim14_static_SHASH, what = "sigma", type = "response")
dat14$nu_hat <- predict(sim14_static_SHASH, what = "nu", type = "response")
dat14$tau_hat <- predict(sim14_static_SHASH, what = "tau", type = "response")


### get sim's sample moments
kurtosis(y14) - 3 # -0.4814689
skewness(y14) # 0.5641898
var(y14) # 1.616886
sd(y14) # 1.271568
mean(y14) # 12.53183

# use model's predicted parameters to simulate 
sim14_mu   <-     12.30494
sim14_sigma <-    4.68395
sim14_nu   <-     4.285368
sim14_tau <-    2.333609

set.seed(5)
sim14_sim <- rSHASH(300, mu = sim14_mu, sigma = sim14_sigma, nu = sim14_nu, tau = sim14_tau)

# get skewness and kurtosis 
kurtosis(sim14_sim) - 3 # -0.2324722
skewness(sim14_sim) # 0.625585
var(sim14_sim) # 1.679167
sd(sim14_sim) # 1.295827
mean(sim14_sim) # 12.54791


##======= getting standard errors =========##

##### simulation sampling variability by bootstrapping 

# get them errors gurll 
set.seed(1)
y14_se <- bootstrap_moments_sim_only(y14) 
sim14_sim_se <- bootstrap_moments_sim_only(sim14_sim) 

y14_se$se
sim14_sim_se$se

sim14_df <- tibble::tibble(
  sims = c("prefit", "postfit"),
  mu = c(12.53183,12.54791),
  mu_se = c(0.07361197, 0.07474666),
  var = c(1.271568,1.295827),
  var_se = c(0.04559036, 0.05182901),
  skew = c(0.5641898,0.625585),
  skew_se = c(0.09292410, 0.10829907),
  kurt = c(-0.4814689,-0.2324722),
  kurt_se = c(0.18500333, 0.24018023)
)

#################################################################################


constant_sim_df <- list(
  sim6  = sim6_df,
  sim7  = sim7_df,
  sim8  = sim8_df,
  sim9  = sim9_df,
  sim14 = sim14_df,
  sim15 = sim15_df
)


constant_sim_long <- bind_rows(constant_sim_df, .id = "sim_id") %>%
  pivot_longer(
    cols = -c(sim_id, sims),
    names_to = c("moment", "type"),
    names_sep = "_",
    values_to = "val"
  ) %>%
  mutate(type = replace_na(type, "value")) %>%
  pivot_wider(names_from = type, values_from = val) %>%
  mutate(
    sim_id = factor(sim_id, levels = c("sim6","sim7","sim8","sim9","sim15","sim14")),
    sims   = factor(sims, levels = c("prefit","postfit")),
    moment_lab = recode(moment,
                        mu   = "Mean",
                        var  = "Standard deviation",
                        skew = "Skewness",
                        kurt = "Kurtosis")
  )


truth_mean_sd <- tibble(
  moment_lab = c("Mean", "Standard deviation"),
  yint       = c(10, 2)
)
truth_sk_k <- tibble(
  moment_lab = c("Skewness", "Kurtosis"),
  yint       = c(0, 0)
)


ggplot(constant_sim_long,
       aes(x = sim_id, y = value, color = sims, group = sims)) +
  geom_point(position = position_dodge(width = 0.55), size = 2) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                position = position_dodge(width = 0.55), width = 0.25) +
  facet_wrap(~ moment_lab, scales = "free_y") +
  geom_hline(data = truth_mean_sd, aes(yintercept = yint), linetype = 2) +
  geom_hline(data = truth_sk_k,   aes(yintercept = yint), linetype = 2, color = "darkred") +
  scale_color_manual(values = c(prefit = "black", postfit = "azure4"),
                     breaks = c("prefit","postfit"),
                     labels = c("prefit","postfit")) +
  scale_x_discrete(labels = c(
    sim6  = "sim 1",
    sim7  = "sim 2",
    sim8  = "sim 3",
    sim9  = "sim 4",
    sim15 = "sim 5",
    sim14 = "sim 6"
  )) +
  labs(x = "Simulations", y = "Sample distributional moments", color = "") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

