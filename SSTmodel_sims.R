# =========================
# sim set up
# =========================

n        <- 800
t_index  <- seq_len(n)
t_scaled <- (t_index - 1)/(n - 1)
mu0       <- 20
sigma0   <- 1
nu0      <- 3  # For SST, 1 = symmetrical, > 1 = right skewed, < 1 = left skewed
tau0     <- 5  # for TF2, nu = 10 means close to normal, the lower the number the heavier the tails
# for SST, large tau (e.g 7) thinner tails, close to normal. Small tau means fatter tails 

# down-samples 
idx <- seq(1, n, by = 10)
idx2 <- seq(1, n, by = 2)
idx5 <- seq(1, n, by = 5)

# Slopes (linear change rates)
beta <- 2   # scale increases by 2 over full range of t_scaled
beta_nu    <- 0.3   # skewness increases linearly
beta_tau   <- -1    # kurtosis decreases linearly

# changes with time
mu_t <- pmax(0, mu0 + beta * t_scaled)
sigma_t <- pmax(0, sigma0 + beta * t_scaled)
nu_t    <- nu0    + beta_nu  * t_scaled
tau_t   <- tau0   + beta_tau * t_scaled


# Using SST distribution 
set.seed(24)
y <- rSST(n, mu = mu0, sigma = sigma_t, nu = nu0, tau = tau0)
y_down10 <- y[idx] # only get every 10th observation
y_down2 <- y[idx2] 
y_down5 <- y[idx5] 


# =========================
# create the sim
# =========================

set.seed(42)
n_sim <- 10
results <- vector("list", n_sim)

for (i in 1:n_sim){
  
  # Using SST distribution
  y <- rSST(n, mu = mu0, sigma = sigma0, nu = nu0, tau = tau0)
  
  # the sim 
  dat <- data.frame(
    t = t_index,
    t_scaled = t_scaled,
    y = y,
    sd_true = sigma0,
    mu_true = mu0,
    nu_true = nu0,
    tau_true = tau0,
    samp = "NO"
  )
  
  # fit a model
  SSTdat <- gamlss(y ~ 1, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1, 
                   family = SST(), data = dat, 
                   method = mixed(5, 200),
                   control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE)
  )
  
  results[[i]] <- coefAll(SSTdat)

}


View(results)



# normal temporal resolution of 800 observations 
dat <- data.frame(
  t = t_index,
  t_scaled = t_scaled,
  y = y,
  sd_true = sigma_t,
  mu_true = mu0,
  nu_true = nu0,
  tau_true = tau0,
  samp = "NO"
)

# only for every 2nd observation
dat_down2 <- data.frame(
  t = t_index[idx2],
  t_scaled = t_scaled[idx2],
  y = y_down2,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = nu0,
  tau_true = tau0,
  samp = "2"
)

# only for every 5th observation
dat_down5 <- data.frame(
  t = t_index[idx5],
  t_scaled = t_scaled[idx5],
  y = y_down5,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = nu0,
  tau_true = tau0,
  samp = "5"
)

# only for every 10th observation
dat_down10 <- data.frame(
  t = t_index[idx],
  t_scaled = t_scaled[idx],
  y = y_down10,
  sd_true = sigma_t,
  mu_true = mu0,
  nu_true = nu0,
  tau_true = tau0,
  samp = "10"
)


# =========================
# plot the sim
# =========================

#  time series plot for normal temporal resolution
ggplot(dat, aes(x = t)) +
  geom_line(
    aes(y = y),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)


#  time series plot for every 2nd observation
ggplot(dat_down2, aes(x = t)) +
  geom_line(
    aes(y = y),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)

#  time series plot for every 5th observation
ggplot(dat_down5, aes(x = t)) +
  geom_line(
    aes(y = y),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)


#  time series plot for every 10th observation
ggplot(dat_down10, aes(x = t)) +
  geom_line(
    aes(y = y),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)



dat_all <- bind_rows(dat, dat_down2, dat_down5, dat_down10)

# density plots
ggplot(dat_all, aes(y = y)) +
  geom_density(aes(color = samp),
               linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


# =========================
# SST models on different sims
# =========================

SSTdat <- gamlss(y ~ 1, sigma.fo = ~ t_scaled, nu.fo = ~ 1, tau.fo = ~ 1, family = SST(), data = dat, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE)
)

SSTdat_down10 <- gamlss(y ~ 1, sigma.fo = ~ t_scaled, nu.fo = ~ 1, tau.fo = ~ 1, family = SST(), data = dat_down10, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE)
)


SSTdat_down2 <- gamlss(y ~ 1, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1, family = SST(), data = dat_down2, 
                        method = mixed(5, 200),
                        control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE)
)

SSTdat_down5 <- gamlss(y ~ 1, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1, family = SST(), data = dat_down5, 
                       method = mixed(5, 200),
                       control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE)
)

summary(SSTdat)
summary(SSTdatc)
summary(SSTdat_down10)
summary(SSTdat_down2)
summary(SSTdat_down5)


# =========================
# parameter recovery for each model
# =========================

param_recov <- data.frame(
  param = c("mean", "SD","skewness","kurtosis"), 
  trueparam = (c(mu0, sigma0, nu0, tau0)),
  normal   = c(coefAll(SSTdat)$mu, exp(coefAll(SSTdat)$sigma), exp(coefAll(SSTdat)$nu),
              (exp(coefAll(SSTdat)$tau) + 2)), 
  down2   = c(coefAll(SSTdat_down2)$mu, exp(coefAll(SSTdat_down2)$sigma), exp(coefAll(SSTdat_down2)$nu),
              (exp(coefAll(SSTdat_down2)$tau) + 2)),
  down5   = c(coefAll(SSTdat_down5)$mu, exp(coefAll(SSTdat_down5)$sigma), exp(coefAll(SSTdat_down5)$nu),
              (exp(coefAll(SSTdat_down5)$tau) + 2)), 
  down10   = c(coefAll(SSTdat_down10)$mu, exp(coefAll(SSTdat_down10)$sigma), exp(coefAll(SSTdat_down10)$nu),
              (exp(coefAll(SSTdat_down10)$tau) + 2))
)
print(param_recov)

tail(sigma_t, 50)


    