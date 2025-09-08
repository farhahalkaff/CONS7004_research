# =========================
# sim set up
# =========================

n        <- 800
t_index  <- seq_len(n)
t_scaled <- (t_index - 1)/(n - 1)
mu0       <- 10
sigma0   <- 2
nu0      <- 1  # For SST, 1 = symmetrical, > 1 = right skewed, < 1 = left skewed
tau0     <- 10  # for TF2, nu = 10 means close to normal, the lower the number the heavier                    the tails
# for SST, large tau (e.g 7) thinner tails, close to normal. Small tau means                  fatter tails 

# Slopes (linear change rates)
beta_sigma <- 0.5   # scale increases by 0.5 over full range of t_scaled
beta_nu    <- 0.3   # skewness increases linearly
beta_tau   <- -1    # kurtosis decreases linearly

# changes with time
sigma_t <- pmax(0, sigma0 + beta_sigma * t_scaled)
nu_t    <- nu0    + beta_nu  * t_scaled
tau_t   <- tau0   + beta_tau * t_scaled

y <- rSST(n, mu = mu0, sigma = sigma0, nu = nu0, tau = tau0)

dat <- data.frame(
  t = t_index,
  t_scaled = t_scaled,
  y = y,
  sd_true = sigma0,
  mu_true = mu0,
  nu_true = nu0,
  tau_true = tau0,
  dist = "SST"
)


# faceted time series plot
ggplot(dat_all, aes(x = t)) +
  geom_line(
    aes(y = y),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  facet_wrap(~dist)

# density plots
ggplot(dat_all, aes(y = y)) +
  geom_density(
    aes(color = dist),
    linewidth = 0.8
  ) +
  labs(
    x = "Density",
    y = "Value"
  ) +
  theme_minimal(base_size = 14) +
  coord_flip()


SSTdat <- gamlss(y ~ 1, sigma.fo = ~ 1, nu.fo = ~ 1, tau.fo = ~ 1, family = SST(), data = dat, 
                 method = mixed(5, 200),
                 control = gamlss.control(n.cyc = 400, c.crit = 0.01, trace = FALSE)
)