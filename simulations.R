# =========================
# Sim 1: time-varying variance (robust, no purrr)
# =========================

library(gamlss)
library(gamlss.dist)
library(ggplot2)
library(dplyr)

set.seed(123)

# ---- 1) SIMULATE ----
n          <- 800
t_index    <- seq_len(n)
t_scaled   <- (t_index - 1)/(n - 1)   # 0..1
mu0        <- 0
sigma0     <- 2
amp        <- 0.7                     # strength of variance modulation

sigma_t <- sigma0 * (1 + amp * sin(2*pi*t_scaled) + 0.2 * t_scaled)
sigma_t <- pmax(sigma_t, 0.2)         # keep positive
y <- rnorm(n, mean = mu0, sd = sigma_t)

dat <- data.frame(
  t = t_index,
  t_scaled = t_scaled,
  y = y,
  sd_true = sigma_t,
  mu_true = mu0
)

value <- rnorm(n, mean = mu0, sd = sigma0)

# plot the sim 
ggplot(dat, aes(x = t)) +
  geom_line(
    aes(y = value),
    color = "darkred",
    linewidth = 0.8
  ) +
  labs(
    x = "Time",
    y = "Value"
  ) +
  theme_minimal(base_size = 14)

# ---- 2) FIT FULL-SERIES MODELS ----
# M0: constant variance
m0 <- gamlss(
  y ~ 1,
  sigma.fo = ~ 1,   # constant variance 
  family = NO(),
  data = dat,
  trace = FALSE
)

# M1: sigma varies smoothly with time
m1 <- gamlss(
  y ~ 1,                     # Sim 1: mean is constant
  sigma.fo = ~ pb(t_scaled),   # sigma(t)
  family = NO(),
  data = dat,
  trace = FALSE
)

aic_tbl <- data.frame(
  model = c("M0: const sigma", "M1: sigma(t)"),
  AIC   = c(GAIC(m0, k = 2), GAIC(m1, k = 2))
)
print(aic_tbl)

# ---- 3) ROLLING ONE-STEP-AHEAD EVALUATION ----
# Helper: CRPS for Normal (no extra pkg needed)
crps_norm_manual <- function(y, mean, sd) {
  # Gneiting & Raftery (2007) closed form
  if (sd <= 0) return(NA_real_)
  z <- (y - mean)/sd
  sd * (z*(2*pnorm(z)-1) + 2*dnorm(z) - 1/sqrt(pi))
}

k_start <- 200      # start after we have enough training data
k_end   <- n - 1

# Preallocate
K <- k_end - k_start + 1
logscore_M0 <- numeric(K); logscore_M1 <- numeric(K)
crps_M0     <- numeric(K); crps_M1     <- numeric(K)
cov95_M0    <- logical(K); cov95_M1    <- logical(K)

i <- 0
for (k in k_start:k_end) {
  i <- i + 1
  
  # Train/test split
  train <- dat[1:k, , drop = FALSE]
  test1 <- dat[k + 1L, , drop = FALSE]
  
  # Fit M0 safely
  m0_k <- try(
    gamlss(y ~ 1, sigma.fo = ~ 1, family = NO(), data = train, trace = FALSE),
    silent = TRUE
  )
  # Fit M1 safely
  m1_k <- try(
    gamlss(y ~ 1, sigma.fo = ~ pb(t_scaled), family = NO(), data = train, trace = FALSE),
    silent = TRUE
  )
  
  # If either fails (should be rare), carry NA and continue
  if (inherits(m0_k, "try-error") || inherits(m1_k, "try-error")) {
    logscore_M0[i] <- NA_real_; logscore_M1[i] <- NA_real_
    crps_M0[i]     <- NA_real_; crps_M1[i]     <- NA_real_
    cov95_M0[i]    <- NA;       cov95_M1[i]    <- NA
    next
  }
  
  # Predict params for the next point
  p0 <- predictAll(m0_k, newdata = test1, type = "response")
  p1 <- predictAll(m1_k, newdata = test1, type = "response")
  
  mu0_hat <- as.numeric(p0$mu);  sd0_hat <- as.numeric(p0$sigma)
  mu1_hat <- as.numeric(p1$mu);  sd1_hat <- as.numeric(p1$sigma)
  y_obs   <- test1$y
  
  # Log score (log predictive density under Normal)
  logscore_M0[i] <- dNO(y_obs, mu = mu0_hat, sigma = sd0_hat, log = TRUE)
  logscore_M1[i] <- dNO(y_obs, mu = mu1_hat, sigma = sd1_hat, log = TRUE)
  
  # CRPS
  crps_M0[i] <- crps_norm_manual(y_obs, mu0_hat, sd0_hat)
  crps_M1[i] <- crps_norm_manual(y_obs, mu1_hat, sd1_hat)
  
  # 95% coverage
  q0_lo <- qNO(0.025, mu = mu0_hat, sigma = sd0_hat)
  q0_hi <- qNO(0.975, mu = mu0_hat, sigma = sd0_hat)
  q1_lo <- qNO(0.025, mu = mu1_hat, sigma = sd1_hat)
  q1_hi <- qNO(0.975, mu = mu1_hat, sigma = sd1_hat)
  
  cov95_M0[i] <- (y_obs >= q0_lo) && (y_obs <= q0_hi)
  cov95_M1[i] <- (y_obs >= q1_lo) && (y_obs <= q1_hi)
}

summary_tbl <- data.frame(
  model          = c("M0: const sigma", "M1: sigma(t)"),
  mean_log_score = c(mean(logscore_M0, na.rm = TRUE), mean(logscore_M1, na.rm = TRUE)), 
  mean_crps      = c(mean(crps_M0,     na.rm = TRUE), mean(crps_M1,     na.rm = TRUE)),
  cover95_rate   = c(mean(cov95_M0,    na.rm = TRUE), mean(cov95_M1,    na.rm = TRUE))
)
print(summary_tbl)

# WHAT DO THEY MEAN?
## mean_log_score = Average log predictive density at the observed values. Closer to 0 (less negative) 
### is better, means that the model put more probability mass on what actually happened. So here, M1 is 
### closer to 0, therefore that model gives systematically higher likelihood to observed outcomes.
## mean_crps = Continuous Ranked Probability Score, it measure the accuracy of the whole predictive distribution 
### (mean + spread), and lower is better. Here, M1 is lower, hence its predictive distribution are sharper and 
### better aligned with reality 
## cover95_rate = its the coverage of 95% interval, looking at how often the realized value fell inside the 
### 95% predictive interval. we want it to be ~ 0.95 (predictive intervals are well calibrated). Here, M0 is ~ 97%, 
### while M1 is ~ 93%. M0 is playing is too safe by predicting larger interval that rarely miss, but waste sharpness 
### (overestimating uncertainty). M1 is too bold by predicting narrower interval that miss a bit more often than they 
### should. Both are close to ~95% so they're well calibrated 

## -- Therefore overall, M1 is the better model -- ##


# ---- 4) VISUAL CHECKS ----
# True vs fitted sigma (full-series M1)
sig_hat <- as.numeric(fitted(m1, "sigma"))

print(ggplot(dat, aes(t)) +
        geom_line(aes(y = sd_true)) +
        geom_line(aes(y = sig_hat), colour = "red") +
        labs(title = "Sim 1: SD true (black) vs fitted sigma(t) (red)",
             x = "time", y = "SD"))

# PIT histogram for M1 (quick calibration feel)
pa_full <- predictAll(m1, type = "response")
u <- pNO(dat$y, mu = as.numeric(pa_full$mu), sigma = as.numeric(pa_full$sigma))

print(ggplot(data.frame(u = u), aes(u)) +
        geom_histogram(bins = 20) +
        labs(title = "PIT histogram (M1)", x = "u", y = "count"))




# =========================
# Sim 2: time-varying mean + variance
# =========================

set.seed(456)

# ---- 1) SIMULATE ----
n          <- 800
t_index    <- seq_len(n)
t_scaled   <- (t_index - 1)/(n - 1)

# Mean varies smoothly (sinusoidal + drift)
mu_t     <- 2 * sin(2*pi*t_scaled) + 0.5 * t_scaled

# Variance varies with time
sigma_t  <- 1.5 + 0.5*sin(4*pi*t_scaled) + 0.3*t_scaled
sigma_t  <- pmax(sigma_t, 0.2)  # keep >0

y <- rnorm(n, mean = mu_t, sd = sigma_t)

dat <- data.frame(
  t = t_index,
  t_scaled = t_scaled,
  y = y,
  mu_true = mu_t,
  sd_true = sigma_t
)

# plot the sim 
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

# ---- 2) FIT MODELS ON FULL SERIES ----
# M0: constant mean + constant variance
m0 <- gamlss(y ~ 1, sigma.fo = ~ 1, family = NO(), data = dat, trace = FALSE)

# M1: time-varying mean only
m1 <- gamlss(y ~ pb(t_scaled), sigma.fo = ~ 1, family = NO(), data = dat, trace = FALSE)

# M2: time-varying mean + time-varying variance
m2 <- gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled), family = NO(), data = dat, trace = FALSE)

# M3: constant mean + time-varying variance 
m3 <- gamlss(y ~ 1, sigma.fo = ~ pb(t_scaled), family = NO(), data = dat, trace = FALSE)

aic_tbl <- data.frame(
  model = c("M0: const mu,sigma", "M1: mu(t), const sigma", "M2: mu(t), sigma(t)", "M3: const mu, sigma(t)"),
  AIC   = c(GAIC(m0, k = 2), GAIC(m1, k = 2), GAIC(m2, k = 2), GAIC(m3, k = 2))
)
print(aic_tbl)

# ---- 3) ROLLING ONE-STEP-AHEAD EVALUATION ----
crps_norm_manual <- function(y, mean, sd) {
  if (sd <= 0) return(NA_real_)
  z <- (y - mean)/sd
  sd * (z*(2*pnorm(z)-1) + 2*dnorm(z) - 1/sqrt(pi))
}

k_start <- 200
k_end   <- n - 1
K <- k_end - k_start + 1

results <- data.frame(
  k = integer(),
  model = character(),
  log_score = numeric(),
  crps = numeric(),
  cover95 = logical(),
  stringsAsFactors = FALSE
)

for (k in k_start:k_end) {
  train <- dat[1:k, ]
  test1 <- dat[k+1, , drop=FALSE]
  
  # Fit models on training data
  m0_k <- try(gamlss(y ~ 1, sigma.fo = ~ 1, family = NO(), data = train, trace = FALSE), silent = TRUE)
  m1_k <- try(gamlss(y ~ pb(t_scaled), sigma.fo = ~ 1, family = NO(), data = train, trace = FALSE), silent = TRUE)
  m2_k <- try(gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled), family = NO(), data = train, trace = FALSE), silent = TRUE)
  m3_k <- try(gamlss(y ~ 1, sigma.fo = ~ pb(t_scaled), family = NO(), data = train, trace = FALSE), silent = TRUE)
  
  
  fits <- list(M0 = m0_k, M1 = m1_k, M2 = m2_k, M3 = m3_k)
  
  for (mod in names(fits)) {
    fit <- fits[[mod]]
    if (inherits(fit, "try-error")) next
    
    pa <- suppressWarnings(predictAll(fit, newdata = test1, type = "response"))
    mu_hat <- as.numeric(pa$mu); sd_hat <- as.numeric(pa$sigma)
    y_obs  <- test1$y
    
    ls  <- dNO(y_obs, mu = mu_hat, sigma = sd_hat, log = TRUE)
    cr  <- crps_norm_manual(y_obs, mu_hat, sd_hat)
    qlo <- qNO(0.025, mu = mu_hat, sigma = sd_hat)
    qhi <- qNO(0.975, mu = mu_hat, sigma = sd_hat)
    cov <- (y_obs >= qlo) && (y_obs <= qhi)
    
    results <- rbind(results, data.frame(
      k = k, model = mod, log_score = ls, crps = cr, cover95 = cov, stringsAsFactors = FALSE
    ))
  }
}

summary_tbl <- results %>%
  group_by(model) %>%
  summarise(
    mean_log_score = mean(log_score, na.rm = TRUE),
    mean_crps      = mean(crps, na.rm = TRUE),
    cover95_rate   = mean(cover95, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_tbl)

# WHAT DO THEY MEAN?
## mean_log_score = Average log predictive density at the observed values. Closer to 0 (less negative) 
### is better, means that the model put more probability mass on what actually happened. So here, M2 is 
### closer to 0, therefore that model gives systematically higher likelihood to observed outcomes.
## mean_crps = Continuous Ranked Probability Score, it measure the accuracy of the whole predictive distribution 
### (mean + spread), and lower is better. Here, M2 is lower, hence its predictive distribution are sharper and 
### better aligned with reality 
## cover95_rate = its the coverage of 95% interval, looking at how often the realized value fell inside the 
### 95% predictive interval. we want it to be ~ 0.95 (predictive intervals are well calibrated). Here, M1 and M3
### are ~95%, means that they are well calibrated (the predictive intervals are the right 'width'). M0 and M2 are 
### under-covering with ~93%. The intervals are slightly too tight (too confident). Makes sense for M0 has it ignore 
### both mean and variance changes. M2 is interesting cause the smooth penalties in pb() can make it slightly under-dispersed 
### in finite samples


# ---- 4) VISUAL CHECK ----
# True vs fitted mean (M2)
mu_hat <- fitted(m2, "mu")

ggplot(dat, aes(t)) +
  geom_line(aes(y = mu_true), colour="black") +
  geom_line(aes(y = mu_hat), colour="red") +
  labs(title="Sim 2: true mean (black) vs fitted mu(t) (red)", y="mean")

# True vs fitted sigma (M2)
sig_hat <- fitted(m2, "sigma")

ggplot(dat, aes(t)) +
  geom_line(aes(y = sd_true), colour="black") +
  geom_line(aes(y = sig_hat), colour="blue") +
  labs(title="Sim 2: true sd (black) vs fitted sigma(t) (blue)", y="sd")




# =========================
# Sim 3 (robust): mu(t), sigma(t) vary; skew & kurtosis constant
# Auto-detect skew-t generator; fall back to t if needed
# =========================
set.seed(1234)

# ---- 1) SIMULATE ----
n        <- 800
t_idx    <- seq_len(n)
t_scaled <- (t_idx - 1)/(n - 1)

mu_t    <- 1.2*sin(2*pi*t_scaled) + 0.5*t_scaled
sigma_t <- 1.2 + 0.4*sin(4*pi*t_scaled + 0.6) + 0.25*t_scaled
sigma_t <- pmax(sigma_t, 0.15)

nu_const  <- 2     # skew parameter (constant)
tau_const <- 7     # tail/heaviness (constant)

y <- rSHASHo(length(t_scaled), mu = mu_t, sigma = sigma_t, nu = nu_const, tau = tau_const)
y <- rTF(length(t_scaled), mu = mu_t, sigma = sigma_t, nu = nu_const)


dat <- data.frame(
  t = t_idx, t_scaled = t_scaled, y = y,
  mu_true = mu_t, sd_true = sigma_t, nu_true = nu_const, tau_true = tau_const
)

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


# ---- 2) FULL-SERIES FITS ----
m0 <- gamlss(y ~ 1, sigma.fo = ~ 1, family = NO(), data = dat, trace = FALSE)
m1 <- gamlss(y ~ pb(t_scaled), sigma.fo = ~ 1, family = NO(), data = dat, trace = FALSE)
m2 <- gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled), family = NO(), data = dat, trace = FALSE)
m3 <- gamlss(y ~ 1, sigma.fo = ~ pb(t_scaled), family = NO(), data = dat, trace = FALSE)
m4 <- gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled), nu.fo = ~ 1, tau.fo = ~ 1, family = SHASHo(), 
             data = dat, trace = FALSE)
m5 <- gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled), nu.fo = ~ 1, family = TF(),
  data = dat, trace = FALSE)


aic_tbl <- data.frame(
  model = c("M0: NO const mu,sigma", 
            "M1: NO mu(t), const sigma",
            "M2: NO mu(t), sigma(t)",
            "M3: NO const mu, sigma(t)",
            "M4: SHASHo mu(t), sigma(t), const nu, const tau",
            "M5: TF mu(t), sigma(t), const nu"),
  AIC   = c(GAIC(m0, k = 2), GAIC(m1, k = 2), GAIC(m2,k = 2), GAIC(m3,k = 2), GAIC(m4, k = 2), GAIC(m5, k =2))
)
print(aic_tbl)


# ---- 3) ROLLING ONE-STEP-AHEAD EVALUATION ----
crps_mc <- function(y, samp) {
  if (length(samp) < 2) return(NA_real_)
  e1 <- mean(abs(samp - y))
  e2 <- mean(abs(samp - sample(samp, length(samp), replace = TRUE)))
  e1 - 0.5*e2
}

draw_pred <- function(fam, pars, n_samp = 400L) {
  if (fam == "NO") {
    rNO(n_samp, mu = pars$mu, sigma = pars$sigma)
  } else if (fam == "ST") {
    rST(n_samp, mu = pars$mu, sigma = pars$sigma, nu = pars$nu, tau = pars$tau)
  } else stop("Unknown family")
}

score_one <- function(fit, newrow, fam) {
  pa <- suppressWarnings(predictAll(fit, newdata = newrow, type = "response"))
  y_obs <- newrow$y
  
  if (fam == "NO") {
    mu <- as.numeric(pa$mu); sg <- as.numeric(pa$sigma)
    ls <- dNO(y_obs, mu = mu, sigma = sg, log = TRUE)
    q  <- qNO(c(0.025, 0.975), mu = mu, sigma = sg)
    samp <- draw_pred("NO", list(mu = mu, sigma = sg))
  } else { # ST
    mu <- as.numeric(pa$mu); sg <- as.numeric(pa$sigma)
    nu <- as.numeric(pa$nu); ta <- as.numeric(pa$tau)
    ls <- dST(y_obs, mu = mu, sigma = sg, nu = nu, tau = ta, log = TRUE)
    q  <- qST(c(0.025, 0.975), mu = mu, sigma = sg, nu = nu, tau = ta)
    samp <- draw_pred("ST", list(mu = mu, sigma = sg, nu = nu, tau = ta))
  }
  
  list(log_score = ls,
       crps      = crps_mc(y_obs, samp),
       cover95   = (y_obs >= q[1]) && (y_obs <= q[2]))
}

k_start <- 200
k_end   <- n - 1

res <- data.frame(k = integer(), model = character(),
                  log_score = numeric(), crps = numeric(), cover95 = logical(),
                  stringsAsFactors = FALSE)

for (k in k_start:k_end) {
  train <- dat[1:k, ]
  test1 <- dat[k+1, , drop = FALSE]
  
  fits <- list(
    M0 = try(gamlss(y ~ 1, sigma.fo = ~ 1, family = NO(), data = train, trace = FALSE), silent = TRUE),
    M1 = try(gamlss(y ~ pb(t_scaled), sigma.fo = ~ 1, family = NO(), data = train, trace = FALSE), silent = TRUE),
    M2 = try(gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled),  family = NO(), data = train, trace = FALSE), silent = TRUE),
    M3 = try(gamlss(y ~ 1, sigma.fo = ~ pb(t_scaled),  family = NO(), data = train, trace = FALSE), silent = TRUE),
    M4 = try(gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled),  nu.fo = ~ 1,   tau.fo = ~ 1,   family = SHASHo(), data = train, trace = FALSE), silent = TRUE),
    M5 = try(gamlss(y ~ pb(t_scaled), sigma.fo = ~ pb(t_scaled),  nu.fo = ~ pb(t_scaled), tau.fo = ~ 1, family = TF(), data = train, trace = FALSE), silent = TRUE),
)
  
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (inherits(fit, "try-error")) next
    fam <- ifelse(nm %in% c("M0","M1","M2","M3"), "NO", "ST")
    sc  <- score_one(fit, test1, fam)
    res <- rbind(res, data.frame(k = k, model = nm,
                                 log_score = sc$log_score,
                                 crps = sc$crps,
                                 cover95 = sc$cover95,
                                 stringsAsFactors = FALSE))
  }
}

summary_tbl <- res %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(
    mean_log_score = mean(log_score, na.rm = TRUE),
    mean_crps      = mean(crps,      na.rm = TRUE),
    cover95_rate   = mean(cover95,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(mean_log_score))
print(summary_tbl)

# ---- 4) VISUAL CHECKS (full-series fits) ----
# True vs fitted mu, sigma for the truth-like model (M4)
mu_hat  <- as.numeric(fitted(m4, "mu"))
sig_hat <- as.numeric(fitted(m4, "sigma"))

print(ggplot(dat, aes(t)) +
        geom_line(aes(y = mu_true)) +
        geom_line(aes(y = mu_hat), colour = "red") +
        labs(title = "Sim 3: true mu (black) vs fitted mu(t) in M4 (red)", y = "mu"))

print(ggplot(dat, aes(t)) +
        geom_line(aes(y = sd_true)) +
        geom_line(aes(y = sig_hat), colour = "blue") +
        labs(title = "Sim 3: true sd (black) vs fitted sigma(t) in M4 (blue)", y = "sd"))

?SN1
?TF
?SHASHo
?SST




################ SIMULATION FOR STATIC DISTRIBUTION #####################


## SIM 1: symmetrical and mesokurtotic --> NO()

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 30 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015   # scale increases by 0.5 over full range of t_scaled

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <-  2 + beta_sigma * t_index



# # observation randomly taken from normal distribution 
y <- rNO(t, mu = mu_t, sigma = sigma_t)
y <- trend + seasonal_effect + y

# create dataset
dat <- data.frame(
  t = t_index,
  month = factor(month),
  y = y,
  sd_true = sigma_t,
  mu_true = mu_t,
  dist = "NO"
)

# plot them to see
ggplot(dat, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 
  

# look at density plot
ggplot(dat, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat$y), linetype = "dashed", col = "grey") +
  theme_minimal()

# histogram and non-parametric density estimate
sim1_p1 <- ggplot(dat, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  #xlim(x = c(0,60)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

mean(dat$y)
median(dat$y)
sd(dat$y)
skewness(dat$y)
kurtosis(dat$y)

### fit models of different family distribution ###

# NO
sim1_static_NO <- gamlss(y ~ 1, family = NO(), data = dat) # two-parameter
# NO2
sim1_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat)
# GU
sim1_static_GU <- gamlss(y ~ 1, family = GU(), data = dat)
# LO
sim1_static_LO <- gamlss(y ~ 1, family = LO(), data = dat)
# RG
sim1_static_RG <- gamlss(y ~ 1, family = RG(), data = dat)
# exGAUS
sim1_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim1_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat)
# PE
sim1_static_PE <- gamlss(y ~ 1, family = PE(), data = dat,
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim1_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim1_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat)
# SN2
sim1_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim1_static_TF <- gamlss(y ~ 1, family = TF(), data = dat)
# TF2
sim1_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim1_static_GT <- gamlss(y ~ 1, family = GT(), data = dat, 
                       method = mixed(10,200),
                       control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim1_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat, 
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim1_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat,
                         mu.start = mean(dat$y), sigma.start = sd(dat$y),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim1_static_NET <- gamlss(y ~ 1, family = NET(), data = dat)
# SHASH
sim1_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim1_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim1_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim1_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim1_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat,
                         mu.start = mean(dat$y), sigma.start = sd(dat$y), 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim1_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat,
                         mu.start = mean(dat$y), sigma.start = sd(dat$y),
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim1_static_SST <- gamlss(y ~ 1, family = SST(), data = dat,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim1_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim1_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim1_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim1_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat,
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim1_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat,
                        mu.start = mean(dat$y), sigma.start = sd(dat$y),
                        method = mixed(10,200),
                        control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim1_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim1_static_NO, k = 2),GAIC(sim1_static_NO2, k = 2), GAIC(sim1_static_GU, k = 2), GAIC(sim1_static_LO, k = 2), 
          GAIC(sim1_static_RG, k = 2), GAIC(sim1_static_exGAUS, k = 2), GAIC(sim1_static_NOF, k = 2), GAIC(sim1_static_PE, k = 2), 
          GAIC(sim1_static_PE2, k = 2),GAIC(sim1_static_SN1, k = 2), GAIC(sim1_static_TF, k = 2), GAIC(sim1_static_TF2, k = 2),
          GAIC(sim1_static_GT, k = 2),GAIC(sim1_static_JSU, k = 2), GAIC(sim1_static_JSUo, k = 2), GAIC(sim1_static_NET, k = 2),
          GAIC(sim1_static_SHASH, k = 2),GAIC(sim1_static_SHASHo, k = 2), GAIC(sim1_static_SHASHo2, k = 2), GAIC(sim1_static_SEP2, k = 2), 
          GAIC(sim1_static_SEP3, k = 2), GAIC(sim1_static_SEP4, k = 2), GAIC(sim1_static_SST, k = 2), GAIC(sim1_static_ST1, k = 2),
          GAIC(sim1_static_ST2, k = 2),GAIC(sim1_static_ST3, k = 2), GAIC(sim1_static_ST4, k = 2), GAIC(sim1_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
                "lepto","meso", "both", "both", "lepto", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto","lepto" ,"lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim1_pAIC_1a <- ggplot(sim1_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


summary(sim1_static_SHASH)
summary(sim1b_static_SHASH)

## SIM 1a: symmetrical, mesokurtotic, constant var --> NO() ========================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015   # scale increases by 0.5 over full range of t_scaled

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma0 <- 2


# # observation randomly taken from normal distribution 
y1a <- rNO(t, mu = mu_t, sigma = sigma0)
y1a <- trend + seasonal_effect + y1a

# create dataset
dat1a <- data.frame(
  t = t_index,
  month = factor(month),
  y = y1a,
  sd_true = sigma0,
  mu_true = mu_t,
  dist = "NO"
)

# plot them to see
ggplot(dat1a, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat1a, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat1a$y), linetype = "dashed", col = "grey") +
  theme_minimal()

# histogram and non-parametric density estimate
sim1a_p1 <- ggplot(dat1a, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim1a_static_NO <- gamlss(y ~ 1, family = NO(), data = dat1a) # two-parameter
# NO2
sim1a_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat1a)
# GU
sim1a_static_GU <- gamlss(y ~ 1, family = GU(), data = dat1a)
# LO
sim1a_static_LO <- gamlss(y ~ 1, family = LO(), data = dat1a)
# RG
sim1a_static_RG <- gamlss(y ~ 1, family = RG(), data = dat1a)
# exGAUS
sim1a_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat1a, # error
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim1a_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat1a)
# PE
sim1a_static_PE <- gamlss(y ~ 1, family = PE(), data = dat1a,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim1a_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat1a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim1a_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat1a)
# SN2
sim1a_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat1a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim1a_static_TF <- gamlss(y ~ 1, family = TF(), data = dat1a)
# TF2
sim1a_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat1a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim1a_static_GT <- gamlss(y ~ 1, family = GT(), data = dat1a, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim1a_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat1a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim1a_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat1a,
                           mu.start = mean(dat1a$y), sigma.start = sd(dat1a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim1a_static_NET <- gamlss(y ~ 1, family = NET(), data = dat1a)
# SHASH
sim1a_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat1a,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim1a_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat1a,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim1a_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat1a,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim1a_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat1a,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim1a_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat1a,
                           mu.start = mean(dat1a$y), sigma.start = sd(dat1a$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim1a_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat1a,
                           mu.start = mean(dat1a$y), sigma.start = sd(dat1a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim1a_static_SST <- gamlss(y ~ 1, family = SST(), data = dat1a, # not working 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim1a_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat1a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim1a_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat1a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim1a_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat1a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim1a_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat1a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim1a_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat1a,
                          mu.start = mean(dat1a$y), sigma.start = sd(dat1a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim1a_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim1a_static_NO, k = 2),GAIC(sim1a_static_NO2, k = 2), GAIC(sim1a_static_GU, k = 2), GAIC(sim1a_static_LO, k = 2), 
          GAIC(sim1a_static_RG, k = 2), GAIC(sim1a_static_NOF, k = 2), GAIC(sim1a_static_PE, k = 2), 
          GAIC(sim1a_static_PE2, k = 2),GAIC(sim1a_static_SN1, k = 2), GAIC(sim1a_static_TF, k = 2), GAIC(sim1a_static_TF2, k = 2),
          GAIC(sim1a_static_GT, k = 2),GAIC(sim1a_static_JSU, k = 2), GAIC(sim1a_static_JSUo, k = 2), GAIC(sim1a_static_NET, k = 2),
          GAIC(sim1a_static_SHASH, k = 2),GAIC(sim1a_static_SHASHo, k = 2), GAIC(sim1a_static_SHASHo2, k = 2), GAIC(sim1a_static_SEP2, k = 2), 
          GAIC(sim1a_static_SEP3, k = 2), GAIC(sim1a_static_SEP4, k = 2), GAIC(sim1a_static_ST1, k = 2),
          GAIC(sim1a_static_ST2, k = 2),GAIC(sim1a_static_ST3, k = 2), GAIC(sim1a_static_ST4, k = 2), GAIC(sim1a_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "sym", "sym", "sym", "both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "meso", "both", "both", "lepto", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim1a_pAIC_1a <- ggplot(sim1a_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso,sym,var0",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# compare
sim1_pAIC_1a + sim1a_pAIC_1a



## SIM 1b: symmetrical, mesokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015   # scale increases by 0.5 over full range of t_scaled

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <-  2 + beta_sigma * t_index



# # observation randomly taken from normal distribution 
y1b <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 1, tau = 1000)
y1b <- trend + seasonal_effect + y1b

# create dataset
dat1b <- data.frame(
  t = t_index,
  month = factor(month),
  y = y1b,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 1,
  tau_true = 1000,
  dist = "SST"
)

# plot them to see
ggplot(dat1b, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat1b, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat1b$y), linetype = "dashed", col = "grey") +
  theme_minimal()

# histogram and non-parametric density estimate
sim1b_p1 <- ggplot(dat1b, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

mean(dat1b$y)
median(dat1b$y)
sd(dat1b$y)
skewness(dat1b$y)
kurtosis(dat1b$y)

### fit models of different family distribution ###

# NO
sim1b_static_NO <- gamlss(y ~ 1, family = NO(), data = dat1b) # two-parameter
# NO2
sim1b_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat1b)
# GU
sim1b_static_GU <- gamlss(y ~ 1, family = GU(), data = dat1b)
# LO
sim1b_static_LO <- gamlss(y ~ 1, family = LO(), data = dat1b)
# RG
sim1b_static_RG <- gamlss(y ~ 1, family = RG(), data = dat1b)
# exGAUS
sim1b_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat1b, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim1b_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat1b)
# PE
sim1b_static_PE <- gamlss(y ~ 1, family = PE(), data = dat1b,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim1b_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat1b, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim1b_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat1b)
# SN2
sim1b_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat1b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim1b_static_TF <- gamlss(y ~ 1, family = TF(), data = dat1b)
# TF2
sim1b_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat1b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim1b_static_GT <- gamlss(y ~ 1, family = GT(), data = dat1b, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim1b_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat1b, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim1b_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat1b,
                           mu.start = mean(dat1b$y), sigma.start = sd(dat1b$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim1b_static_NET <- gamlss(y ~ 1, family = NET(), data = dat1b)
# SHASH
sim1b_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat1b,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim1b_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat1b,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim1b_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat1b,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim1b_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat1b,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim1b_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat1b,
                           mu.start = mean(dat1b$y), sigma.start = sd(dat1b$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim1b_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat1b,
                           mu.start = mean(dat1b$y), sigma.start = sd(dat1b$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim1b_static_SST <- gamlss(y ~ 1, family = SST(), data = dat1b, # not working
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim1b_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat1b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim1b_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat1b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim1b_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat1b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim1b_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat1b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim1b_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat1b,
                          mu.start = mean(dat1b$y), sigma.start = sd(dat1b$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim1b_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim1b_static_NO, k = 2),GAIC(sim1b_static_NO2, k = 2), GAIC(sim1b_static_GU, k = 2), GAIC(sim1b_static_LO, k = 2), 
          GAIC(sim1b_static_RG, k = 2), GAIC(sim1b_static_exGAUS, k = 2), GAIC(sim1b_static_NOF, k = 2), GAIC(sim1b_static_PE, k = 2), 
          GAIC(sim1b_static_PE2, k = 2),GAIC(sim1b_static_SN1, k = 2), GAIC(sim1b_static_TF, k = 2), GAIC(sim1b_static_TF2, k = 2),
          GAIC(sim1b_static_GT, k = 2),GAIC(sim1b_static_JSU, k = 2), GAIC(sim1b_static_JSUo, k = 2), GAIC(sim1b_static_NET, k = 2),
          GAIC(sim1b_static_SHASH, k = 2),GAIC(sim1b_static_SHASHo, k = 2), GAIC(sim1b_static_SHASHo2, k = 2), GAIC(sim1b_static_SEP2, k = 2), 
          GAIC(sim1b_static_SEP3, k = 2), GAIC(sim1b_static_SEP4, k = 2), GAIC(sim1b_static_ST1, k = 2),
          GAIC(sim1b_static_ST2, k = 2),GAIC(sim1b_static_ST3, k = 2), GAIC(sim1b_static_ST4, k = 2), GAIC(sim1b_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both", "both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto" ,"lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim1b_pAIC_1a <- ggplot(sim1b_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso,sym,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


sim1_pAIC_1a + sim1a_pAIC_1a + sim1b_pAIC_1a

# look at paramter summary 
summary(sim1b_static_NO)
summary(sim1b_static_SHASH)



## SIM 2: +ve skewness, mesokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015   # scale increases by 0.5 over full range of t_scaled

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)



# # observation randomly taken from normal distribution 
y2 <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 3, tau = 1000)
y2 <- trend + seasonal_effect + y2

# create dataset
dat2 <- data.frame(
  t = t_index,
  month = factor(month),
  y = y2,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 3,
  tau_true = 1000,
  dist = "SST"
)

# plot them to see
ggplot(dat2, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat2, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim2_p1 <- ggplot(dat2, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim2_static_NO <- gamlss(y ~ 1, family = NO(), data = dat2) # two-parameter
# NO2
sim2_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat2)
# GU
sim2_static_GU <- gamlss(y ~ 1, family = GU(), data = dat2)
# LO
sim2_static_LO <- gamlss(y ~ 1, family = LO(), data = dat2)
# RG
sim2_static_RG <- gamlss(y ~ 1, family = RG(), data = dat2)
# exGAUS
sim2_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat2, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim2_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat2)
# PE
sim2_static_PE <- gamlss(y ~ 1, family = PE(), data = dat2,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim2_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat2, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim2_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat2)
# SN2
sim2_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim2_static_TF <- gamlss(y ~ 1, family = TF(), data = dat2)
# TF2
sim2_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim2_static_GT <- gamlss(y ~ 1, family = GT(), data = dat2, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim2_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat2, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim2_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat2,
                            mu.start = mean(dat2$y), sigma.start = sd(dat2$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim2_static_NET <- gamlss(y ~ 1, family = NET(), data = dat2)
# SHASH
sim2_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat2,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim2_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat2,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim2_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat2,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim2_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat2, # not working
                           u.start = mean(dat2$y), sigma.start = sd(dat2$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim2_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat2,
                            mu.start = mean(dat2$y), sigma.start = sd(dat2$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim2_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat2,  # convergence issue
                            mu.start = mean(dat2$y), sigma.start = sd(dat2$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim2_static_SST <- gamlss(y ~ 1, family = SST(), data = dat2, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim2_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim2_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat2,
                          mu.start = mean(dat2$y), sigma.start = sd(dat2$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim2_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim2_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat2,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim2_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat2,
                           mu.start = mean(dat2$y), sigma.start = sd(dat2$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim2_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN1", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP3","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim2_static_NO, k = 2),GAIC(sim2_static_NO2, k = 2), GAIC(sim2_static_GU, k = 2), GAIC(sim2_static_LO, k = 2), 
          GAIC(sim2_static_RG, k = 2), GAIC(sim2_static_exGAUS, k = 2), GAIC(sim2_static_NOF, k = 2), GAIC(sim2_static_PE, k = 2), 
          GAIC(sim2_static_PE2, k = 2),GAIC(sim2_static_SN1, k = 2), GAIC(sim2_static_TF, k = 2), GAIC(sim2_static_TF2, k = 2),
          GAIC(sim2_static_GT, k = 2),GAIC(sim2_static_JSU, k = 2), GAIC(sim2_static_JSUo, k = 2), GAIC(sim2_static_NET, k = 2),
          GAIC(sim2_static_SHASH, k = 2),GAIC(sim2_static_SHASHo, k = 2), GAIC(sim2_static_SHASHo2, k = 2),  
          GAIC(sim2_static_SEP3, k = 2), GAIC(sim2_static_SST, k = 2),GAIC(sim2_static_ST1, k = 2),
          GAIC(sim2_static_ST2, k = 2),GAIC(sim2_static_ST3, k = 2), GAIC(sim2_static_ST4, k = 2), GAIC(sim2_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "both", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","lepto" ,"lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim2_pAIC_1a <- ggplot(sim2_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso,+ve skew,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


## SIM 2a: +ve skewness(time), mesokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015   # scale increases by 0.5 over full range of t_scaled
beta_nu <- 0.0005

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)
nu_t <- pmax(0, 1 + beta_nu * t_index)


# # observation randomly taken from normal distribution 
y2a <- rSST(t, mu = mu_t, sigma = sigma_t, nu = nu_t, tau = 1000)
y2a <- trend + seasonal_effect + y2a

# create dataset
dat2a <- data.frame(
  t = t_index,
  month = factor(month),
  y = y2a,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = nu_t,
  tau_true = 1000,
  dist = "SST"
)

# plot them to see
ggplot(dat2a, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat2a, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim2a_p1 <- ggplot(dat2a, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim2a_static_NO <- gamlss(y ~ 1, family = NO(), data = dat2a) # two-parameter
# NO2
sim2a_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat2a)
# GU
sim2a_static_GU <- gamlss(y ~ 1, family = GU(), data = dat2a)
# LO
sim2a_static_LO <- gamlss(y ~ 1, family = LO(), data = dat2a)
# RG
sim2a_static_RG <- gamlss(y ~ 1, family = RG(), data = dat2a)
# exGAUS
sim2a_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat2a, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim2a_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat2a)
# PE
sim2a_static_PE <- gamlss(y ~ 1, family = PE(), data = dat2a,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim2a_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat2a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim2a_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat2a) # error
# SN2
sim2a_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat2a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim2a_static_TF <- gamlss(y ~ 1, family = TF(), data = dat2a)
# TF2
sim2a_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat2a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim2a_static_GT <- gamlss(y ~ 1, family = GT(), data = dat2a, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim2a_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat2a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim2a_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat2a,
                           mu.start = mean(dat2a$y), sigma.start = sd(dat2a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim2a_static_NET <- gamlss(y ~ 1, family = NET(), data = dat2a)
# SHASH
sim2a_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat2a,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim2a_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat2a,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim2a_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat2a,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim2a_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat2a, # not working
                           u.start = mean(dat2a$y), sigma.start = sd(dat2a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim2a_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat2a,
                           mu.start = mean(dat2a$y), sigma.start = sd(dat2a$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim2a_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat2a,  # convergence issue
                           mu.start = mean(dat2a$y), sigma.start = sd(dat2a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim2a_static_SST <- gamlss(y ~ 1, family = SST(), data = dat2a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim2a_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat2a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim2a_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat2a,
                          mu.start = mean(dat2a$y), sigma.start = sd(dat2a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim2a_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat2a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim2a_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat2a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim2a_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat2a,
                          mu.start = mean(dat2a$y), sigma.start = sd(dat2a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim2a_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP3","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim2a_static_NO, k = 2),GAIC(sim2a_static_NO2, k = 2), GAIC(sim2a_static_GU, k = 2), GAIC(sim2a_static_LO, k = 2), 
          GAIC(sim2a_static_RG, k = 2), GAIC(sim2a_static_exGAUS, k = 2), GAIC(sim2a_static_NOF, k = 2), GAIC(sim2a_static_PE, k = 2), 
          GAIC(sim2a_static_PE2, k = 2),GAIC(sim2a_static_TF, k = 2), GAIC(sim2a_static_TF2, k = 2),
          GAIC(sim2a_static_GT, k = 2),GAIC(sim2a_static_JSU, k = 2), GAIC(sim2a_static_JSUo, k = 2), GAIC(sim2a_static_NET, k = 2),
          GAIC(sim2a_static_SHASH, k = 2),GAIC(sim2a_static_SHASHo, k = 2), GAIC(sim2a_static_SHASHo2, k = 2),  
          GAIC(sim2a_static_SEP3, k = 2), GAIC(sim2a_static_SST, k = 2),GAIC(sim2a_static_ST1, k = 2),
          GAIC(sim2a_static_ST2, k = 2),GAIC(sim2a_static_ST3, k = 2), GAIC(sim2a_static_ST4, k = 2), GAIC(sim2a_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4","4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym", "sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","lepto" ,"lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim2a_pAIC_1a <- ggplot(sim2a_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso,+ve skew(time),var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# compare
sim2_pAIC_1a + sim2a_pAIC_1a



## SIM 3: -ve skewness, mesokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015   # scale increases by 0.5 over full range of t_scaled


# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)


# # observation randomly taken from normal distribution 
y3 <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 0.3, tau = 1000)
y3 <- trend + seasonal_effect + y3

# create dataset
dat3 <- data.frame(
  t = t_index,
  month = factor(month),
  y = y3,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 0.3,
  tau_true = 1000,
  dist = "SST"
)

# plot them to see
ggplot(dat3, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat3, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim3_p1 <- ggplot(dat3, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim3_static_NO <- gamlss(y ~ 1, family = NO(), data = dat3) # two-parameter
# NO2
sim3_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat3)
# GU
sim3_static_GU <- gamlss(y ~ 1, family = GU(), data = dat3)
# LO
sim3_static_LO <- gamlss(y ~ 1, family = LO(), data = dat3)
# RG
sim3_static_RG <- gamlss(y ~ 1, family = RG(), data = dat3)
# exGAUS
sim3_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat3, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim3_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat3)
# PE
sim3_static_PE <- gamlss(y ~ 1, family = PE(), data = dat3,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim3_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat3, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim3_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat3) 
# SN2
sim3_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat3,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim3_static_TF <- gamlss(y ~ 1, family = TF(), data = dat3)
# TF2
sim3_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat3,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim3_static_GT <- gamlss(y ~ 1, family = GT(), data = dat3, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim3_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat3, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim3_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat3,
                            mu.start = mean(dat3$y), sigma.start = sd(dat3$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim3_static_NET <- gamlss(y ~ 1, family = NET(), data = dat3)
# SHASH
sim3_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat3,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim3_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat3,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim3_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat3,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim3_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat3, 
                            u.start = mean(dat3$y), sigma.start = sd(dat3$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim3_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat3,
                            mu.start = mean(dat3$y), sigma.start = sd(dat3$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim3_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat3,  
                            mu.start = mean(dat3$y), sigma.start = sd(dat3$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim3_static_SST <- gamlss(y ~ 1, family = SST(), data = dat3, # not working
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim3_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat3,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim3_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat3,
                           mu.start = mean(dat3$y), sigma.start = sd(dat3$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim3_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat3,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim3_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat3,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim3_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat3,
                           mu.start = mean(dat3$y), sigma.start = sd(dat3$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim3_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2","SN1", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim3_static_NO, k = 2),GAIC(sim3_static_NO2, k = 2), GAIC(sim3_static_GU, k = 2), GAIC(sim3_static_LO, k = 2), 
          GAIC(sim3_static_RG, k = 2), GAIC(sim3_static_exGAUS, k = 2), GAIC(sim3_static_NOF, k = 2), GAIC(sim3_static_PE, k = 2), 
          GAIC(sim3_static_PE2, k = 2),GAIC(sim3_static_SN1, k = 2),GAIC(sim3_static_SN2, k = 2),GAIC(sim3_static_TF, k = 2), GAIC(sim3_static_TF2, k = 2),
          GAIC(sim3_static_GT, k = 2),GAIC(sim3_static_JSU, k = 2), GAIC(sim3_static_JSUo, k = 2), GAIC(sim3_static_NET, k = 2),
          GAIC(sim3_static_SHASH, k = 2),GAIC(sim3_static_SHASHo, k = 2), GAIC(sim3_static_SHASHo2, k = 2),  
          GAIC(sim3_static_SEP2, k = 2),GAIC(sim3_static_SEP3, k = 2),GAIC(sim3_static_SEP4, k = 2),GAIC(sim3_static_ST1, k = 2),
          GAIC(sim3_static_ST2, k = 2),GAIC(sim3_static_ST3, k = 2), GAIC(sim3_static_ST4, k = 2), GAIC(sim3_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3", "3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym","both" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both" ,"both","both","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim3_pAIC_1a <- ggplot(sim3_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso,-ve skew,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

#compare
sim2_pAIC_1a + sim3_pAIC_1a





## SIM 3a: -ve to +ve skew, mesokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015  
beta_nu <- 0.0005

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)
nu_t <- 0.3   + beta_nu * t_index

# # observation randomly taken from normal distribution 
y3a <- rSST(t, mu = mu_t, sigma = sigma_t, nu = nu_t, tau = 1000)
y3a <- trend + seasonal_effect + y3a

# create dataset
dat3a <- data.frame(
  t = t_index,
  month = factor(month),
  y = y3a,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = nu_t,
  tau_true = 1000,
  dist = "SST"
)

# plot them to see
ggplot(dat3a, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat3a, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim3a_p1 <- ggplot(dat3a, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim3a_static_NO <- gamlss(y ~ 1, family = NO(), data = dat3a) # two-parameter
# NO2
sim3a_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat3a)
# GU
sim3a_static_GU <- gamlss(y ~ 1, family = GU(), data = dat3a)
# LO
sim3a_static_LO <- gamlss(y ~ 1, family = LO(), data = dat3a)
# RG
sim3a_static_RG <- gamlss(y ~ 1, family = RG(), data = dat3a)
# exGAUS
sim3a_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat3a, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim3a_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat3a)
# PE
sim3a_static_PE <- gamlss(y ~ 1, family = PE(), data = dat3a,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim3a_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat3a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim3a_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat3a) 
# SN2
sim3a_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat3a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim3a_static_TF <- gamlss(y ~ 1, family = TF(), data = dat3a)
# TF2
sim3a_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat3a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim3a_static_GT <- gamlss(y ~ 1, family = GT(), data = dat3a, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim3a_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat3a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim3a_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat3a,
                           mu.start = mean(dat3a$y), sigma.start = sd(dat3a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim3a_static_NET <- gamlss(y ~ 1, family = NET(), data = dat3a)
# SHASH
sim3a_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat3a,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim3a_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat3a,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim3a_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat3a,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim3a_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat3a, 
                           u.start = mean(dat3a$y), sigma.start = sd(dat3a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim3a_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat3,
                           mu.start = mean(dat3a$y), sigma.start = sd(dat3a$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim3a_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat3a,  # convergence issue
                           mu.start = mean(dat3a$y), sigma.start = sd(dat3a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim3a_static_SST <- gamlss(y ~ 1, family = SST(), data = dat3a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim3a_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat3a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim3a_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat3a,
                          mu.start = mean(dat3a$y), sigma.start = sd(dat3a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim3a_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat3a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim3a_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat3a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim3a_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat3a,
                          mu.start = mean(dat3a$y), sigma.start = sd(dat3a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim3a_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2","SN1", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim3a_static_NO, k = 2),GAIC(sim3a_static_NO2, k = 2), GAIC(sim3a_static_GU, k = 2), GAIC(sim3a_static_LO, k = 2), 
          GAIC(sim3a_static_RG, k = 2), GAIC(sim3a_static_exGAUS, k = 2), GAIC(sim3a_static_NOF, k = 2), GAIC(sim3a_static_PE, k = 2), 
          GAIC(sim3a_static_PE2, k = 2),GAIC(sim3a_static_SN1, k = 2),GAIC(sim3a_static_SN2, k = 2),GAIC(sim3a_static_TF, k = 2), GAIC(sim3a_static_TF2, k = 2),
          GAIC(sim3a_static_GT, k = 2),GAIC(sim3a_static_JSU, k = 2), GAIC(sim3a_static_JSUo, k = 2), GAIC(sim3a_static_NET, k = 2),
          GAIC(sim3a_static_SHASH, k = 2),GAIC(sim3a_static_SHASHo, k = 2), GAIC(sim3a_static_SHASHo2, k = 2),  
          GAIC(sim3a_static_SEP2, k = 2),GAIC(sim3a_static_SEP3, k = 2),GAIC(sim3a_static_SST, k = 2),GAIC(sim3a_static_ST1, k = 2),
          GAIC(sim3a_static_ST2, k = 2),GAIC(sim3a_static_ST3, k = 2), GAIC(sim3a_static_ST4, k = 2), GAIC(sim3a_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3", "3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym","both" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both" ,"both","lepto","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim3a_pAIC_1a <- ggplot(sim3a_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso,-ve to +ve skew,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# compare plots
sim2_pAIC_1a + sim3_pAIC_1a + sim3a_pAIC_1a




## SIM 4: symmetrical, leptokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015  


# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)


# # observation randomly taken from normal distribution 
y4 <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 1, tau = 4)
y4 <- trend + seasonal_effect + y4

# create dataset
dat4 <- data.frame(
  t = t_index,
  month = factor(month),
  y = y4,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 1,
  tau_true = 4,
  dist = "SST"
)

# plot them to see
ggplot(dat4, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat4, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim4_p1 <- ggplot(dat4, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim4_static_NO <- gamlss(y ~ 1, family = NO(), data = dat4) # two-parameter
# NO2
sim4_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat4)
# GU
sim4_static_GU <- gamlss(y ~ 1, family = GU(), data = dat4) # error
# LO
sim4_static_LO <- gamlss(y ~ 1, family = LO(), data = dat4)
# RG
sim4_static_RG <- gamlss(y ~ 1, family = RG(), data = dat4)
# exGAUS
sim4_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat4, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim4_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat4)
# PE
sim4_static_PE <- gamlss(y ~ 1, family = PE(), data = dat4,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim4_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat4, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim4_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat4) 
# SN2
sim4_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat4,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim4_static_TF <- gamlss(y ~ 1, family = TF(), data = dat4)
# TF2
sim4_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat4,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim4_static_GT <- gamlss(y ~ 1, family = GT(), data = dat4, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim4_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat4, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim4_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat4,
                            mu.start = mean(dat4$y), sigma.start = sd(dat4$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim4_static_NET <- gamlss(y ~ 1, family = NET(), data = dat4)
# SHASH
sim4_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat4,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim4_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat4,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim4_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat4,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim4_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat4, 
                            u.start = mean(dat4$y), sigma.start = sd(dat4$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim4_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat4,
                            mu.start = mean(dat4$y), sigma.start = sd(dat4$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim4_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat4,  
                            mu.start = mean(dat4$y), sigma.start = sd(dat4$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim4_static_SST <- gamlss(y ~ 1, family = SST(), data = dat4, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim4_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat4,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim4_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat4,
                           mu.start = mean(dat4$y), sigma.start = sd(dat4$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim4_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat4,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim4_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat4,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim4_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat4,
                           mu.start = mean(dat4$y), sigma.start = sd(dat4$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim4_static_AIC <- data.frame(
  models = c("NO", "NO2", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2","SN1", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim4_static_NO, k = 2),GAIC(sim4_static_NO2, k = 2), GAIC(sim4_static_LO, k = 2), 
          GAIC(sim4_static_RG, k = 2), GAIC(sim4_static_exGAUS, k = 2), GAIC(sim4_static_NOF, k = 2), GAIC(sim4_static_PE, k = 2), 
          GAIC(sim4_static_PE2, k = 2),GAIC(sim4_static_SN1, k = 2),GAIC(sim4_static_SN2, k = 2),GAIC(sim4_static_TF, k = 2), GAIC(sim4_static_TF2, k = 2),
          GAIC(sim4_static_GT, k = 2),GAIC(sim4_static_JSU, k = 2), GAIC(sim4_static_JSUo, k = 2), GAIC(sim4_static_NET, k = 2),
          GAIC(sim4_static_SHASH, k = 2),GAIC(sim4_static_SHASHo, k = 2), GAIC(sim4_static_SHASHo2, k = 2),  
          GAIC(sim4_static_SEP2, k = 2),GAIC(sim4_static_SEP3, k = 2),GAIC(sim4_static_SEP4, k = 2),GAIC(sim4_static_SST, k = 2),GAIC(sim4_static_ST1, k = 2),
          GAIC(sim4_static_ST2, k = 2),GAIC(sim4_static_ST3, k = 2), GAIC(sim4_static_ST4, k = 2), GAIC(sim4_static_ST5, k = 2)),
  param = c("2", "2", "2", "2",  
            "3","3","3","3","3","3", "3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4","4"),
  skewness = c("sym", "sym", "sym", "+ve", 
               "+ve","sym", "sym", "sym","both" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both" ,"both","both","lepto","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim4_pAIC_1a <- ggplot(sim4_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "lepto,sym,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



## SIM 4a: symmetrical, meso to leptokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015  
beta_tau   <- -0.0025    

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)
tau_t   <- 18   + beta_tau * t_index

# # observation randomly taken from normal distribution 
y4a <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 1, tau = tau_t)
y4a <- trend + seasonal_effect + y4a

# create dataset
dat4a <- data.frame(
  t = t_index,
  month = factor(month),
  y = y4a,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 1,
  tau_true = tau_t,
  dist = "SST"
)

# plot them to see
ggplot(dat4a, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat4a, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim4a_p1 <- ggplot(dat4a, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,40)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim4a_static_NO <- gamlss(y ~ 1, family = NO(), data = dat4a) # two-parameter
# NO2
sim4a_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat4a)
# GU
sim4a_static_GU <- gamlss(y ~ 1, family = GU(), data = dat4a) # error
# LO
sim4a_static_LO <- gamlss(y ~ 1, family = LO(), data = dat4a)
# RG
sim4a_static_RG <- gamlss(y ~ 1, family = RG(), data = dat4a)
# exGAUS
sim4a_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat4a, 
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim4a_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat4a)
# PE
sim4a_static_PE <- gamlss(y ~ 1, family = PE(), data = dat4a,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim4a_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat4a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim4a_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat4a) 
# SN2
sim4a_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat4a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim4a_static_TF <- gamlss(y ~ 1, family = TF(), data = dat4a)
# TF2
sim4a_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat4a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim4a_static_GT <- gamlss(y ~ 1, family = GT(), data = dat4a, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim4a_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat4a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim4a_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat4a,
                           mu.start = mean(dat4a$y), sigma.start = sd(dat4a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim4a_static_NET <- gamlss(y ~ 1, family = NET(), data = dat4a)
# SHASH
sim4a_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat4a,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim4a_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat4a,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim4a_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat4a,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim4a_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat4a, 
                           u.start = mean(dat4a$y), sigma.start = sd(dat4a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim4a_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat4a,
                           mu.start = mean(dat4a$y), sigma.start = sd(dat4a$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim4a_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat4a,  
                           mu.start = mean(dat4a$y), sigma.start = sd(dat4a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim4a_static_SST <- gamlss(y ~ 1, family = SST(), data = dat4a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim4a_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat4a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim4a_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat4a,
                          mu.start = mean(dat4a$y), sigma.start = sd(dat4a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim4a_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat4a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim4a_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat4a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim4a_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat4a,
                          mu.start = mean(dat4a$y), sigma.start = sd(dat4a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim4a_static_AIC <- data.frame(
  models = c("NO", "NO2","GU" ,"LO", "RG",
             "exGAUS","NOF", "PE", "PE2","SN1", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim4a_static_NO, k = 2),GAIC(sim4a_static_NO2, k = 2), GAIC(sim4a_static_GU, k = 2),GAIC(sim4a_static_LO, k = 2), 
          GAIC(sim4a_static_RG, k = 2), GAIC(sim4a_static_exGAUS, k = 2), GAIC(sim4a_static_NOF, k = 2), GAIC(sim4a_static_PE, k = 2), 
          GAIC(sim4a_static_PE2, k = 2),GAIC(sim4a_static_SN1, k = 2),GAIC(sim4a_static_SN2, k = 2),GAIC(sim4a_static_TF, k = 2), GAIC(sim4a_static_TF2, k = 2),
          GAIC(sim4a_static_GT, k = 2),GAIC(sim4a_static_JSU, k = 2), GAIC(sim4a_static_JSUo, k = 2), GAIC(sim4a_static_NET, k = 2),
          GAIC(sim4a_static_SHASH, k = 2),GAIC(sim4a_static_SHASHo, k = 2), GAIC(sim4a_static_SHASHo2, k = 2),  
          GAIC(sim4a_static_SEP2, k = 2),GAIC(sim4a_static_SEP3, k = 2),GAIC(sim4a_static_SEP4, k = 2),GAIC(sim4a_static_SST, k = 2),GAIC(sim4a_static_ST1, k = 2),
          GAIC(sim4a_static_ST2, k = 2),GAIC(sim4a_static_ST3, k = 2), GAIC(sim4a_static_ST4, k = 2), GAIC(sim4a_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3", "3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4","4"),
  skewness = c("sym", "sym", "-ve","sym", "+ve", 
               "+ve","sym", "sym", "sym","both" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto","lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both" ,"both","both","lepto","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim4a_pAIC_1a <- ggplot(sim4a_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso to lepto,sym,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


# compare plots
sim1_pAIC_1a + sim4_pAIC_1a + sim4a_pAIC_1a




## SIM 5: +ve skew, leptokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015  

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)

# # observation randomly taken from normal distribution 
y5 <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 3, tau = 4)
y5 <- trend + seasonal_effect + y5

# create dataset
dat5 <- data.frame(
  t = t_index,
  month = factor(month),
  y = y5,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 3,
  tau_true = 4,
  dist = "SST"
)

# plot them to see
ggplot(dat5, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat5, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim5_p1 <- ggplot(dat5, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim5_static_NO <- gamlss(y ~ 1, family = NO(), data = dat5) # two-parameter
# NO2
sim5_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat5)
# GU
sim5_static_GU <- gamlss(y ~ 1, family = GU(), data = dat5) # error
# LO
sim5_static_LO <- gamlss(y ~ 1, family = LO(), data = dat5)
# RG
sim5_static_RG <- gamlss(y ~ 1, family = RG(), data = dat5)
# exGAUS
sim5_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat5, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim5_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat5)
# PE
sim5_static_PE <- gamlss(y ~ 1, family = PE(), data = dat5,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim5_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat5, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim5_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat5) # not working
# SN2
sim5_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat5,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim5_static_TF <- gamlss(y ~ 1, family = TF(), data = dat5)
# TF2
sim5_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat5,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim5_static_GT <- gamlss(y ~ 1, family = GT(), data = dat5, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim5_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat5, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim5_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat5,
                            mu.start = mean(dat5$y), sigma.start = sd(dat5$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim5_static_NET <- gamlss(y ~ 1, family = NET(), data = dat5)
# SHASH
sim5_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat5,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim5_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat5,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim5_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat5,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim5_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat5, 
                            u.start = mean(dat5$y), sigma.start = sd(dat5$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim5_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat5,
                            mu.start = mean(dat5$y), sigma.start = sd(dat5$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim5_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat5,  # convergence issue
                            mu.start = mean(dat5$y), sigma.start = sd(dat5$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim5_static_SST <- gamlss(y ~ 1, family = SST(), data = dat5, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim5_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat5,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim5_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat5,
                           mu.start = mean(dat5$y), sigma.start = sd(dat5$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim5_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat5,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim5_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat5,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim5_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat5,
                           mu.start = mean(dat5$y), sigma.start = sd(dat5$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim5_static_AIC <- data.frame(
  models = c("NO", "NO2","GU" ,"LO", "RG",
             "exGAUS","NOF", "PE", "PE2", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim4a_static_NO, k = 2),GAIC(sim4a_static_NO2, k = 2), GAIC(sim4a_static_GU, k = 2),GAIC(sim4a_static_LO, k = 2), 
          GAIC(sim4a_static_RG, k = 2), GAIC(sim4a_static_exGAUS, k = 2), GAIC(sim4a_static_NOF, k = 2), GAIC(sim4a_static_PE, k = 2), 
          GAIC(sim4a_static_PE2, k = 2),GAIC(sim4a_static_SN2, k = 2),GAIC(sim4a_static_TF, k = 2), GAIC(sim4a_static_TF2, k = 2),
          GAIC(sim4a_static_GT, k = 2),GAIC(sim4a_static_JSU, k = 2), GAIC(sim4a_static_JSUo, k = 2), GAIC(sim4a_static_NET, k = 2),
          GAIC(sim4a_static_SHASH, k = 2),GAIC(sim4a_static_SHASHo, k = 2), GAIC(sim4a_static_SHASHo2, k = 2),  
          GAIC(sim4a_static_SEP2, k = 2),GAIC(sim4a_static_SEP3, k = 2),GAIC(sim4a_static_SST, k = 2),GAIC(sim4a_static_ST1, k = 2),
          GAIC(sim4a_static_ST2, k = 2),GAIC(sim4a_static_ST3, k = 2), GAIC(sim4a_static_ST4, k = 2), GAIC(sim4a_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve","sym", "+ve", 
               "+ve","sym", "sym", "sym" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto","lepto", "lepto",
               "lepto","meso", "both", "both", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both" ,"both","lepto","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim5_pAIC_1a <- ggplot(sim5_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "lepto,+ve skew,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )



## SIM 5a: -ve skew, leptokurtotic, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015  

# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)

# # observation randomly taken from normal distribution 
y5a <- rSST(t, mu = mu_t, sigma = sigma_t, nu = 0.3, tau = 4)
y5a <- trend + seasonal_effect + y5a

# create dataset
dat5a <- data.frame(
  t = t_index,
  month = factor(month),
  y = y5a,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = 0.3,
  tau_true = 4,
  dist = "SST"
)

# plot them to see
ggplot(dat5a, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat5a, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim5a_p1 <- ggplot(dat5a, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim5a_static_NO <- gamlss(y ~ 1, family = NO(), data = dat5a) # two-parameter
# NO2
sim5a_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat5a)
# GU
sim5a_static_GU <- gamlss(y ~ 1, family = GU(), data = dat5a) 
# LO
sim5a_static_LO <- gamlss(y ~ 1, family = LO(), data = dat5a)
# RG
sim5a_static_RG <- gamlss(y ~ 1, family = RG(), data = dat5a) # error sigma
# exGAUS
sim5a_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat5a, # error sigma
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim5a_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat5a)
# PE
sim5a_static_PE <- gamlss(y ~ 1, family = PE(), data = dat5a,
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim5a_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat5a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim5a_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat5a) # not working
# SN2
sim5a_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat5a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim5a_static_TF <- gamlss(y ~ 1, family = TF(), data = dat5a)
# TF2
sim5a_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat5a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim5a_static_GT <- gamlss(y ~ 1, family = GT(), data = dat5a, 
                         method = mixed(10,200),
                         control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim5a_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat5a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim5a_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat5a,
                           mu.start = mean(dat5a$y), sigma.start = sd(dat5a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim5a_static_NET <- gamlss(y ~ 1, family = NET(), data = dat5a)
# SHASH
sim5a_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat5a,
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim5a_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat5a,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim5a_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat5a,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim5a_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat5a, 
                           u.start = mean(dat5a$y), sigma.start = sd(dat5a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim5a_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat5a,
                           mu.start = mean(dat5a$y), sigma.start = sd(dat5a$y), 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim5a_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat5a,  
                           mu.start = mean(dat5a$y), sigma.start = sd(dat5a$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim5a_static_SST <- gamlss(y ~ 1, family = SST(), data = dat5a, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim5a_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat5a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim5a_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat5a,
                          mu.start = mean(dat5a$y), sigma.start = sd(dat5a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim5a_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat5a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim5a_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat5a,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim5a_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat5a,
                          mu.start = mean(dat5a$y), sigma.start = sd(dat5a$y),
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim5a_static_AIC <- data.frame(
  models = c("NO", "NO2","GU" ,"LO",
             "NOF", "PE", "PE2", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SEP4","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim4a_static_NO, k = 2),GAIC(sim4a_static_NO2, k = 2), GAIC(sim4a_static_GU, k = 2),GAIC(sim4a_static_LO, k = 2), 
          GAIC(sim4a_static_NOF, k = 2), GAIC(sim4a_static_PE, k = 2), 
          GAIC(sim4a_static_PE2, k = 2),GAIC(sim4a_static_SN2, k = 2),GAIC(sim4a_static_TF, k = 2), GAIC(sim4a_static_TF2, k = 2),
          GAIC(sim4a_static_GT, k = 2),GAIC(sim4a_static_JSU, k = 2), GAIC(sim4a_static_JSUo, k = 2), GAIC(sim4a_static_NET, k = 2),
          GAIC(sim4a_static_SHASH, k = 2),GAIC(sim4a_static_SHASHo, k = 2), GAIC(sim4a_static_SHASHo2, k = 2),  
          GAIC(sim4a_static_SEP2, k = 2),GAIC(sim4a_static_SEP3, k = 2),GAIC(sim4a_static_SEP4, k = 2),GAIC(sim4a_static_SST, k = 2),GAIC(sim4a_static_ST1, k = 2),
          GAIC(sim4a_static_ST2, k = 2),GAIC(sim4a_static_ST3, k = 2), GAIC(sim4a_static_ST4, k = 2), GAIC(sim4a_static_ST5, k = 2)),
  param = c("2", "2", "2", "2",
            "3","3","3","3","3","3", 
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4", "4"),
  skewness = c("sym", "sym", "-ve","sym",
               "sym", "sym", "sym" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto","lepto",
               "meso", "both", "both", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both","both","both","lepto","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim5a_pAIC_1a <- ggplot(sim5a_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "lepto,+ve skew,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


## SIM 5b: -ve to +ve skew, meso to lepto, var(time) --> SST() ===============================================================

set.seed(123)
t <- 6000
t_index <- seq_len(t)
n <- 360 # 25 years
month <- rep(1:12, times = t/12)
seasonal_effect <- sin(2 * pi * month / 12) * 5 # Example sinusoidal seasonality
trend <- 0.0025 * (1:t)


beta_mu <- 0.001
beta_sigma <- 0.0015  
beta_nu <- 0.0005
beta_tau   <- -0.0025    


# changes with time
mu_t <- 10   + beta_mu * t_index
sigma_t <- pmax(0, 2 + beta_sigma * t_index)
nu_t <- 0.3   + beta_nu * t_index
tau_t   <- 18   + beta_tau * t_index

# # observation randomly taken from normal distribution 
y5b <- rSST(t, mu = mu_t, sigma = sigma_t, nu = nu_t, tau = tau_t)
y5b <- trend + seasonal_effect + y5b

# create dataset
dat5b <- data.frame(
  t = t_index,
  month = factor(month),
  y = y5b,
  sd_true = sigma_t,
  mu_true = mu_t,
  nu_true = nu_t,
  tau_true = tau_t,
  dist = "SST"
)

# plot them to see
ggplot(dat5b, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat5b, aes(x = y)) +
  geom_density()

# histogram and non-parametric density estimate
sim5b_p1 <- ggplot(dat5b, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 2) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,150)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 


### fit models of different family distribution ###

# NO
sim5b_static_NO <- gamlss(y ~ 1, family = NO(), data = dat5b) # two-parameter
# NO2
sim5b_static_NO2 <- gamlss(y ~ 1, family = NO2(), data = dat5b)
# GU
sim5b_static_GU <- gamlss(y ~ 1, family = GU(), data = dat5b)
# LO
sim5b_static_LO <- gamlss(y ~ 1, family = LO(), data = dat5b)
# RG
sim5b_static_RG <- gamlss(y ~ 1, family = RG(), data = dat5b)
# exGAUS
sim5b_static_exGAUS <- gamlss(y ~ 1, family = exGAUS(), data = dat5b, 
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # three-parameter
# NOF
sim5b_static_NOF <- gamlss(y ~ 1, family = NOF(), data = dat5b)
# PE
sim5b_static_PE <- gamlss(y ~ 1, family = PE(), data = dat5b,
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# PE2
sim5b_static_PE2 <- gamlss(y ~ 1, family = PE2(), data = dat5b, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SN1
sim5b_static_SN1 <- gamlss(y ~ 1, family = SN1(), data = dat5b) # not working
# SN2
sim5b_static_SN2 <- gamlss(y ~ 1, family = SN2(), data = dat5b,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# TF
sim5b_static_TF <- gamlss(y ~ 1, family = TF(), data = dat5b)
# TF2
sim5b_static_TF2 <- gamlss(y ~ 1, family = TF2(), data = dat5b,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# GT
sim5b_static_GT <- gamlss(y ~ 1, family = GT(), data = dat5b, 
                          method = mixed(10,200),
                          control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) # four-parameter
# JSU
sim5b_static_JSU <- gamlss(y ~ 1, family = JSU(), data = dat5b, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# JSUo
sim5b_static_JSUo <- gamlss(y ~ 1, family = JSUo(), data = dat5b,
                            mu.start = mean(dat5b$y), sigma.start = sd(dat5b$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# NET
sim5b_static_NET <- gamlss(y ~ 1, family = NET(), data = dat5b)
# SHASH
sim5b_static_SHASH <- gamlss(y ~ 1, family = SHASH(), data = dat5b,
                             method = mixed(10,200),
                             control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo
sim5b_static_SHASHo <- gamlss(y ~ 1, family = SHASHo(), data = dat5b,
                              method = mixed(10,200),
                              control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SHASHo2
sim5b_static_SHASHo2 <- gamlss(y ~ 1, family = SHASHo2(), data = dat5b,
                               method = mixed(10,200),
                               control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP2
sim5b_static_SEP2 <- gamlss(y ~ 1, family = SEP2(), data = dat5b, 
                            u.start = mean(dat5b$y), sigma.start = sd(dat5b$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# SEP3
sim5b_static_SEP3 <- gamlss(y ~ 1, family = SEP3(), data = dat5b,
                            mu.start = mean(dat5b$y), sigma.start = sd(dat5b$y), 
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SEP4
sim5b_static_SEP4 <- gamlss(y ~ 1, family = SEP4(), data = dat5b,  # convergence issue
                            mu.start = mean(dat5b$y), sigma.start = sd(dat5b$y),
                            method = mixed(10,200),
                            control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE)) 
# SST
sim5b_static_SST <- gamlss(y ~ 1, family = SST(), data = dat5b, 
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))  
# ST1
sim5b_static_ST1 <- gamlss(y ~ 1, family = ST1(), data = dat5b,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST2
sim5b_static_ST2 <- gamlss(y ~ 1, family = ST2(), data = dat5b,
                           mu.start = mean(dat5b$y), sigma.start = sd(dat5b$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST3
sim5b_static_ST3 <- gamlss(y ~ 1, family = ST3(), data = dat5b,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST4
sim5b_static_ST4 <- gamlss(y ~ 1, family = ST4(), data = dat5b,
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))
# ST5
sim5b_static_ST5 <- gamlss(y ~ 1, family = ST5(), data = dat5b,
                           mu.start = mean(dat5b$y), sigma.start = sd(dat5b$y),
                           method = mixed(10,200),
                           control = gamlss.control(n.cyc = 200, c.crit = 0.01, trace = FALSE))



# scale the AIC 
sim5b_static_AIC <- data.frame(
  models = c("NO", "NO2", "GU", "LO", "RG",
             "exGAUS","NOF", "PE", "PE2","SN1", "SN2", "TF", "TF2", 
             "GT", "JSU", "JSUo", "NET", "SHASH", "SHASHo", "SHASHo2",
             "SEP2","SEP3","SST", "ST1", "ST2", "ST3", "ST4", "ST5"),
  AIC = c(GAIC(sim5b_static_NO, k = 2),GAIC(sim5b_static_NO2, k = 2), GAIC(sim5b_static_GU, k = 2), GAIC(sim5b_static_LO, k = 2), 
          GAIC(sim5b_static_RG, k = 2), GAIC(sim5b_static_exGAUS, k = 2), GAIC(sim5b_static_NOF, k = 2), GAIC(sim5b_static_PE, k = 2), 
          GAIC(sim5b_static_PE2, k = 2),GAIC(sim5b_static_SN1, k = 2),GAIC(sim5b_static_SN2, k = 2),GAIC(sim5b_static_TF, k = 2), GAIC(sim5b_static_TF2, k = 2),
          GAIC(sim5b_static_GT, k = 2),GAIC(sim5b_static_JSU, k = 2), GAIC(sim5b_static_JSUo, k = 2), GAIC(sim5b_static_NET, k = 2),
          GAIC(sim5b_static_SHASH, k = 2),GAIC(sim5b_static_SHASHo, k = 2), GAIC(sim5b_static_SHASHo2, k = 2),  
          GAIC(sim5b_static_SEP2, k = 2),GAIC(sim5b_static_SEP3, k = 2),GAIC(sim5b_static_SST, k = 2),GAIC(sim5b_static_ST1, k = 2),
          GAIC(sim5b_static_ST2, k = 2),GAIC(sim5b_static_ST3, k = 2), GAIC(sim5b_static_ST4, k = 2), GAIC(sim5b_static_ST5, k = 2)),
  param = c("2", "2", "2", "2", "2", 
            "3","3","3","3","3","3", "3", "3",
            "4","4","4","4","4","4","4","4","4","4","4","4","4", "4", "4"),
  skewness = c("sym", "sym", "-ve", "sym", "+ve", 
               "+ve","sym", "sym", "sym","both" ,"both","sym", "sym",
               "sym", "both", "both", "sym", "both", "both", "both",
               "both","both","both","both","both", "both", "both", "both"),
  kurtosis = c("meso", "meso", "lepto", "lepto", "lepto",
               "lepto","meso", "both", "both", "lepto", "lepto","lepto", "lepto",
               "both", "lepto", "lepto", "lepto", "both", "both", "both",
               "both" ,"both","lepto","lepto", "lepto", "lepto", "lepto", "lepto")
)

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim5b_pAIC_1a <- ggplot(sim5b_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "meso to lepto,-ve to +ve skew,var(time); SST()",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )





##### SIM 6: simple, sym, meso --> SST() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
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

# plot them to see
ggplot(dat6, aes(x = t)) +
  geom_line(aes(y = y), color = "darkred", linewidth = 0.8) +
  labs(x = "Time", y = "Value") +
  theme_minimal(base_size = 14) 


# look at density plot
ggplot(dat6, aes(x = y)) +
  geom_density() +
  geom_vline(xintercept = mean(dat6$y), linetype = "dashed", col = "grey") +
  theme_minimal()

# histogram and non-parametric density estimate
sim6_p1 <- ggplot(dat6, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  xlab("value") + 
  ylab("Density") +
  xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

mean(dat6$y)
median(dat6$y)
sd(dat6$y)
skewness(dat6$y)
kurtosis(dat6$y)

### fit models of different family distribution ###

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



# scale the AIC 
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

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim6_pAIC_1a <- ggplot(sim6_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 6: meso, sym (shapes)",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

summary(sim6_static_ST3)
summary(sim6_static_SN2)
summary(sim6_static_TF)

param_color <- c("4" = "blue4", "3k" = "black", "4k" = "black", "3s" = "gray48", "2" = "grey88",  "3" = "grey88")

sim6_pAIC <- ggplot(sim6_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = param_color, name = "param") +
  #scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 6: meso, sym",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


sim6_pAIC_1a + sim6_pAIC

# param recovery 
param_table <- data.frame(
  models = c("true","NO", "SN2", "PE", "exGAUS", "TF", "SEP3"),
  param = c("-", "2", "3", "3", "3", "3", "4"),
  AIC = c("-",round(AIC(sim6_static_NO),2), round(AIC(sim6_static_SN2),2), round(AIC(sim6_static_PE),2),
          round(AIC(sim6_static_exGAUS),2), round(AIC(sim6_static_TF),2), round(AIC(sim6_static_SEP3),2)),
  deviance = c("-",round(deviance(sim6_static_NO),2), round(deviance(sim6_static_SN2),2), round(deviance(sim6_static_PE),2),
               round(deviance(sim6_static_exGAUS),2), round(deviance(sim6_static_TF),2), round(deviance(sim6_static_SEP3),2)),
  mu = c("10",round(sim6_static_NO$mu.coefficients,3), round(sim6_static_SN2$mu.coefficients,3), round(sim6_static_PE$mu.coefficients,3),
         round(sim6_static_exGAUS$mu.coefficients,3), round(sim6_static_TF$mu.coefficients,3), round(sim6_static_SEP3$mu.coefficients,3)),
  sigma = c("2",round(exp(sim6_static_NO$sigma.coefficients),3), round(exp(sim6_static_SN2$sigma.coefficients),3), round(exp(sim6_static_PE$sigma.coefficients),3),
            round(exp(sim6_static_exGAUS$sigma.coefficients),3), round(exp(sim6_static_TF$sigma.coefficients),3), round(exp(sim6_static_SEP3$sigma.coefficients),3)),
  nu = c("1","-",round(exp(sim6_static_SN2$nu.coefficients),3), "-",
         "-", "-", round(exp(sim6_static_SEP3$nu.coefficients),3)),
  tau = c("1000","-", "-", round(exp(sim6_static_PE$nu.coefficients),3), 
          "-", round(exp(sim6_static_TF$nu.coefficients),3), round(exp(sim6_static_SEP3$tau.coefficients),3))
)

library(knitr)
param_table %>%
  kable(caption = "Parameters recovery for sim6")





##### SIM 7: simple, +ve skew, meso --> SST() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
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

# histogram and non-parametric density estimate
sim7_p1 <- ggplot(dat7, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  xlab("value") + 
  ylab("Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

mean(dat6$y)
median(dat6$y)
sd(dat6$y)
skewness(dat6$y)
kurtosis(dat6$y)

### fit models of different family distribution ###

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



# scale the AIC 
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
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim7_pAIC_1a <- ggplot(sim7_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim7: meso, +ve skew (shape)",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


param_color <- c("4" = "blue4", "3k" = "black", "3s" = "gray48",
                 "2" = "grey88",  "3" = "grey88", "4k" = "black")

sim7_pAIC <- ggplot(sim7_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = param_color, name = "param") +
  #scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 7: meso, +ve skew",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# combineation 
sim7_pAIC_1a + sim7_pAIC


# param recovery 
param_table <- data.frame(
  models = c("true","SN2", "SEP3", "exGAUS", "TF", "NO2"),
  param = c("-", "3", "4", "3", "3", "2"),
  AIC = c("-",round(AIC(sim7_static_SN2),2), round(AIC(sim7_static_SEP3),2), round(AIC(sim7_static_exGAUS),2),
          round(AIC(sim7_static_TF),2), round(AIC(sim7_static_NO2),2)),
  deviance = c("-",round(deviance(sim7_static_SN2),2), round(deviance(sim7_static_SEP3),2), round(deviance(sim7_static_exGAUS),2),
          round(deviance(sim7_static_TF),2), round(deviance(sim7_static_NO2),2)),
  mu = c("10",round(sim7_static_SN2$mu.coefficients,3), round(sim7_static_SEP3$mu.coefficients,3), round(sim7_static_exGAUS$mu.coefficients,3),
         round(sim7_static_TF$mu.coefficients,3), round(sim7_static_NO2$mu.coefficients,3)),
  sigma = c("2",round(sim7_static_SN2$sigma.coefficients,3), round(sim7_static_SEP3$sigma.coefficients,3), round(sim7_static_exGAUS$mu.coefficients,3),
            round(sim7_static_TF$sigma.coefficients,3), round(sim7_static_NO2$sigma.coefficients,3)),
  ##nu = c("3",round(sim7_static_SN2$sigma.coefficients,3), round(sim7_static_SEP3$sigma.coefficients,3), round(sim7_static_exGAUS$mu.coefficients,3),
            round(sim7_static_TF$sigma.coefficients,3), round(sim7_static_NO2$sigma.coefficients,3)),
  tau = c("1000","-", "-", round(exp(sim6_static_PE$nu.coefficients),3), 
          "-", round(exp(sim6_static_TF$nu.coefficients),3), round(exp(sim6_static_SEP3$tau.coefficients),3))
)

param_table %>%
  kable(caption = "Predicted parameters for sim7")

# SEP3 says its leptokurtic slightly (< 2)



##### SIM 8: simple, sym, lepto --> SST() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
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

# histogram and non-parametric density estimate
sim8_p1 <- ggplot(dat8, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  xlab("value") + 
  ylab("Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

mean(dat6$y)
median(dat6$y)
sd(dat6$y)
skewness(dat8$y)
kurtosis(dat8$y)

### fit models of different family distribution ###

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



# scale the AIC 
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

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim8_pAIC_1a <- ggplot(sim8_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 8: lepto, sym (shape)",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

param_color <- c("4" = "blue4", "3k" = "black", "3s" = "gray48", "2" = "grey88",  "3" = "grey88")

sim8_pAIC <- ggplot(sim8_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = param_color, name = "param") +
  #scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 8: lepto, sym",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# combeanation 
sim8_pAIC_1a + sim8_pAIC


# param recovery 
param_table <- data.frame(
  models = c("true","SST", "ST3", "JSU", "TF"),
  mu = c("10",sim8_static_SST$mu.coefficients, sim8_static_ST3$mu.coefficients, sim8_static_JSU$mu.coefficients, sim8_static_TF$mu.coefficients),
  sigma = c("2",exp(sim8_static_SST$sigma.coefficients), exp(sim8_static_ST3$sigma.coefficients), exp(sim8_static_JSU$sigma.coefficients), exp(sim8_static_TF$sigma.coefficients)),
  nu = c("1", exp(sim8_static_SST$nu.coefficients), exp(sim8_static_ST3$nu.coefficients), sim8_static_JSU$nu.coefficients, exp(sim8_static_TF$nu.coefficients)),
  tau = c("4",exp(sim8_static_SST$tau.coefficients) + 2, exp(sim8_static_ST3$tau.coefficients), exp(sim8_static_JSU$tau.coefficients), "-")
)

param_table %>%
  kable(caption = "Predicted parameters for sim8")




##### SIM 9: simple, +ve, lepto --> SST() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
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

# histogram and non-parametric density estimate
sim9_p1 <- ggplot(dat9, aes(x = y)) +
  geom_histogram(aes(y = ..density..), colour = "black", fill = "white", binwidth = 0.5) + 
  xlab("value") + 
  ylab("Density") +
  #xlim(x = c(0,18)) +
  geom_density(colour = "black", fill = "grey", alpha = 0.5) +
  theme_bw() 

mean(dat6$y)
median(dat6$y)
sd(dat6$y)
skewness(dat6$y)
kurtosis(dat6$y)

### fit models of different family distribution ###

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



# scale the AIC 
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

# plot colored by skewness and shape by kurtosis 
skew_color <- c("sym" = "blue4", "both" = "black", "+ve" = "gray48", "-ve" = "grey88")

sim9_pAIC_1a <- ggplot(sim9_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = skewness, shape = kurtosis)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = skew_color, name = "param") +
  scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 9: lepto, +ve skew (shape)",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )


param_color <- c("4" = "blue4", "3k" = "black", "3s" = "gray48", "2" = "grey88",  "3" = "grey88")

sim9_pAIC <- ggplot(sim9_static_AIC, aes(x = AIC, y = reorder(models, AIC), color = param)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = param_color, name = "param") +
  #scale_shape_manual(values = c("lepto" = 19, "both" = 15, "meso" = 17)) +
  labs(title = "sim 9: lepto, +ve skew",x = "AIC", y = "value") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_line(color = "grey90"),
    legend.position = "none"
  )

# combeanation 
sim9_pAIC_1a + sim9_pAIC


# param recovery 
param_table <- data.frame(
  models = c("true","SST", "ST3", "JSU", "TF"),
  mu = c("10",sim8_static_SST$mu.coefficients, sim8_static_ST3$mu.coefficients, sim8_static_JSU$mu.coefficients, sim8_static_TF$mu.coefficients),
  sigma = c("2",exp(sim8_static_SST$sigma.coefficients), exp(sim8_static_ST3$sigma.coefficients), exp(sim8_static_JSU$sigma.coefficients), exp(sim8_static_TF$sigma.coefficients)),
  nu = c("1", exp(sim8_static_SST$nu.coefficients), exp(sim8_static_ST3$nu.coefficients), sim8_static_JSU$nu.coefficients, exp(sim8_static_TF$nu.coefficients)),
  tau = c("4",exp(sim8_static_SST$tau.coefficients) + 2, exp(sim8_static_ST3$tau.coefficients), exp(sim8_static_JSU$tau.coefficients), "-")
)

param_table %>%
  kable(caption = "Predicted parameters for sim8")


# compare
sim6_pAIC_1a + sim7_pAIC_1a + sim8_pAIC_1a + sim9_pAIC_1a





##### SIM 10: simple, +ve, lepto --> SST() ==================================================================

set.seed(123)
t2 <- 300
t_index2 <- seq_len(t2)

# constant
mu0 <- 10
sigma0 <- 2

# # observation randomly taken from SST
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
