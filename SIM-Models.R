# =============================================================================
# Simulation Study: model fitting and comparison
# "Analysis of Survival Studies with Cure Rates using Soft BART"
# Mandal, Sinha & Pal
#
# Methods compared (matching Table 1):
#   SBART-Cure  -- proposed (multiplicative Gamma frailty)        [hand-coded MCMC]
#   RSF         -- random survival forest (no cure, no clustering)[randomForestSRC]
#   WGF         -- Weibull PH + Gamma frailty (Sahu et al. 1997)  [parfm, MLE]
#   SGF         -- semiparametric Cox + Gamma frailty (Sargent 98)[frailtyEM, EM]
#   M-CIS       -- multivariate promotion-time cure (Chen 2002)   [hand-coded MCMC]
#
# Off-the-shelf competitors (WGF, SGF, RSF) are fit with standard packages and
# their conditional survival S(t | W_i, x_ij) is reconstructed from the fitted
# baseline, regression coefficients, and per-cluster frailty estimates W_hat_i.
# WGF and SGF are pure PH models with no cure fraction -> MSEcure = NA.
#
# IMPORTANT: source the data-generation script first. It must define:
#   alpha_true, theta_true, b_true(), beta_true, beta_stable, beta_add,
#   lambda_cens_1, lambda_cens_2, lambda_cens_3, and generate_data().
# =============================================================================

rm(list = ls())
set.seed(123)

library(survival)
library(randomForestSRC)
library(SoftBart)
library(msm)
library(parfm)
library(frailtyEM)
library(parallel)

# source("SIM-DATAgen.R")   # <-- run this first

# =============================================================================
# SECTION 1: TRUE survival and cure functions
# =============================================================================

# --- Multiplicative frailty (Sims 1, 2, 4) -----------------------------------
# S(t | W, x) = exp(-W * beta * integral_0^t Phi(b(u,x)) f_lambda(u) du)
true_survival_mult <- function(t_grid, W, x1, x2, beta, n_mc = 4000) {
  T_mc <- rgamma(n_mc, shape = alpha_true, rate = theta_true)
  phi  <- pnorm(b_true(T_mc, x1, x2))
  sapply(t_grid, function(tg) if (tg <= 0) 1.0 else exp(-W * beta * mean(phi * (T_mc <= tg))))
}
true_cure_mult <- function(W, x1, x2, beta, n_mc = 4000) {
  T_mc <- rgamma(n_mc, shape = alpha_true, rate = theta_true)
  exp(-W * beta * mean(pnorm(b_true(T_mc, x1, x2))))
}

# --- Additive frailty (Sim 3) ------------------------------------------------
# h(t | W, x) = (W + beta * Phi(b(t,x))) f_lambda(t)
# S(t | W, x) = exp(-W * F_lambda(t) - beta * J(t,x)),  J(t,x)=int_0^t Phi(b) f_lambda
# pi(x | W)   = exp(-W - beta * I(x)),  I(x) = J(inf, x)
true_survival_add <- function(t_grid, W, x1, x2, beta, n_mc = 4000) {
  T_mc <- rgamma(n_mc, shape = alpha_true, rate = theta_true)
  phi  <- pnorm(b_true(T_mc, x1, x2))
  sapply(t_grid, function(tg) {
    if (tg <= 0) return(1.0)
    F_t <- mean(T_mc <= tg)
    J_t <- mean(phi * (T_mc <= tg))
    exp(-W * F_t - beta * J_t)
  })
}
true_cure_add <- function(W, x1, x2, beta, n_mc = 4000) {
  T_mc <- rgamma(n_mc, shape = alpha_true, rate = theta_true)
  exp(-W - beta * mean(pnorm(b_true(T_mc, x1, x2))))
}

# =============================================================================
# SECTION 2: Metric machinery
#
# The truth is computed ONCE per data set (fixed MC draws) so every method is
# scored against the same gold standard. All methods share one evaluation grid:
#   t_grid = seq(0, 99th-pct of observed y, length = G), G = 100.
# =============================================================================

make_truth <- function(dat, beta, G = 100, sim_type = "multiplicative", n_mc = 4000) {
  t_grid <- seq(0, quantile(dat$y, 0.99), length.out = G)
  tsf <- if (sim_type == "additive") true_survival_add else true_survival_mult
  tcf <- if (sim_type == "additive") true_cure_add     else true_cure_mult
  n <- nrow(dat)
  Sg <- matrix(NA_real_, n, G); So <- numeric(n); pc <- numeric(n)
  for (s in 1:n) {
    S <- tsf(c(t_grid, dat$y[s]), dat$W_true[s], dat$x1[s], dat$x2[s], beta, n_mc)
    Sg[s, ] <- S[1:G]; So[s] <- S[G + 1]
    pc[s]   <- tcf(dat$W_true[s], dat$x1[s], dat$x2[s], beta, n_mc)
  }
  list(t_grid = t_grid, Sg = Sg, So = So, pc = pc)
}

# Generic metric for any method exposing S_hat_fn(t_vec, s) and (optionally) cure_hat_fn(s).
metric_model <- function(dat, truth, S_hat_fn, cure_hat_fn = NULL) {
  t_grid <- truth$t_grid; G <- length(t_grid); n <- nrow(dat)
  mc <- mo <- numeric(n); mk <- rep(NA_real_, n)
  for (s in 1:n) {
    Sh <- S_hat_fn(c(t_grid, dat$y[s]), s)
    mc[s] <- mean((Sh[1:G] - truth$Sg[s, ])^2)
    mo[s] <- (Sh[G + 1] - truth$So[s])^2
    if (!is.null(cure_hat_fn)) mk[s] <- (cure_hat_fn(s) - truth$pc[s])^2
  }
  c(MSEcurve = mean(mc), MSEobs = mean(mo),
    MSEcure  = if (all(is.na(mk))) NA_real_ else mean(mk, na.rm = TRUE))
}

# =============================================================================
# SECTION 3: SBART-Cure (proposed)  -- unchanged from working version
#   Two time scales: t_orig for ALL Gamma computations, t_scaled (0,1) for the
#   forest input only. Posterior survival/cure stored per post-burn iteration.
# =============================================================================

scale01_fixed <- function(x, x_min, x_max) (x - x_min) / (x_max - x_min)

run_sbart_mcmc <- function(dat, num_iter = 5000, num_burn = 2500,
                           shape_prior_b = 0.1, rate_prior_b = 0.1,
                           shape_prior_alpha = 0.1, rate_prior_alpha = 0.1,
                           shape_prior_theta = 0.1, rate_prior_theta = 0.1,
                           shape_prior_eta = 4, rate_prior_eta = 0.01) {

  cat("  Running SBART MCMC...\n")
  t_orig  <- dat$y;  status <- dat$delta;  cluster <- dat$cluster
  X_cov   <- cbind(dat$x1, dat$x2)
  n <- nrow(dat);  p <- ncol(X_cov);  N_clust <- length(unique(cluster))
  no_exact <- sum(status == 1)
  t_min <- min(t_orig);  t_max <- max(t_orig)
  t_scaled <- scale01_fixed(t_orig, t_min, t_max)

  beta <- 1.0;  alpha <- 1.0;  theta <- 1.0;  eta <- 2.0
  w <- rgamma(N_clust, shape = eta, rate = eta)

  hypers <- Hypers(cbind(t_scaled, X_cov), t_scaled, num_tree = 200)
  opts   <- Opts(cache_trees = FALSE, update_sigma = TRUE, update_s = FALSE,
                 update_alpha = FALSE, update_sigma_mu = FALSE, update_tau = TRUE,
                 num_burn = num_burn, num_thin = 1, num_save = num_iter - num_burn)
  my_forest <- MakeForest(hypers, opts)

  beta_save <- alpha_save <- theta_save <- eta_save <- numeric(num_iter)
  w_save <- matrix(NA, num_iter, N_clust)
  q_vec  <- integer(n);  m_vec <- integer(n)

  G_pred <- 100;  n_post <- num_iter - num_burn
  surv_post     <- array(NA, dim = c(n_post, n, G_pred))
  cure_post     <- matrix(NA, nrow = n_post, ncol = n)
  surv_obs_post <- matrix(NA, nrow = n_post, ncol = n)
  t_grid_pred   <- seq(0, quantile(t_orig, 0.99), length.out = G_pred)
  n_mc_pred     <- 1000

  slice_sample <- function(logf, x0, w_step = 1.0, m = 10, lower = 1e-6, upper = Inf) {
    logy <- logf(x0) - rexp(1)
    u <- runif(1, 0, w_step);  L <- x0 - u;  R <- x0 + (w_step - u)
    J <- floor(m * runif(1));  K <- (m - 1) - J
    while (J > 0 && L > lower && logf(L) > logy) { L <- L - w_step; J <- J - 1 }
    while (K > 0 && R < upper && logf(R) > logy) { R <- R + w_step; K <- K - 1 }
    L <- max(L, lower);  R <- min(R, upper)
    for (i in 1:100) {
      x1 <- runif(1, L, R)
      if (logf(x1) >= logy) return(x1)
      if (x1 < x0) L <- x1 else R <- x1
    }
    return(x0)
  }

  start_time <- Sys.time()
  for (iter in 1:num_iter) {
    z_all <- NULL;  X_store <- NULL;  time_points <- NULL;  G_all_orig <- NULL

    for (i in 1:n) {
      g_orig <- NULL;  g_scaled <- NULL;  w_i <- w[cluster[i]]
      pois_rate <- w_i * beta * pgamma(t_orig[i], shape = alpha, rate = theta)
      q_vec[i]  <- rpois(1, pois_rate)
      if (q_vec[i] == 0) {
        m_vec[i] <- 0
      } else {
        c_vals        <- runif(q_vec[i], 0, pois_rate)
        a_vals_orig   <- qgamma(c_vals / (w_i * beta), shape = alpha, rate = theta)
        a_vals_scaled <- scale01_fixed(a_vals_orig, t_min, t_max)
        a_X_mat <- cbind(a_vals_scaled,
                         matrix(rep(X_cov[i, ], length(a_vals_orig)),
                                nrow = length(a_vals_orig), ncol = p, byrow = TRUE))
        l_vals   <- my_forest$do_predict(a_X_mat)
        keep     <- runif(q_vec[i]) < (1 - pnorm(l_vals))
        g_orig   <- a_vals_orig[keep];  g_scaled <- a_vals_scaled[keep]
        m_vec[i] <- length(g_orig)
      }

      if (status[i] == 1) {
        z_t <- msm::rtnorm(1, mean = my_forest$do_predict(
                 matrix(c(t_scaled[i], X_cov[i, ]), nrow = 1)),
                 sd = 1, lower = 0, upper = Inf)
        if (m_vec[i] > 0) {
          z_g <- numeric(m_vec[i])
          for (j in 1:m_vec[i])
            z_g[j] <- msm::rtnorm(1, mean = my_forest$do_predict(
                        matrix(c(g_scaled[j], X_cov[i, ]), nrow = 1)),
                        sd = 1, lower = -Inf, upper = 0)
          z_store <- c(z_g, z_t)
        } else z_store <- z_t
        z_all       <- c(z_all, z_store)
        X_store     <- rbind(X_store, matrix(rep(X_cov[i, ], m_vec[i] + 1),
                                             nrow = m_vec[i] + 1, ncol = p, byrow = TRUE))
        time_points <- c(time_points, g_scaled, t_scaled[i])
        G_all_orig  <- c(G_all_orig, g_orig)
      } else {
        if (m_vec[i] > 0) {
          z_g <- numeric(m_vec[i])
          for (j in 1:m_vec[i])
            z_g[j] <- msm::rtnorm(1, mean = my_forest$do_predict(
                        matrix(c(g_scaled[j], X_cov[i, ]), nrow = 1)),
                        sd = 1, lower = -Inf, upper = 0)
          z_store <- z_g
          X_i     <- matrix(rep(X_cov[i, ], m_vec[i]), nrow = m_vec[i], ncol = p, byrow = TRUE)
        } else { z_store <- NULL;  X_i <- NULL }
        z_all       <- c(z_all, z_store)
        X_store     <- rbind(X_store, X_i)
        time_points <- c(time_points, g_scaled)
        G_all_orig  <- c(G_all_orig, g_orig)
      }
    }

    if (!is.null(z_all) && length(z_all) >= 100 &&
        length(unique(round(z_all, 6))) > 10) {
      my_forest$do_gibbs(cbind(time_points, X_store), z_all,
                         cbind(time_points, X_store), 1)
    }

    w_vec       <- w[cluster]
    event_times <- t_orig[status == 1]

    loglik_alpha <- function(a) {
      if (a <= 0) return(-Inf)
      logl <- sum(dgamma(event_times, shape = a, rate = theta, log = TRUE))
      if (length(G_all_orig) > 0)
        logl <- logl + sum(dgamma(G_all_orig, shape = a, rate = theta, log = TRUE))
      logl <- logl - beta * sum(w_vec * pgamma(t_orig, shape = a, rate = theta))
      logl <- logl + dgamma(a, shape = shape_prior_alpha, rate = rate_prior_alpha, log = TRUE)
      if (is.nan(logl) || is.na(logl)) return(-Inf);  logl
    }
    alpha <- slice_sample(loglik_alpha, alpha, w_step = 0.5, lower = 1e-4, upper = 50)

    loglik_theta <- function(th) {
      if (th <= 0) return(-Inf)
      logl <- sum(dgamma(event_times, shape = alpha, rate = th, log = TRUE))
      if (length(G_all_orig) > 0)
        logl <- logl + sum(dgamma(G_all_orig, shape = alpha, rate = th, log = TRUE))
      logl <- logl - beta * sum(w_vec * pgamma(t_orig, shape = alpha, rate = th))
      logl <- logl + dgamma(th, shape = shape_prior_theta, rate = rate_prior_theta, log = TRUE)
      if (is.nan(logl) || is.na(logl)) return(-Inf);  logl
    }
    theta <- slice_sample(loglik_theta, theta, w_step = 0.5, lower = 1e-4, upper = 50)

    shape_b <- shape_prior_b + sum(m_vec) + no_exact
    rate_b  <- rate_prior_b  + sum(w_vec * pgamma(t_orig, shape = alpha, rate = theta))
    beta    <- rgamma(1, shape = shape_b, rate = rate_b)

    for (k in 1:N_clust) {
      idx_k <- which(cluster == k)
      w[k]  <- rgamma(1,
                 shape = eta + sum(status[idx_k]) + sum(m_vec[idx_k]),
                 rate  = eta + beta * sum(pgamma(t_orig[idx_k], shape = alpha, rate = theta)))
    }
    w_vec <- w[cluster]

    loglik_eta <- function(e) {
      if (e <= 0) return(-Inf)
      logl <- sum(dgamma(w, shape = e, rate = e, log = TRUE)) +
              dgamma(e, shape = shape_prior_eta, rate = rate_prior_eta, log = TRUE)
      if (is.nan(logl) || is.na(logl)) return(-Inf);  logl
    }
    eta <- slice_sample(loglik_eta, eta, w_step = 1.0, lower = 1e-4, upper = 1000)

    beta_save[iter] <- beta; alpha_save[iter] <- alpha
    theta_save[iter] <- theta; eta_save[iter] <- eta; w_save[iter, ] <- w

    if (iter > num_burn) {
      s_idx <- iter - num_burn
      T_mc_orig_pred   <- rgamma(n_mc_pred, shape = alpha, rate = theta)
      T_mc_scaled_pred <- scale01_fixed(T_mc_orig_pred, t_min, t_max)
      X_rep <- matrix(0, nrow = n * n_mc_pred, ncol = p + 1)
      for (i in 1:n) {
        idx <- ((i - 1) * n_mc_pred + 1):(i * n_mc_pred)
        X_rep[idx, 1] <- T_mc_scaled_pred
        X_rep[idx, 2:(p + 1)] <- matrix(rep(X_cov[i, ], n_mc_pred), nrow = n_mc_pred, byrow = TRUE)
      }
      phi_big <- pnorm(as.numeric(my_forest$do_predict(X_rep)))
      phi_mat <- matrix(phi_big, nrow = n, ncol = n_mc_pred, byrow = TRUE)
      for (i in 1:n) {
        W_i_cur <- w[cluster[i]]
        surv_post[s_idx, i, ] <- sapply(t_grid_pred, function(tg)
          if (tg <= 0) 1.0 else exp(-W_i_cur * beta * mean(phi_mat[i, ] * (T_mc_orig_pred <= tg))))
        cure_post[s_idx, i] <- exp(-W_i_cur * beta * mean(phi_mat[i, ]))
        y_i <- t_orig[i]
        surv_obs_post[s_idx, i] <- if (y_i <= 0) 1.0 else
          exp(-W_i_cur * beta * mean(phi_mat[i, ] * (T_mc_orig_pred <= y_i)))
      }
    }

    if (iter %% 200 == 0)
      cat(sprintf("    Iter %d/%d (%.1f min) beta=%.3f alpha=%.3f theta=%.3f\n",
          iter, num_iter, as.numeric(difftime(Sys.time(), start_time, units = "mins")),
          beta, alpha, theta))
  }
  cat("  SBART MCMC done.\n")

  list(surv_post = surv_post, cure_post = cure_post, surv_obs_post = surv_obs_post,
       t_grid_pred = t_grid_pred)
}

# SBART metric: posterior means of the stored survival/cure arrays vs shared truth.
# Requires sbart_fit$t_grid_pred == truth$t_grid (both use the same dat and G=100).
sbart_metric <- function(dat, sbart_fit, truth) {
  stopifnot(isTRUE(all.equal(sbart_fit$t_grid_pred, truth$t_grid)))
  Sg_hat <- apply(sbart_fit$surv_post, c(2, 3), mean)
  So_hat <- colMeans(sbart_fit$surv_obs_post)
  pc_hat <- colMeans(sbart_fit$cure_post)
  n <- nrow(dat);  mc <- mo <- mk <- numeric(n)
  for (s in 1:n) {
    mc[s] <- mean((Sg_hat[s, ] - truth$Sg[s, ])^2)
    mo[s] <- (So_hat[s] - truth$So[s])^2
    mk[s] <- (pc_hat[s] - truth$pc[s])^2
  }
  c(MSEcurve = mean(mc), MSEobs = mean(mo), MSEcure = mean(mk))
}

# =============================================================================
# SECTION 4: M-CIS  (Chen et al. 2002 promotion-time cure, Gamma frailty)
#
#   S(t | W, x) = exp(-W e^{x'b} F(t)),  F = Weibull CDF,  cure = exp(-W e^{x'b})
#
# Chen-Ibrahim-Sinha (1999) data augmentation with TOTAL latent causes
#   N_ij ~ Poisson(W_i e^{x'b}); complete-data factor S(y)^{N-delta}(N f(y))^delta.
# Full conditionals use the survival S(y)=1-F(y) (NOT F), and the Poisson
# normaliser e^{-W e^{x'b}} (NOT e^{-W e^{x'b} F}).
# =============================================================================

weibull_cdf <- function(t, a, s) 1 - exp(-(t / s)^a)
weibull_pdf <- function(t, a, s) (a / s) * (t / s)^(a - 1) * exp(-(t / s)^a)

logpost_weibull_cis <- function(log_a, log_s, N_vec, y_vec, delta_vec) {
  a <- exp(log_a);  s <- exp(log_s)
  F_vec <- weibull_cdf(y_vec, a, s);  f_vec <- weibull_pdf(y_vec, a, s);  S_vec <- 1 - F_vec
  ll <- sum(delta_vec * log(pmax(f_vec, 1e-300))) +
        sum((N_vec - delta_vec) * log(pmax(S_vec, 1e-300)))   # S(y)^(N - delta)
  ll + dnorm(log_a, 0, 1, log = TRUE) + dnorm(log_s, 0, 2, log = TRUE) + log_a + log_s
}

run_cis_mcmc <- function(dat, num_iter = 5000, num_burn = 2500) {
  N <- length(unique(dat$cluster));  n_subj <- nrow(dat);  p <- 3
  X_mat <- cbind(1, dat$x1, dat$x2);  clust <- dat$cluster
  beta_cur <- c(0, 0, 0);  alpha_W_cur <- 1.5;  sigma_cur <- 1.0
  eta_cur <- 2.0;  W_cur <- rgamma(N, shape = eta_cur, rate = eta_cur)

  n_post <- num_iter - num_burn
  beta_post <- matrix(NA, n_post, p)
  alpha_W_post <- sigma_post <- eta_post <- numeric(n_post)
  W_post <- matrix(NA, n_post, N)
  mh_beta_sd <- 0.05;  mh_weib_sd <- 0.05;  acc_beta <- 0;  acc_weib <- 0
  cat("  Running M-CIS MCMC...\n")

  for (iter in 1:num_iter) {
    F_y   <- weibull_cdf(dat$y, alpha_W_cur, sigma_cur);  S_y <- 1 - F_y
    theta <- as.numeric(exp(X_mat %*% beta_cur));  W_subj <- W_cur[clust]

    # Step 1: latent total causes  N_i = delta_i + Poisson(W theta S(y))
    N_vec <- rpois(n_subj, W_subj * theta * S_y) + dat$delta

    # Step 2: W_i ~ Gamma(eta + sum N, eta + sum theta)        [theta, not theta*F]
    for (i in 1:N) {
      idx <- which(clust == i)
      W_cur[i] <- rgamma(1, shape = eta_cur + sum(N_vec[idx]),
                            rate  = eta_cur + sum(theta[idx]))
    }
    W_subj <- W_cur[clust]

    # Step 3: beta (RW-MH)   loglik = sum(N log theta - W theta)   [no F]
    beta_prop  <- beta_cur + rnorm(p, 0, mh_beta_sd)
    theta_prop <- as.numeric(exp(X_mat %*% beta_prop))
    ll_cur  <- sum(N_vec * log(pmax(theta, 1e-300))      - W_subj * theta)      +
               sum(dnorm(beta_cur,  0, 10, log = TRUE))
    ll_prop <- sum(N_vec * log(pmax(theta_prop, 1e-300)) - W_subj * theta_prop) +
               sum(dnorm(beta_prop, 0, 10, log = TRUE))
    if (log(runif(1)) < ll_prop - ll_cur) {
      beta_cur <- beta_prop;  theta <- theta_prop;  acc_beta <- acc_beta + 1
    }

    # Step 4: Weibull (RW-MH)  loglik = sum(delta log f) + sum((N-delta) log S)
    la <- log(alpha_W_cur);  ls <- log(sigma_cur)
    la_p <- la + rnorm(1, 0, mh_weib_sd);  ls_p <- ls + rnorm(1, 0, mh_weib_sd)
    lp_cur  <- logpost_weibull_cis(la,   ls,   N_vec, dat$y, dat$delta)
    lp_prop <- logpost_weibull_cis(la_p, ls_p, N_vec, dat$y, dat$delta)
    if (log(runif(1)) < lp_prop - lp_cur) {
      alpha_W_cur <- exp(la_p);  sigma_cur <- exp(ls_p);  acc_weib <- acc_weib + 1
    }

    # Step 5: eta (RW-MH on log scale)
    log_fc_eta <- function(le) {
      et <- exp(le)
      sum(dgamma(W_cur, shape = et, rate = et, log = TRUE)) +
        dgamma(et, shape = 4, rate = 0.01, log = TRUE) + le
    }
    le_cur <- log(eta_cur);  le_prop <- le_cur + rnorm(1, 0, 0.2)
    if (log(runif(1)) < log_fc_eta(le_prop) - log_fc_eta(le_cur)) eta_cur <- exp(le_prop)

    if (iter > num_burn) {
      k <- iter - num_burn
      beta_post[k, ] <- beta_cur;  alpha_W_post[k] <- alpha_W_cur
      sigma_post[k]  <- sigma_cur; eta_post[k]     <- eta_cur;  W_post[k, ] <- W_cur
    }
    if (iter %% 100 == 0 && iter <= num_burn) {
      rb <- acc_beta / iter;  rw <- acc_weib / iter
      if (rb > 0.44) mh_beta_sd <- mh_beta_sd * 1.2
      if (rb < 0.23) mh_beta_sd <- mh_beta_sd * 0.8
      if (rw > 0.44) mh_weib_sd <- mh_weib_sd * 1.2
      if (rw < 0.23) mh_weib_sd <- mh_weib_sd * 0.8
    }
  }
  cat(sprintf("  M-CIS done. beta acc=%.2f, weibull acc=%.2f\n",
              acc_beta / num_iter, acc_weib / num_iter))
  list(beta_post = beta_post, alpha_W_post = alpha_W_post, sigma_post = sigma_post,
       eta_post = eta_post, W_post = W_post, n_post = n_post)
}

cis_survival <- function(t_grid, s, dat, fit) {
  x_s <- c(1, dat$x1[s], dat$x2[s]);  cs <- dat$cluster[s]
  S <- sapply(1:fit$n_post, function(r) {
    th <- exp(sum(x_s * fit$beta_post[r, ]));  Wr <- fit$W_post[r, cs]
    exp(-Wr * th * weibull_cdf(t_grid, fit$alpha_W_post[r], fit$sigma_post[r]))
  })
  if (is.matrix(S)) rowMeans(S) else mean(S)
}
cis_cure <- function(s, dat, fit) {
  x_s <- c(1, dat$x1[s], dat$x2[s]);  cs <- dat$cluster[s]
  mean(sapply(1:fit$n_post, function(r) exp(-fit$W_post[r, cs] * exp(sum(x_s * fit$beta_post[r, ])))))
}

# =============================================================================
# SECTION 5: WGF  -- Weibull PH + Gamma frailty via parfm (Sahu et al. 1997)
#
#   h(t | W, x) = W * rho*lambda*t^(rho-1) * exp(x'b)   (parfm Weibull, NO intercept)
#   Lambda_0(t) = lambda * t^rho
#   S(t | W, x) = exp(-W * lambda * t^rho * exp(x'b))
# parfm gives MLEs; predict() returns empirical-Bayes per-cluster frailties.
# No cure fraction -> MSEcure = NA.
# =============================================================================

wgf_fit <- function(dat) {
  fit <- parfm(Surv(y, delta) ~ x1 + x2, cluster = "cluster", data = dat,
               dist = "weibull", frailty = "gamma")
  est <- fit[, "ESTIMATE"]                 # named: rho, lambda, theta(frailty var), x1, x2
  pw    <- predict(fit)                    # per-cluster EB frailties (mean-1 scale)
  W_hat <- setNames(as.numeric(pw), names(pw))
  list(rho = unname(est["rho"]), lambda = unname(est["lambda"]),
       b1 = unname(est["x1"]),   b2 = unname(est["x2"]),
       W_hat = W_hat)
}
wgf_survival <- function(t_grid, s, dat, P) {
  cl <- as.character(dat$cluster[s]);  Wi <- P$W_hat[[cl]]
  lp <- P$b1 * dat$x1[s] + P$b2 * dat$x2[s]
  Lam0 <- P$lambda * t_grid^P$rho;  Lam0[t_grid <= 0] <- 0
  exp(-Wi * Lam0 * exp(lp))
}

# =============================================================================
# SECTION 6: SGF  -- semiparametric Cox + Gamma frailty via frailtyEM (Sargent 98)
#
#   lambda(t | W, x) = W * lambda_0(t) * exp(x'b),  lambda_0 = Breslow baseline
#   S(t | W, x) = exp(-W * Lambda_0(t) * exp(x'b))
# emfrail gives EM estimates; summary()$frail gives EB per-cluster frailties.
# Baseline Lambda_0(t) is obtained once via predict() at covariates = 0.
# No cure fraction -> MSEcure = NA.
# =============================================================================

sgf_fit <- function(dat) {
  fit <- emfrail(Surv(y, delta) ~ x1 + x2 + cluster(cluster), data = dat,
                 distribution = emfrail_dist("gamma"))
  fr  <- summary(fit)$frail                # data.frame with cluster id and frailty 'z'
  id_col <- if ("id" %in% names(fr)) "id" else names(fr)[1]
  z_col  <- if ("z"  %in% names(fr)) "z"  else
            names(fr)[which(sapply(fr, is.numeric))[1]]
  W_hat  <- setNames(fr[[z_col]], as.character(fr[[id_col]]))
  beta   <- coef(fit)                       # named x1, x2 (log-HR); coef(), NOT $coefficients
  if (is.null(beta) || length(beta) < 2) stop("could not extract coefficients from emfrail fit")
  # Baseline cumulative hazard Lambda_0(t) (conditional, frailty=1, lp=0)
  base <- predict(fit, newdata = data.frame(x1 = 0, x2 = 0),
                  quantity = "cumhaz", type = "conditional")
  tcol <- if ("time" %in% names(base)) "time" else names(base)[1]
  hcol <- names(base)[grep("cumhaz", names(base))[1]]
  list(beta = beta, W_hat = W_hat,
       base_time = base[[tcol]], base_cumhaz = base[[hcol]])
}
sgf_survival <- function(t_grid, s, dat, P) {
  cl <- as.character(dat$cluster[s]);  Wi <- P$W_hat[[cl]]
  lp <- sum(P$beta * c(dat$x1[s], dat$x2[s]))
  L0 <- approx(P$base_time, P$base_cumhaz, xout = t_grid, method = "constant", rule = 2)$y
  L0[t_grid <= 0] <- 0
  exp(-Wi * L0 * exp(lp))
}

# =============================================================================
# SECTION 7: RSF  -- random survival forest (no cure, no clustering)
#   Per the paper, clustering is NOT used (cluster id is not a predictor).
#   Cure is a surrogate: S evaluated at a large time point.
# =============================================================================

fit_rsf <- function(dat) {
  rfsrc(Surv(y, delta) ~ x1 + x2,
        data = data.frame(y = dat$y, delta = dat$delta, x1 = dat$x1, x2 = dat$x2),
        ntree = 500, nodesize = 15, nsplit = 10, importance = FALSE)
}
rsf_survival <- function(t_grid, s, dat, rf_fit) {
  S <- approx(rf_fit$time.interest, rf_fit$survival[s, ],
              xout = t_grid, method = "constant", rule = 2)$y
  S[t_grid <= 0] <- 1.0;  S
}
rsf_cure <- function(s, dat, rf_fit) {
  t_large <- qgamma(0.9999, shape = alpha_true, rate = theta_true)
  rsf_survival(t_large, s, dat, rf_fit)
}

# =============================================================================
# SECTION 8: One-replicate comparison
# =============================================================================

compare_models <- function(dat, beta, sim_type = "multiplicative", G = 100,
                           sbart_iter = 5000, sbart_burn = 2500,
                           cis_iter = 5000, cis_burn = 2500) {

  truth <- make_truth(dat, beta, G = G, sim_type = sim_type)
  out   <- list()

  cat("Fitting SBART-Cure...\n")
  sb <- run_sbart_mcmc(dat, num_iter = sbart_iter, num_burn = sbart_burn)
  out[["SBART_Cure"]] <- sbart_metric(dat, sb, truth)

  cat("Fitting RSF...\n")
  rf <- fit_rsf(dat)
  out[["RSF"]] <- metric_model(dat, truth,
                    function(t, s) rsf_survival(t, s, dat, rf),
                    function(s)    rsf_cure(s, dat, rf))

  cat("Fitting WGF (parfm)...\n")
  Pw <- tryCatch(wgf_fit(dat), error = function(e) { message("  WGF failed: ", e$message); NULL })
  out[["WGF"]] <- if (is.null(Pw)) c(MSEcurve = NA, MSEobs = NA, MSEcure = NA) else
                    metric_model(dat, truth, function(t, s) wgf_survival(t, s, dat, Pw), NULL)

  cat("Fitting SGF (frailtyEM)...\n")
  Ps <- tryCatch(sgf_fit(dat), error = function(e) { message("  SGF failed: ", e$message); NULL })
  out[["SGF"]] <- if (is.null(Ps)) c(MSEcurve = NA, MSEobs = NA, MSEcure = NA) else
                    metric_model(dat, truth, function(t, s) sgf_survival(t, s, dat, Ps), NULL)

  cat("Fitting M-CIS...\n")
  mc <- run_cis_mcmc(dat, num_iter = cis_iter, num_burn = cis_burn)
  out[["M_CIS"]] <- metric_model(dat, truth,
                      function(t, s) cis_survival(t, s, dat, mc),
                      function(s)    cis_cure(s, dat, mc))

  do.call(rbind, lapply(names(out), function(nm)
    data.frame(Model = nm,
               MSEcurve = out[[nm]]["MSEcurve"],
               MSEobs   = out[[nm]]["MSEobs"],
               MSEcure  = out[[nm]]["MSEcure"], row.names = NULL)))
}

# =============================================================================
# SECTION 9: Parallel simulation study (4 scenarios x R replicates)
# =============================================================================

run_one_replicate <- function(r, sim) {
  set.seed(sim * 1000 + r)
  settings <- list(
    `1` = list(N = 50, lambda_ni = 20, frailty = "gamma",    beta_d = beta_true,   lam = lambda_cens_1, sim_type = "multiplicative"),
    `2` = list(N = 50, lambda_ni = 20, frailty = "stable",   beta_d = beta_stable, lam = lambda_cens_2, sim_type = "multiplicative"),
    `3` = list(N = 50, lambda_ni = 20, frailty = "additive", beta_d = beta_add,    lam = lambda_cens_3, sim_type = "additive"),
    `4` = list(N = 50, lambda_ni = 3,  frailty = "gamma",    beta_d = beta_true,   lam = lambda_cens_1, sim_type = "multiplicative")
  )
  s   <- settings[[as.character(sim)]]
  dat <- generate_data(N = s$N, lambda_ni = s$lambda_ni, frailty = s$frailty,
                       beta_for_data = s$beta_d, lambda_cens = s$lam)
  tryCatch(
    compare_models(dat, beta = s$beta_d, sim_type = s$sim_type),
    error = function(e) { cat(sprintf("  ERROR sim=%d rep=%d: %s\n", sim, r, e$message)); NULL })
}

run_simulation_parallel <- function(R, n_cores = detectCores() - 2) {
  jobs <- expand.grid(r = 1:R, sim = 1:4)
  cat(sprintf("Total jobs: %d (%d reps x 4 sims) on %d cores\n", nrow(jobs), R, n_cores))
  all_results <- mclapply(seq_len(nrow(jobs)), function(j) {
    r <- jobs$r[j];  sim <- jobs$sim[j]
    res <- run_one_replicate(r, sim)
    if (!is.null(res)) { res$Replicate <- r; res$Simulation <- sim }
    list(sim = sim, r = r, result = res)
  }, mc.cores = n_cores, mc.silent = FALSE)

  agg <- function(sim_num, label) {
    rl <- Filter(Negate(is.null),
                 lapply(Filter(function(x) x$sim == sim_num, all_results), function(x) x$result))
    if (length(rl) == 0) return(NULL)
    combined <- do.call(rbind, rl)
    a <- aggregate(cbind(MSEcurve, MSEobs, MSEcure) ~ Model, data = combined,
                   FUN = function(v) mean(v, na.rm = TRUE), na.action = na.pass)
    a$Simulation <- label;  a
  }
  list(
    table_sim1 = agg(1, "Sim1_MultGamma"),
    table_sim2 = agg(2, "Sim2_MultStable"),
    table_sim3 = agg(3, "Sim3_AddGamma"),
    table_sim4 = agg(4, "Sim4_SmallClust")
  )
}

# --- Run -----------------------------------------------------------------------
R <- 100
sim_results <- run_simulation_parallel(R = R)
table1 <- do.call(rbind, sim_results)
print(table1, row.names = FALSE)
write.csv(table1, "table1_simulation_results.csv", row.names = FALSE)
