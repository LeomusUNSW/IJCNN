set.seed(2025)
suppressPackageStartupMessages({
  library(mvtnorm)
  library(stats)
  library(pROC)
  library(ggplot2)
  library(dplyr)
})

## =========================
## 0) Preparation: Data generation + MAR missingness
## =========================

# ---- Parameters ----
n <- 200   # sample size
d <- 2     # dimension
pi <- c(0.5, 0.5)  # mixing proportions

# Parameters for two components of a 2D Gamma mixture
# (Each component has independent Gamma distributions on x1 and x2 with different parameters)
param1 <- list(shape = c(3, 3), scale = c(1, 1))  # Cluster 1 gamma parameters for (x1, x2)
param2 <- list(shape = c(5, 5), scale = c(1, 1))  # Cluster 2 gamma parameters for (x1, x2)

# ---- Posterior probability function (for data filtering only) under Gamma mixture ----
posterior_prob <- function(y, pi, param1, param2) {
  # Compute the mixture component densities for point y = (x1, x2)
  dens1 <- dgamma(y[1], shape = param1$shape[1], scale = param1$scale[1]) *
    dgamma(y[2], shape = param1$shape[2], scale = param1$scale[2])
  dens2 <- dgamma(y[1], shape = param2$shape[1], scale = param2$scale[1]) *
    dgamma(y[2], shape = param2$shape[2], scale = param2$scale[2])
  # Weighted by mixing proportions
  f1 <- dens1 * pi[1]
  f2 <- dens2 * pi[2]
  f1 / (f1 + f2)
}

# ---- Margin confidence function (based on posterior probability) ----
margin_conf <- function(y) {
  tau1 <- posterior_prob(y, pi, param1, param2)
  abs(2*tau1 - 1)
}

# ---- Generate data with margin confidence in [0, 0.6] (hard set) ----
sample_data <- function(n, pi, param1, param2) {
  out <- matrix(NA, nrow = 0, ncol = 3)
  while (nrow(out) < n) {
    z <- sample(1:2, 1, prob = pi)  # sample component label
    # Generate a 2D point from the corresponding Gamma distribution
    if (z == 1) {
      y1 <- rgamma(1, shape = param1$shape[1], scale = param1$scale[1])
      y2 <- rgamma(1, shape = param1$shape[2], scale = param1$scale[2])
    } else {
      y1 <- rgamma(1, shape = param2$shape[1], scale = param2$scale[1])
      y2 <- rgamma(1, shape = param2$shape[2], scale = param2$scale[2])
    }
    y <- c(y1, y2)
    mc <- margin_conf(y)
    # Keep the sample only if margin confidence is between 0 and 0.6
    if (mc >= 0 && mc <= 0.6) {
      out <- rbind(out, c(y, z))
    }
  }
  colnames(out) <- c("x1", "x2", "label")
  as.data.frame(out)
}

# Generate the data
data <- sample_data(n, pi, param1, param2)

# ---- Add margin confidence to data ----
data$mc <- apply(data[, c("x1", "x2")], 1, margin_conf)

# ---- Aranda–Ordaz missingness (asymmetric form) ----
ao_inv <- function(eta, lambda) {
  q <- 1 - (1 + lambda * exp(eta))^(-1/lambda)
  pmin(pmax(q, 1e-12), 1 - 1e-12)
}
q_aranda_ordaz <- function(delta2, alpha0, alpha1, lambda) {
  ao_inv(alpha0 + alpha1 * delta2, lambda)
}
calibrate_alpha0 <- function(delta2, alpha1, lambda, target) {
  f <- function(a0) { mean(q_aranda_ordaz(delta2, a0, alpha1, lambda)) - target }
  uniroot(f, c(-50, 50))$root
}

delta2 <- data$mc^2
lambda0 <- 2
alpha1_0 <- -6
target_missing <- 0.7
alpha0_0 <- calibrate_alpha0(delta2, alpha1_0, lambda0, target_missing)
q_miss <- q_aranda_ordaz(delta2, alpha0_0, alpha1_0, lambda0)

# ---- Store true labels & apply missingness ----
data$label_true   <- as.integer(data$label)          # True labels (keep for evaluation)
data$label_missing <- rbinom(nrow(data), 1, q_miss)  # Missing indicator (1 = label missing)
data$q_miss        <- q_miss
data$label[data$label_missing == 1] <- NA_integer_   # Observed labels (with NA for missing)

## =========================
## 1) Building blocks for ECM (strict route)
## =========================

# Vectorized MVN density (for Gaussian mixture calculations in ECM algorithm)
dmvprod <- function(Y, mu, Sigma) { mvtnorm::dmvnorm(Y, mean = mu, sigma = Sigma) }

# Responsibilities and margin delta (K = 2)
compute_tau_delta <- function(Y, pis, mus, Sigmas) {
  f1 <- dmvprod(Y, mus[[1]], Sigmas[[1]]) * pis[1]
  f2 <- dmvprod(Y, mus[[2]], Sigmas[[2]]) * pis[2]
  denom <- pmax(f1 + f2, 1e-300)
  tau1 <- f1 / denom
  tau2 <- f2 / denom
  delta <- abs(2*tau1 - 1)
  list(tau = cbind(tau1, tau2), delta = delta)
}

# Incomplete-data log-likelihood for mixture
ig_loglik <- function(Y, pis, mus, Sigmas, z_obs) {
  f1 <- dmvprod(Y, mus[[1]], Sigmas[[1]]) * pis[1]
  f2 <- dmvprod(Y, mus[[2]], Sigmas[[2]]) * pis[2]
  ll <- numeric(nrow(Y))
  obs_id <- which(!is.na(z_obs))
  mis_id <- which(is.na(z_obs))
  if (length(obs_id) > 0) {
    k1 <- which(z_obs[obs_id] == 1)
    k2 <- which(z_obs[obs_id] == 2)
    if (length(k1) > 0) ll[obs_id[k1]] <- log(pmax(f1[obs_id[k1]], 1e-300))
    if (length(k2) > 0) ll[obs_id[k2]] <- log(pmax(f2[obs_id[k2]], 1e-300))
  }
  if (length(mis_id) > 0) {
    ll[mis_id] <- log(pmax(f1[mis_id] + f2[mis_id], 1e-300))
  }
  sum(ll)
}

# Missing-mechanism negative log-likelihood for alpha (given lambda, delta2) -- AO link
miss_negloglik_alpha <- function(par_alpha, delta2, m, lambda) {
  q <- ao_inv(par_alpha[1] + par_alpha[2] * delta2, lambda)
  -sum((1-m)*log(1-q) + m*log(q))
}

## ---- pack/unpack theta (unconstrained) for strict CM1 ----
# theta layout: [s | mu1(d) | mu2(d) | L1(lower incl diag) | L2(lower incl diag)]
pack_theta <- function(pis, mus, Sigmas) {
  d <- length(mus[[1]])
  s <- qlogis(pis[1])  # pis[2] = 1 - pis[1]
  L1 <- chol(Sigmas[[1]])
  L2 <- chol(Sigmas[[2]])
  L1v <- c(L1[lower.tri(L1, diag=TRUE)])
  L2v <- c(L2[lower.tri(L2, diag=TRUE)])
  # log the diagonal entries for positivity
  log_diag <- function(Lv) {
    idx <- cumsum(1:d)  # diagonal positions
    Lv[idx] <- log(Lv[idx])
    Lv
  }
  L1v <- log_diag(L1v)
  L2v <- log_diag(L2v)
  c(s, mus[[1]], mus[[2]], L1v, L2v)
}

unpack_theta <- function(theta_vec, d) {
  s <- theta_vec[1]
  off <- 2
  mu1 <- theta_vec[off:(off+d-1)]; off <- off + d
  mu2 <- theta_vec[off:(off+d-1)]; off <- off + d
  Llen <- d*(d+1)/2
  L1v <- theta_vec[off:(off+Llen-1)]; off <- off + Llen
  L2v <- theta_vec[off:(off+Llen-1)]
  
  toL <- function(Lv) {
    L <- matrix(0, d, d)
    L[lower.tri(L, diag=TRUE)] <- Lv
    for (i in 1:d) L[i,i] <- exp(L[i,i])
    L
  }
  L1 <- toL(L1v); L2 <- toL(L2v)
  S1 <- L1 %*% t(L1); S2 <- L2 %*% t(L2)
  pi1 <- plogis(s); pis <- c(pi1, 1 - pi1)
  list(pis = pis, mus = list(mu1, mu2), Sigmas = list(S1, S2))
}

# Joint objective (strict CM1): -(ig_loglik + ll_miss) w.r.t theta (AO link)
joint_negloglik_theta <- function(theta_vec, Y, z_obs, m, alpha0, alpha1, lambda) {
  d <- ncol(Y)
  pars <- unpack_theta(theta_vec, d)
  comp <- compute_tau_delta(Y, pars$pis, pars$mus, pars$Sigmas)
  tau <- comp$tau
  delta <- comp$delta
  if (any(m == 0L)) {
    tau[cbind(which(m == 0L), z_obs[m == 0L])] <- 1
    tau[cbind(which(m == 0L), 3 - z_obs[m == 0L])] <- 0
  }
  ll_ig <- ig_loglik(Y, pars$pis, pars$mus, pars$Sigmas, z_obs)
  q_now <- ao_inv(alpha0 + alpha1 * (delta^2), lambda)
  ll_miss <- sum((1 - m)*log(1 - q_now) + m*log(q_now))
  -(ll_ig + ll_miss)
}

## =========================
## 2) ECM with strict CM1 + lambda estimation (AO link)
## =========================
ECM_GMM_MAR <- function(Y, z_obs, lambda_init = 0.5, estimate_lambda = TRUE,
                        max_iter = 200, tol = 1e-6, ridge = 1e-6, verbose = TRUE) {
  n <- nrow(Y); d <- ncol(Y); K <- 2
  stopifnot(K == 2)
  m <- ifelse(is.na(z_obs), 1L, 0L)  # missing indicator (1 = missing, 0 = observed)
  
  # ---- init theta via kmeans ----
  km <- kmeans(Y, centers = K, nstart = 10)
  pis <- as.numeric(table(km$cluster) / n)
  mus <- list(colMeans(Y[km$cluster == 1, , drop=FALSE]),
              colMeans(Y[km$cluster == 2, , drop=FALSE]))
  Sigmas <- list(cov(Y[km$cluster == 1, , drop=FALSE]) + diag(ridge, d),
                 cov(Y[km$cluster == 2, , drop=FALSE]) + diag(ridge, d))
  
  # ---- init missingness params ----
  lambda <- lambda_init
  comp0 <- compute_tau_delta(Y, pis, mus, Sigmas)
  delta2 <- comp0$delta^2
  alpha1 <- -4
  target <- mean(m)
  cal_fn <- function(a0) mean(ao_inv(a0 + alpha1*delta2, lambda)) - target
  alpha0 <- tryCatch(uniroot(cal_fn, c(-50, 50))$root, error = function(e) 0)
  
  loglik_trace <- numeric(max_iter)
  for (t in 1:max_iter) {
    ## ======== E-step ========
    comp <- compute_tau_delta(Y, pis, mus, Sigmas)
    tau <- comp$tau
    delta <- comp$delta
    delta2 <- delta^2
    if (any(m == 0L)) {
      tau[cbind(which(m == 0L), z_obs[m == 0L])] <- 1
      tau[cbind(which(m == 0L), 3 - z_obs[m == 0L])] <- 0
    }
    
    ## ======== CM-step 1 (STRICT): maximize ig_loglik + ll_miss w.r.t. theta ========
    theta0 <- pack_theta(pis, mus, Sigmas)
    opt_theta <- optim(
      par = theta0,
      fn = joint_negloglik_theta, Y = Y, z_obs = z_obs, m = m,
      alpha0 = alpha0, alpha1 = alpha1, lambda = lambda,
      method = "BFGS", control = list(maxit = 300, reltol = 1e-7)
    )
    pars_upd <- unpack_theta(opt_theta$par, d)
    pis <- pars_upd$pis; mus <- pars_upd$mus; Sigmas <- pars_upd$Sigmas
    
    ## ======== CM-step 2: update (alpha, lambda) ========
    comp <- compute_tau_delta(Y, pis, mus, Sigmas)
    delta2 <- comp$delta^2
    
    opt_alpha <- optim(par = c(alpha0, alpha1),
                       fn = miss_negloglik_alpha,
                       delta2 = delta2, m = m, lambda = lambda,
                       method = "BFGS", control = list(maxit = 300))
    alpha0_new <- opt_alpha$par[1]; alpha1_new <- opt_alpha$par[2]
    
    if (estimate_lambda) {
      prof_fn_t <- function(tval) {
        lmb <- exp(tval)
        op <- optim(c(alpha0_new, alpha1_new), miss_negloglik_alpha,
                    delta2 = delta2, m = m, lambda = lmb,
                    method = "BFGS", control = list(maxit = 300))
        op$value
      }
      grid_t <- seq(log(0.05), log(8), length.out = 31)
      vals <- sapply(grid_t, prof_fn_t)
      t0 <- grid_t[which.min(vals)]
      opt_t <- optimize(prof_fn_t, c(t0 - 1, t0 + 1))
      lambda <- exp(opt_t$minimum)
      
      opt_alpha <- optim(par = c(alpha0_new, alpha1_new),
                         fn = miss_negloglik_alpha,
                         delta2 = delta2, m = m, lambda = lambda,
                         method = "BFGS", control = list(maxit = 300))
      alpha0_new <- opt_alpha$par[1]; alpha1_new <- opt_alpha$par[2]
    }
    
    alpha0 <- alpha0_new; alpha1 <- alpha1_new
    
    ## ======== log-likelihood trace ========
    ll_ig <- ig_loglik(Y, pis, mus, Sigmas, z_obs)
    q_now <- ao_inv(alpha0 + alpha1*delta2, lambda)
    ll_miss <- sum((1-m)*log(1 - q_now) + m*log(q_now))
    loglik_trace[t] <- ll_ig + ll_miss
    
    if (verbose && (t %% 5 == 0)) {
      message(sprintf("AO Iter %d | loglik = %.4f | lambda = %.4f | alpha=(%.4f, %.4f) | pi=(%.3f, %.3f)",
                      t, loglik_trace[t], lambda, alpha0, alpha1, pis[1], pis[2]))
    }
    if (t > 1 && abs(loglik_trace[t] - loglik_trace[t-1]) < tol) {
      loglik_trace <- loglik_trace[1:t]
      break
    }
  }
  
  comp <- compute_tau_delta(Y, pis, mus, Sigmas)
  tau <- comp$tau
  delta <- comp$delta
  z_impute <- apply(tau, 1, which.max)
  list(pis = pis, mus = mus, Sigmas = Sigmas,
       alpha = c(alpha0 = alpha0, alpha1 = alpha1),
       lambda = lambda, tau = tau, delta = delta,
       z_impute = z_impute, loglik = loglik_trace)
}

## =========================
## 2b) ECM with strict CM1 (Logit link for missingness)
## =========================

# Logistic inverse link for missingness
logistic_inv <- function(eta) {
  q <- 1 / (1 + exp(-eta))
  pmin(pmax(q, 1e-12), 1 - 1e-12)
}

# Negative log-likelihood for alpha under logit link
miss_negloglik_alpha_logit <- function(par_alpha, delta2, m) {
  q <- logistic_inv(par_alpha[1] + par_alpha[2]*delta2)
  -sum((1-m)*log(1-q) + m*log(q))
}

# Joint objective (strict CM1) for theta under logit link
joint_negloglik_theta_logit <- function(theta_vec, Y, z_obs, m, alpha0, alpha1) {
  d <- ncol(Y)
  pars <- unpack_theta(theta_vec, d)
  comp <- compute_tau_delta(Y, pars$pis, pars$mus, pars$Sigmas)
  tau <- comp$tau
  delta <- comp$delta
  if (any(m == 0L)) {
    tau[cbind(which(m == 0L), z_obs[m == 0L])] <- 1
    tau[cbind(which(m == 0L), 3 - z_obs[m == 0L])] <- 0
  }
  ll_ig <- ig_loglik(Y, pars$pis, pars$mus, pars$Sigmas, z_obs)
  q_now <- logistic_inv(alpha0 + alpha1 * (delta^2))
  ll_miss <- sum((1 - m)*log(1 - q_now) + m*log(q_now))
  -(ll_ig + ll_miss)
}

# ECM main for Logit missingness (no lambda parameter)
ECM_GMM_MAR_LOGIT <- function(Y, z_obs,
                              max_iter = 200, tol = 1e-6, ridge = 1e-6, verbose = TRUE) {
  n <- nrow(Y); d <- ncol(Y); K <- 2
  stopifnot(K == 2)
  m <- ifelse(is.na(z_obs), 1L, 0L)
  
  # init theta via kmeans
  km <- kmeans(Y, centers = K, nstart = 10)
  pis <- as.numeric(table(km$cluster) / n)
  mus <- list(colMeans(Y[km$cluster == 1, , drop=FALSE]),
              colMeans(Y[km$cluster == 2, , drop=FALSE]))
  Sigmas <- list(cov(Y[km$cluster == 1, , drop=FALSE]) + diag(ridge, d),
                 cov(Y[km$cluster == 2, , drop=FALSE]) + diag(ridge, d))
  
  # init alpha by matching mean missing rate with logistic link
  comp0 <- compute_tau_delta(Y, pis, mus, Sigmas)
  delta2 <- comp0$delta^2
  alpha1 <- -4
  target <- mean(m)
  cal_fn <- function(a0) mean(logistic_inv(a0 + alpha1*delta2)) - target
  alpha0 <- tryCatch(uniroot(cal_fn, c(-50, 50))$root, error = function(e) 0)
  
  loglik_trace <- numeric(max_iter)
  for (t in 1:max_iter) {
    # E-step
    comp <- compute_tau_delta(Y, pis, mus, Sigmas)
    tau <- comp$tau
    delta <- comp$delta
    delta2 <- delta^2
    if (any(m == 0L)) {
      tau[cbind(which(m == 0L), z_obs[m == 0L])] <- 1
      tau[cbind(which(m == 0L), 3 - z_obs[m == 0L])] <- 0
    }
    
    # CM-step 1: strict joint update of mixture theta under logit missingness
    theta0 <- pack_theta(pis, mus, Sigmas)
    opt_theta <- optim(
      par = theta0,
      fn = joint_negloglik_theta_logit,
      Y = Y, z_obs = z_obs, m = m,
      alpha0 = alpha0, alpha1 = alpha1,
      method = "BFGS", control = list(maxit = 300, reltol = 1e-7)
    )
    pars_upd <- unpack_theta(opt_theta$par, d)
    pis <- pars_upd$pis; mus <- pars_upd$mus; Sigmas <- pars_upd$Sigmas
    
    # CM-step 2: update alpha (logit link)
    comp <- compute_tau_delta(Y, pis, mus, Sigmas)
    delta2 <- comp$delta^2
    opt_alpha <- optim(par = c(alpha0, alpha1),
                       fn = miss_negloglik_alpha_logit,
                       delta2 = delta2, m = m,
                       method = "BFGS", control = list(maxit = 300))
    alpha0 <- opt_alpha$par[1]; alpha1 <- opt_alpha$par[2]
    
    # loglik trace
    ll_ig <- ig_loglik(Y, pis, mus, Sigmas, z_obs)
    q_now <- logistic_inv(alpha0 + alpha1*delta2)
    ll_miss <- sum((1-m)*log(1 - q_now) + m*log(q_now))
    loglik_trace[t] <- ll_ig + ll_miss
    
    if (verbose && (t %% 5 == 0)) {
      message(sprintf("LOGIT Iter %d | loglik = %.4f | alpha=(%.4f, %.4f) | pi=(%.3f, %.3f)",
                      t, loglik_trace[t], alpha0, alpha1, pis[1], pis[2]))
    }
    if (t > 1 && abs(loglik_trace[t] - loglik_trace[t-1]) < tol) {
      loglik_trace <- loglik_trace[1:t]
      break
    }
  }
  
  comp <- compute_tau_delta(Y, pis, mus, Sigmas)
  tau <- comp$tau
  delta <- comp$delta
  z_impute <- apply(tau, 1, which.max)
  list(pis = pis, mus = mus, Sigmas = Sigmas,
       alpha = c(alpha0 = alpha0, alpha1 = alpha1),
       tau = tau, delta = delta,
       z_impute = z_impute, loglik = loglik_trace)
}

## =========================
## 3) Run ECMs and evaluate
## =========================
Y <- as.matrix(data[, c("x1", "x2")])
z_obs <- data$label           # Observed labels (NA for missing)
z_true <- data$label_true     # True labels
y_true01 <- as.integer(z_true == 1)

# ---- ECM (AO link) ----
fit <- ECM_GMM_MAR(Y, z_obs, lambda_init = 0.5, estimate_lambda = TRUE,
                   max_iter = 500, tol = 1e-6, verbose = FALSE)
p_hat_ecm <- fit$tau[, 1]

roc_ecm <- roc(response = y_true01, predictor = p_hat_ecm, quiet = TRUE)
auc_ecm <- as.numeric(auc(roc_ecm))

metrics_at_thr <- function(p_hat, y_true01, thr) {
  y_pred01 <- as.integer(p_hat >= thr)
  TP <- sum(y_pred01 == 1 & y_true01 == 1)
  TN <- sum(y_pred01 == 0 & y_true01 == 0)
  FP <- sum(y_pred01 == 1 & y_true01 == 0)
  FN <- sum(y_pred01 == 0 & y_true01 == 1)
  ACC <- (TP + TN) / length(y_true01)
  Precision <- ifelse((TP+FP) > 0, TP/(TP+FP), NA)
  Recall <- ifelse((TP+FN) > 0, TP/(TP+FN), NA)
  F1 <- ifelse(is.na(Precision) | is.na(Recall) | (Precision+Recall) == 0,
               NA, 2*Precision*Recall/(Precision+Recall))
  list(ACC = ACC, Precision = Precision, Recall = Recall, F1 = F1,
       CM = matrix(c(TN, FP, FN, TP), nrow = 2, byrow = TRUE,
                   dimnames = list('True' = c('Class0','Class1'),
                                   'Pred' = c('Class0','Class1'))))
}
thr_fixed <- 0.5
thr_best_ecm <- as.numeric(coords(roc_ecm, "best", best.method = "youden", ret = "threshold"))
m_fix_ecm <- metrics_at_thr(p_hat_ecm, y_true01, thr_fixed)
m_best_ecm <- metrics_at_thr(p_hat_ecm, y_true01, thr_best_ecm)

eps <- 1e-12
LogLoss_ecm <- -mean(y_true01*log(pmin(p_hat_ecm, 1-eps)) + (1-y_true01)*log(pmin(1-p_hat_ecm, 1-eps)))
Brier_ecm <- mean((p_hat_ecm - y_true01)^2)

cat("==== ECM (STRICT, AO missingness) Evaluation ====\n")
cat(sprintf("AUC : %.4f\n", auc_ecm))
cat(sprintf("LogLoss : %.4f\n", LogLoss_ecm))
cat(sprintf("Brier : %.4f\n\n", Brier_ecm))
cat(sprintf("[Threshold = %.2f] ACC=%.4f Precision=%.4f Recall=%.4f F1=%.4f\n",
            thr_fixed, m_fix_ecm$ACC, m_fix_ecm$Precision, m_fix_ecm$Recall, m_fix_ecm$F1))
print(m_fix_ecm$CM)
cat(sprintf("\n[Youden best = %.3f] ACC=%.4f Precision=%.4f Recall=%.4f F1=%.4f\n",
            thr_best_ecm, m_best_ecm$ACC, m_best_ecm$Precision, m_best_ecm$Recall, m_best_ecm$F1))
print(m_best_ecm$CM)

cat("\n===== Fitted mixture params (STRICT AO) =====\n")
cat("pi: ", paste(round(fit$pis, 4), collapse=", "), "\n")
cat("mu1:", paste(round(fit$mus[[1]], 4), collapse=", "), "\n")
cat("mu2:", paste(round(fit$mus[[2]], 4), collapse=", "), "\n")
cat("alpha0, alpha1:", paste(round(fit$alpha, 4), collapse=", "), "\n")
cat("lambda:", round(fit$lambda, 4), "\n")

# ---- ECM (Logit link) ----
fit_logit <- ECM_GMM_MAR_LOGIT(Y, z_obs, max_iter = 500, tol = 1e-6, verbose = FALSE)
p_hat_ecm_logit <- fit_logit$tau[, 1]

roc_ecm_logit <- roc(response = y_true01, predictor = p_hat_ecm_logit, quiet = TRUE)
auc_ecm_logit <- as.numeric(auc(roc_ecm_logit))
thr_best_ecm_logit <- as.numeric(coords(roc_ecm_logit, "best", best.method = "youden", ret = "threshold"))
m_fix_ecm_logit <- metrics_at_thr(p_hat_ecm_logit, y_true01, thr_fixed)
m_best_ecm_logit <- metrics_at_thr(p_hat_ecm_logit, y_true01, thr_best_ecm_logit)

LogLoss_ecm_logit <- -mean(y_true01*log(pmin(p_hat_ecm_logit, 1-eps)) + (1-y_true01)*log(pmin(1-p_hat_ecm_logit, 1-eps)))
Brier_ecm_logit <- mean((p_hat_ecm_logit - y_true01)^2)

cat("\n==== ECM (STRICT, Logit missingness) Evaluation ====\n")
cat(sprintf("AUC : %.4f\n", auc_ecm_logit))
cat(sprintf("LogLoss : %.4f\n", LogLoss_ecm_logit))
cat(sprintf("Brier : %.4f\n\n", Brier_ecm_logit))
cat(sprintf("[Threshold = %.2f] ACC=%.4f Precision=%.4f Recall=%.4f F1=%.4f\n",
            thr_fixed, m_fix_ecm_logit$ACC, m_fix_ecm_logit$Precision, m_fix_ecm_logit$Recall, m_fix_ecm_logit$F1))
print(m_fix_ecm_logit$CM)
cat(sprintf("\n[Youden best = %.3f] ACC=%.4f Precision=%.4f Recall=%.4f F1=%.4f\n",
            thr_best_ecm_logit, m_best_ecm_logit$ACC, m_best_ecm_logit$Precision, m_best_ecm_logit$Recall, m_best_ecm_logit$F1))
print(m_best_ecm_logit$CM)

cat("\n===== Fitted mixture params (STRICT Logit) =====\n")
cat("pi: ", paste(round(fit_logit$pis, 4), collapse=", "), "\n")
cat("mu1:", paste(round(fit_logit$mus[[1]], 4), collapse=", "), "\n")
cat("mu2:", paste(round(fit_logit$mus[[2]], 4), collapse=", "), "\n")
cat("alpha0, alpha1:", paste(round(fit_logit$alpha, 4), collapse=", "), "\n")

## =========================
## 4) Logistic regression baseline (ignoring missingness)
## =========================
obs_idx <- which(!is.na(data$label))
Y_train <- data[obs_idx, c("x1", "x2")]
z_train <- factor(data$label_true[obs_idx], levels = c(1, 2))

logit_model <- glm(z_train ~ x1 + x2, data = Y_train, family = binomial)
p_hat_logit <- predict(logit_model, newdata = data[, c("x1", "x2")], type = "response")

roc_logit <- roc(response = y_true01, predictor = p_hat_logit, quiet = TRUE)
auc_logit <- as.numeric(auc(roc_logit))
thr_best_logit <- as.numeric(coords(roc_logit, "best", best.method = "youden", ret = "threshold"))
m_fix_logit <- metrics_at_thr(p_hat_logit, y_true01, thr_fixed)
m_best_logit <- metrics_at_thr(p_hat_logit, y_true01, thr_best_logit)

LogLoss_logit <- -mean(y_true01*log(pmin(p_hat_logit, 1-eps)) + (1-y_true01)*log(pmin(1-p_hat_logit, 1-eps)))
Brier_logit <- mean((p_hat_logit - y_true01)^2)

cat("\n==== Logistic Regression (ignoring missingness) Evaluation ====\n")
cat(sprintf("AUC : %.4f\n", auc_logit))
cat(sprintf("LogLoss : %.4f\n", LogLoss_logit))
cat(sprintf("Brier : %.4f\n\n", Brier_logit))
cat(sprintf("[Threshold = %.2f] ACC=%.4f Precision=%.4f Recall=%.4f F1=%.4f\n",
            thr_fixed, m_fix_logit$ACC, m_fix_logit$Precision, m_fix_logit$Recall, m_fix_logit$F1))
print(m_fix_logit$CM)
cat(sprintf("\n[Youden best = %.3f] ACC=%.4f Precision=%.4f Recall=%.4f F1=%.4f\n",
            thr_best_logit, m_best_logit$ACC, m_best_logit$Precision, m_best_logit$Recall, m_best_logit$F1))
print(m_best_logit$CM)

## =========================
## 4b) One combined ROC plot (all three methods)
## =========================
library(RColorBrewer)

df_roc_ecm <- data.frame(
  fpr = 1 - roc_ecm$specificities,
  tpr = roc_ecm$sensitivities,
  Method = sprintf("ECM-AO (AUC = %.3f)", auc_ecm)
)
df_roc_ecm_logit <- data.frame(
  fpr = 1 - roc_ecm_logit$specificities,
  tpr = roc_ecm_logit$sensitivities,
  Method = sprintf("ECM-Logit (AUC = %.3f)", auc_ecm_logit)
)
df_roc_logit <- data.frame(
  fpr = 1 - roc_logit$specificities,
  tpr = roc_logit$sensitivities,
  Method = sprintf("Logistic Baseline (AUC = %.3f)", auc_logit)
)

df_roc_all <- bind_rows(df_roc_ecm, df_roc_ecm_logit, df_roc_logit)
cols <- brewer.pal(3, "Set1")

print(
  ggplot(df_roc_all, aes(x = fpr, y = tpr, color = Method)) +
    geom_path(size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    scale_color_manual(values = cols) +
    coord_equal() +
    labs(title = "ROC Curves: ECM-AO vs ECM-Logit vs Logistic",
         x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank())
)

## =========================
## 5) Summarize into comparison table
## =========================
extract_metrics <- function(p_hat, y_true, roc_obj, thr_fixed = 0.5) {
  auc_val <- as.numeric(auc(roc_obj))
  eps <- 1e-12
  LogLoss <- -mean(y_true*log(pmax(p_hat, eps)) + (1-y_true)*log(pmax(1-p_hat, eps)))
  Brier <- mean((p_hat - y_true)^2)
  
  # fixed threshold metrics
  y_pred_fix <- ifelse(p_hat >= thr_fixed, 1, 0)
  ACC_fix <- mean(y_pred_fix == y_true)
  Precision_fix <- sum(y_pred_fix == 1 & y_true == 1) / max(1, sum(y_pred_fix == 1))
  Recall_fix <- sum(y_pred_fix == 1 & y_true == 1) / max(1, sum(y_true == 1))
  F1_fix <- ifelse(Precision_fix + Recall_fix == 0, NA,
                   2 * Precision_fix * Recall_fix / (Precision_fix + Recall_fix))
  
  # optimal threshold (Youden's index)
  thr_opt <- as.numeric(coords(roc_obj, "best", best.method = "youden", ret = "threshold"))
  y_pred_opt <- ifelse(p_hat >= thr_opt, 1, 0)
  ACC_opt <- mean(y_pred_opt == y_true)
  Precision_opt <- sum(y_pred_opt == 1 & y_true == 1) / max(1, sum(y_pred_opt == 1))
  Recall_opt <- sum(y_pred_opt == 1 & y_true == 1) / max(1, sum(y_true == 1))
  F1_opt <- ifelse(Precision_opt + Recall_opt == 0, NA,
                   2 * Precision_opt * Recall_opt / (Precision_opt + Recall_opt))
  
  data.frame(
    AUC = auc_val, LogLoss = LogLoss, Brier = Brier,
    ACC_0.5 = ACC_fix, Precision_0.5 = Precision_fix, Recall_0.5 = Recall_fix, F1_0.5 = F1_fix,
    ACC_opt = ACC_opt, Precision_opt = Precision_opt, Recall_opt = Recall_opt, F1_opt = F1_opt,
    row.names = NULL
  )
}

metrics_ecm       <- extract_metrics(p_hat_ecm,       y_true01, roc_ecm)
metrics_ecm_logit <- extract_metrics(p_hat_ecm_logit, y_true01, roc_ecm_logit)
metrics_logit     <- extract_metrics(p_hat_logit,     y_true01, roc_logit)

result_table <- bind_rows(
  cbind(Method = "ECM (with missingness, AO link)",    metrics_ecm),
  cbind(Method = "ECM (with missingness, Logit link)", metrics_ecm_logit),
  cbind(Method = "Logistic (ignoring missingness)",    metrics_logit)
)

cat("\n==== Comparison Table (ECM-AO vs ECM-Logit vs Logistic) ====\n")
print(result_table)
