## =========================
## FINAL SCRIPT (Strict route + joint lambda estimation)
## =========================
set.seed(2025)
suppressPackageStartupMessages({
  library(mvtnorm)
  library(stats)
  library(pROC)
  library(ggplot2)
})

## =========================
## 0) Preparation: Data generation + MAR missingness
## =========================

# ---- Parameters ----
n <- 500   # sample size
d <- 2     # dimension
pi <- c(0.5, 0.5)          # mixing proportions
mu <- list(c(-1,0), c(1,0))  # means of two components
Sigma <- list(diag(2), diag(2)) # covariance matrices

# ---- Posterior probability function (for data filtering only) ----
posterior_prob <- function(y, pi, mu, Sigma){
  f1 <- dmvnorm(y, mean=mu[[1]], sigma=Sigma[[1]]) * pi[1]
  f2 <- dmvnorm(y, mean=mu[[2]], sigma=Sigma[[2]]) * pi[2]
  f1 / (f1 + f2)
}

# ---- Margin confidence function ----
margin_conf <- function(y) {
  tau1 <- posterior_prob(y, pi, mu, Sigma)
  abs(2*tau1 - 1)
}

# ---- Generate data with margin confidence in [0,0.6] (hard set) ----
sample_data <- function(n, pi, mu, Sigma){
  out <- matrix(NA, nrow=0, ncol=3)
  while(nrow(out) < n){
    z <- sample(1:2, 1, prob=pi)
    y <- rmvnorm(1, mean=mu[[z]], sigma=Sigma[[z]])
    mc <- margin_conf(y)
    if(mc >= 0 && mc <= 0.6){
      out <- rbind(out, c(y,z))
    }
  }
  colnames(out) <- c("x1","x2","label")
  as.data.frame(out)
}

data <- sample_data(n, pi, mu, Sigma)

# ---- Add margin confidence ----
data$mc <- apply(data[, c("x1","x2")], 1, margin_conf)

# ---- Aranda–Ordaz missingness (asymmetric form) ----
ao_inv <- function(eta, lambda) {
  q <- 1 - (1 + lambda * exp(eta))^(-1/lambda)
  pmin(pmax(q, 1e-12), 1 - 1e-12)
}

q_aranda_ordaz <- function(delta2, alpha0, alpha1, lambda){
  ao_inv(alpha0 + alpha1 * delta2, lambda)
}

calibrate_alpha0 <- function(delta2, alpha1, lambda, target){
  f <- function(a0){ mean(q_aranda_ordaz(delta2, a0, alpha1, lambda)) - target }
  uniroot(f, c(-50, 50))$root
}

delta2  <- data$mc^2
lambda0 <- 0.5
alpha1_0 <- -6
target_missing <- 0.7

alpha0_0 <- calibrate_alpha0(delta2, alpha1_0, lambda0, target_missing)
q_miss <- q_aranda_ordaz(delta2, alpha0_0, alpha1_0, lambda0)

# ---- Store true labels & apply missingness ----
data$label_true    <- as.integer(data$label)                 # True labels (kept)
data$label_missing <- rbinom(nrow(data), 1, q_miss)          # Missing indicator
data$q_miss        <- q_miss
data$label[data$label_missing == 1] <- NA_integer_           # Observed labels (with NA)

## =========================
## 1) Building blocks for ECM (strict route)
## =========================

# Vectorized MVN density
dmvprod <- function(Y, mu, Sigma) {
  # Y: n x d matrix, mu: d-vector, Sigma: d x d covariance
  mvtnorm::dmvnorm(Y, mean = mu, sigma = Sigma)
}

# Responsibilities and margin delta (K=2)
compute_tau_delta <- function(Y, pis, mus, Sigmas) {
  f1 <- dmvprod(Y, mus[[1]], Sigmas[[1]]) * pis[1]
  f2 <- dmvprod(Y, mus[[2]], Sigmas[[2]]) * pis[2]
  denom <- pmax(f1 + f2, 1e-300)
  tau1 <- f1 / denom
  tau2 <- f2 / denom
  delta <- abs(2*tau1 - 1)  # binary margin confidence
  list(tau = cbind(tau1, tau2), delta = delta)
}

# Incomplete-data log-likelihood for mixture
ig_loglik <- function(Y, pis, mus, Sigmas, z_obs) {
  f1 <- dmvprod(Y, mus[[1]], Sigmas[[1]]) * pis[1]
  f2 <- dmvprod(Y, mus[[2]], Sigmas[[2]]) * pis[2]
  ll <- numeric(nrow(Y))
  obs_id <- which(!is.na(z_obs))
  mis_id <- which(is.na(z_obs))
  if(length(obs_id) > 0){
    k1 <- which(z_obs[obs_id] == 1)
    k2 <- which(z_obs[obs_id] == 2)
    if(length(k1) > 0) ll[obs_id[k1]] <- log(pmax(f1[obs_id[k1]], 1e-300))
    if(length(k2) > 0) ll[obs_id[k2]] <- log(pmax(f2[obs_id[k2]], 1e-300))
  }
  if(length(mis_id) > 0){
    ll[mis_id] <- log(pmax(f1[mis_id] + f2[mis_id], 1e-300))
  }
  sum(ll)
}

# Missing-mechanism negative log-likelihood for alpha (given lambda, delta2)
miss_negloglik_alpha <- function(par_alpha, delta2, m, lambda) {
  q <- ao_inv(par_alpha[1] + par_alpha[2]*delta2, lambda)
  -sum( (1-m)*log(1-q) + m*log(q) )
}

## ---- pack/unpack theta (unconstrained) for strict CM1 ----
# theta layout: [s | mu1(d) | mu2(d) | L1(lower incl diag) | L2(lower incl diag)]
pack_theta <- function(pis, mus, Sigmas) {
  d <- length(mus[[1]])
  s <- qlogis(pis[1])  # pis[2]=1-pis[1]
  L1 <- chol(Sigmas[[1]])
  L2 <- chol(Sigmas[[2]])
  L1v <- c(L1[lower.tri(L1, diag=TRUE)])
  L2v <- c(L2[lower.tri(L2, diag=TRUE)])
  # log the diagonal entries for positivity
  log_diag <- function(Lv){
    # diagonal positions in lower-tri vector are 1, 1+2, 1+2+3, ...
    idx <- cumsum(1:length(diag(matrix(0, nrow=length(mus[[1]]), ncol=length(mus[[1]])))))
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
  # rebuild lower-tri with exp on diagonal
  toL <- function(Lv) {
    L <- matrix(0, d, d)
    L[lower.tri(L, diag=TRUE)] <- Lv
    for (i in 1:d) L[i,i] <- exp(L[i,i])
    L
  }
  L1 <- toL(L1v); L2 <- toL(L2v)
  S1 <- L1 %*% t(L1); S2 <- L2 %*% t(L2)
  pi1 <- plogis(s); pis <- c(pi1, 1 - pi1)
  list(pis=pis, mus=list(mu1,mu2), Sigmas=list(S1,S2))
}

# Joint objective (strict CM1): -(ig_loglik + ll_miss) to minimize w.r.t theta
joint_negloglik_theta <- function(theta_vec, Y, z_obs, m, alpha0, alpha1, lambda) {
  d <- ncol(Y)
  pars <- unpack_theta(theta_vec, d)
  comp <- compute_tau_delta(Y, pars$pis, pars$mus, pars$Sigmas)
  tau  <- comp$tau
  delta <- comp$delta
  if (any(m==0L)) {
    tau[cbind(which(m==0L), z_obs[m==0L])] <- 1
    tau[cbind(which(m==0L), 3 - z_obs[m==0L])] <- 0
  }
  ll_ig   <- ig_loglik(Y, pars$pis, pars$mus, pars$Sigmas, z_obs)
  q_now   <- ao_inv(alpha0 + alpha1 * (delta^2), lambda)
  ll_miss <- sum((1 - m)*log(1 - q_now) + m*log(q_now))
  -(ll_ig + ll_miss)
}

## =========================
## 2) ECM with strict CM1 + lambda estimation
## =========================
ECM_GMM_MAR <- function(Y, z_obs,
                        lambda_init = 0.5, estimate_lambda = TRUE,
                        max_iter = 200, tol = 1e-6, ridge = 1e-6, verbose = TRUE) {
  n <- nrow(Y); d <- ncol(Y); K <- 2
  stopifnot(K == 2)
  m <- ifelse(is.na(z_obs), 1L, 0L)  # missing indicator (1=missing/0=observed)
  
  # ---- init theta via kmeans ----
  km <- kmeans(Y, centers = K, nstart = 10)
  pis <- as.numeric(table(km$cluster)/n)
  mus <- list(colMeans(Y[km$cluster==1, , drop=FALSE]),
              colMeans(Y[km$cluster==2, , drop=FALSE]))
  Sigmas <- list(cov(Y[km$cluster==1, , drop=FALSE]) + diag(ridge, d),
                 cov(Y[km$cluster==2, , drop=FALSE]) + diag(ridge, d))
  
  # ---- init missingness params ----
  lambda <- lambda_init
  comp0 <- compute_tau_delta(Y, pis, mus, Sigmas)
  delta2 <- comp0$delta^2
  alpha1 <- -4
  target <- mean(m)
  cal_fn <- function(a0) mean(ao_inv(a0 + alpha1*delta2, lambda)) - target
  alpha0 <- tryCatch(uniroot(cal_fn, c(-50,50))$root, error=function(e) 0)
  
  loglik_trace <- numeric(max_iter)
  for (t in 1:max_iter) {
    
    ## ======== E-step ========
    comp <- compute_tau_delta(Y, pis, mus, Sigmas)
    tau  <- comp$tau
    delta <- comp$delta
    delta2 <- delta^2
    # enforce one-hot for labeled samples
    if(any(m==0L)){
      tau[cbind(which(m==0L), z_obs[m==0L])] <- 1
      tau[cbind(which(m==0L), 3 - z_obs[m==0L])] <- 0
    }
    
    ## ======== CM-step 1 (STRICT): maximize ig_loglik + ll_miss wrt theta ========
    theta0 <- pack_theta(pis, mus, Sigmas)
    opt_theta <- optim(
      par = theta0,
      fn  = joint_negloglik_theta,
      Y = Y, z_obs = z_obs, m = m,
      alpha0 = alpha0, alpha1 = alpha1, lambda = lambda,
      method = "BFGS",
      control = list(maxit = 300, reltol = 1e-7)
    )
    pars_upd <- unpack_theta(opt_theta$par, d)
    pis    <- pars_upd$pis
    mus    <- pars_upd$mus
    Sigmas <- pars_upd$Sigmas
    
    ## ======== CM-step 2: update (alpha, lambda) ========
    # First re-compute delta2 with updated theta
    comp <- compute_tau_delta(Y, pis, mus, Sigmas)
    delta2 <- comp$delta^2
    
    # optimize alpha given current lambda
    opt_alpha <- optim(par = c(alpha0, alpha1),
                       fn = miss_negloglik_alpha,
                       delta2 = delta2, m = m, lambda = lambda,
                       method = "BFGS", control = list(maxit = 300))
    alpha0_new <- opt_alpha$par[1]; alpha1_new <- opt_alpha$par[2]
    
    if (estimate_lambda) {
      # profile over t = log(lambda) for numerical stability
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
      
      # re-fit alpha with updated lambda
      opt_alpha <- optim(par = c(alpha0_new, alpha1_new),
                         fn = miss_negloglik_alpha,
                         delta2 = delta2, m = m, lambda = lambda,
                         method = "BFGS", control = list(maxit = 300))
      alpha0_new <- opt_alpha$par[1]; alpha1_new <- opt_alpha$par[2]
    }
    alpha0 <- alpha0_new; alpha1 <- alpha1_new
    
    ## ======== log-likelihood trace (for convergence monitor) ========
    ll_ig   <- ig_loglik(Y, pis, mus, Sigmas, z_obs)
    q_now   <- ao_inv(alpha0 + alpha1*delta2, lambda)
    ll_miss <- sum( (1-m)*log(1 - q_now) + m*log(q_now) )
    loglik_trace[t] <- ll_ig + ll_miss
    
    if (verbose && (t %% 5 == 0)) {
      message(sprintf("Iter %d | loglik = %.4f | lambda = %.4f | alpha=(%.4f, %.4f) | pi=(%.3f, %.3f)",
                      t, loglik_trace[t], lambda, alpha0, alpha1, pis[1], pis[2]))
    }
    if (t > 1 && abs(loglik_trace[t] - loglik_trace[t-1]) < tol) {
      loglik_trace <- loglik_trace[1:t]
      break
    }
  }
  
  # Final E for outputs
  comp <- compute_tau_delta(Y, pis, mus, Sigmas)
  tau  <- comp$tau
  delta <- comp$delta
  z_impute <- apply(tau, 1, which.max)
  
  list(pis = pis, mus = mus, Sigmas = Sigmas,
       alpha = c(alpha0 = alpha0, alpha1 = alpha1),
       lambda = lambda,
       tau = tau, delta = delta,
       z_impute = z_impute,
       loglik = loglik_trace)
}

## =========================
## 3) Run ECM and evaluate
## =========================
Y <- as.matrix(data[, c("x1","x2")])
z_obs <- data$label          # With NA
z_true <- data$label_true    # True labels

fit <- ECM_GMM_MAR(Y, z_obs,
                   lambda_init = 0.5,
                   estimate_lambda = TRUE,   # jointly estimate lambda
                   max_iter = 200, tol = 1e-6, verbose = TRUE)

# Posterior prob for class 1
p_hat <- fit$tau[,1]
y_true01 <- as.integer(z_true == 1)

# ROC & AUC
roc_obj <- roc(response = y_true01, predictor = p_hat, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

# Use both fixed threshold and Youden-optimal threshold
thr_fixed <- 0.5
thr_best  <- as.numeric(coords(roc_obj, "best", best.method = "youden", ret = "threshold"))

metrics_at_thr <- function(thr) {
  y_pred01 <- as.integer(p_hat >= thr)
  TP <- sum(y_pred01==1 & y_true01==1)
  TN <- sum(y_pred01==0 & y_true01==0)
  FP <- sum(y_pred01==1 & y_true01==0)
  FN <- sum(y_pred01==0 & y_true01==1)
  ACC <- (TP + TN) / length(y_true01)
  Precision <- ifelse((TP+FP)>0, TP/(TP+FP), NA)
  Recall    <- ifelse((TP+FN)>0, TP/(TP+FN), NA)
  F1 <- ifelse(is.na(Precision)|is.na(Recall)|(Precision+Recall)==0,
               NA, 2*Precision*Recall/(Precision+Recall))
  list(threshold = thr, ACC=ACC, Precision=Precision, Recall=Recall, F1=F1,
       CM = matrix(c(TN, FP, FN, TP), nrow=2, byrow=TRUE,
                   dimnames=list('True' = c('Class0','Class1'),
                                 'Pred' = c('Class0','Class1'))))
}

m_fix  <- metrics_at_thr(thr_fixed)
m_best <- metrics_at_thr(thr_best)

# LogLoss & Brier
eps <- 1e-12
LogLoss <- -mean( y_true01*log(pmin(p_hat,1-eps)) + (1-y_true01)*log(pmin(1-p_hat,1-eps)) )
Brier   <- mean( (p_hat - y_true01)^2 )

cat("==== ECM (STRICT) Evaluation ====\n")
cat(sprintf("AUC               : %.4f\n", auc_val))
cat(sprintf("LogLoss           : %.4f\n", LogLoss))
cat(sprintf("Brier             : %.4f\n\n", Brier))

cat(sprintf("[Threshold = %.2f]  ACC=%.4f  Precision=%.4f  Recall=%.4f  F1=%.4f\n",
            m_fix$threshold, m_fix$ACC, m_fix$Precision, m_fix$Recall, m_fix$F1))
print(m_fix$CM)

cat(sprintf("\n[Youden best = %.3f] ACC=%.4f  Precision=%.4f  Recall=%.4f  F1=%.4f\n",
            m_best$threshold, m_best$ACC, m_best$Precision, m_best$Recall, m_best$F1))
print(m_best$CM)

# ROC curve
ggplot(data.frame(fpr = 1-roc_obj$specificities,
                  tpr = roc_obj$sensitivities),
       aes(x=fpr, y=tpr)) +
  geom_path(size=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title=sprintf("ROC Curve (AUC = %.3f)", auc_val),
       x="False Positive Rate", y="True Positive Rate") +
  theme_minimal(base_size = 12)

## =========================
## 4) Print fitted parameters
## =========================
cat("\n===== Fitted mixture params (STRICT) =====\n")
cat("pi: ", paste(round(fit$pis,4), collapse=", "), "\n")
cat("mu1:", paste(round(fit$mus[[1]],4), collapse=", "), "\n")
cat("mu2:", paste(round(fit$mus[[2]],4), collapse=", "), "\n")
cat("alpha0, alpha1:", paste(round(fit$alpha,4), collapse=", "), "\n")
cat("lambda:", round(fit$lambda,4), "\n")







## =========================
## Baseline: 忽略缺失机制，只用有标签样本训练
## =========================

# 1) 拆分数据
Y <- as.matrix(data[, c("x1","x2")])
z_obs <- data$label          # 含 NA
z_true <- data$label_true    # 真标签

obs_id <- which(!is.na(z_obs))
mis_id <- which(is.na(z_obs))

Y_train <- Y[obs_id, , drop=FALSE]
z_train <- z_obs[obs_id]
Y_test  <- Y[mis_id, , drop=FALSE]
z_test  <- z_true[mis_id]

# 2) 用有标签的数据估计参数
pi_hat <- as.numeric(table(z_train) / length(z_train))
mu_hat <- lapply(1:2, function(k) colMeans(Y_train[z_train==k, , drop=FALSE]))
Sigma_hat <- lapply(1:2, function(k) cov(Y_train[z_train==k, , drop=FALSE]))

# 3) 对所有数据（含缺失标签部分）计算 posterior prob
f1 <- dmvnorm(Y, mean=mu_hat[[1]], sigma=Sigma_hat[[1]]) * pi_hat[1]
f2 <- dmvnorm(Y, mean=mu_hat[[2]], sigma=Sigma_hat[[2]]) * pi_hat[2]
p_hat <- f1 / (f1 + f2)

# 4) 模型评估
y_true01 <- as.integer(z_true == 1)

roc_obj <- roc(response=y_true01, predictor=p_hat, quiet=TRUE)
auc_val <- as.numeric(auc(roc_obj))

thr_fixed <- 0.5
thr_best  <- as.numeric(coords(roc_obj, "best", best.method="youden", ret="threshold"))

metrics_at_thr <- function(thr) {
  y_pred01 <- as.integer(p_hat >= thr)
  TP <- sum(y_pred01==1 & y_true01==1)
  TN <- sum(y_pred01==0 & y_true01==0)
  FP <- sum(y_pred01==1 & y_true01==0)
  FN <- sum(y_pred01==0 & y_true01==1)
  ACC <- (TP+TN)/length(y_true01)
  Precision <- ifelse((TP+FP)>0, TP/(TP+FP), NA)
  Recall    <- ifelse((TP+FN)>0, TP/(TP+FN), NA)
  F1 <- ifelse(is.na(Precision)|is.na(Recall)|(Precision+Recall)==0,
               NA, 2*Precision*Recall/(Precision+Recall))
  list(threshold=thr, ACC=ACC, Precision=Precision, Recall=Recall, F1=F1,
       CM=matrix(c(TN,FP,FN,TP), nrow=2, byrow=TRUE,
                 dimnames=list(True=c("Class0","Class1"),
                               Pred=c("Class0","Class1"))))
}

m_fix  <- metrics_at_thr(thr_fixed)
m_best <- metrics_at_thr(thr_best)

eps <- 1e-12
LogLoss <- -mean(y_true01*log(pmin(p_hat,1-eps)) + (1-y_true01)*log(pmax(1-p_hat,eps)))
Brier   <- mean((p_hat - y_true01)^2)

cat("\n==== Baseline (No missingness modeling) ====\n")
cat(sprintf("AUC               : %.4f\n", auc_val))
cat(sprintf("LogLoss           : %.4f\n", LogLoss))
cat(sprintf("Brier             : %.4f\n\n", Brier))

cat(sprintf("[Threshold = %.2f]  ACC=%.4f  Precision=%.4f  Recall=%.4f  F1=%.4f\n",
            m_fix$threshold, m_fix$ACC, m_fix$Precision, m_fix$Recall, m_fix$F1))
print(m_fix$CM)

cat(sprintf("\n[Youden best = %.3f] ACC=%.4f  Precision=%.4f  Recall=%.4f  F1=%.4f\n",
            m_best$threshold, m_best$ACC, m_best$Precision, m_best$Recall, m_best$F1))
print(m_best$CM)

# ROC 曲线
ggplot(data.frame(fpr=1-roc_obj$specificities,
                  tpr=roc_obj$sensitivities),
       aes(x=fpr,y=tpr)) +
  geom_path(size=1) +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  labs(title=sprintf("ROC Curve (Baseline, AUC=%.3f)", auc_val),
       x="False Positive Rate", y="True Positive Rate") +
  theme_minimal(base_size=12)












simulate_and_eval <- function(reps = 20, estimate_lambda = FALSE, seed = 2025){
  set.seed(seed)
  out <- matrix(NA_real_, nrow = reps, ncol = 6,
                dimnames = list(NULL, c("AUC","ACC","Precision","Recall","F1","LogLoss")))
  for (r in 1:reps){
    data <- sample_data(n, pi, mu, Sigma)
    data$mc <- apply(data[, c("x1","x2")], 1, margin_conf)
    delta2  <- data$mc^2
    alpha0_0 <- calibrate_alpha0(delta2, alpha1_0, lambda0, target_missing)
    q_miss <- q_aranda_ordaz(delta2, alpha0_0, alpha1_0, lambda0)
    data$label_true <- as.integer(data$label)
    data$label_missing <- rbinom(nrow(data), 1, q_miss)
    data$q_miss <- q_miss
    data$label[data$label_missing==1] <- NA_integer_
    
    Y <- as.matrix(data[, c("x1","x2")])
    z_obs <- data$label
    z_true <- data$label_true
    
    fit <- ECM_GMM_MAR(Y, z_obs, lambda_init = lambda0,
                       estimate_lambda = estimate_lambda,
                       max_iter = 300, tol = 1e-6, verbose = FALSE)
    
    p_hat <- fit$tau[,1]
    y_true01 <- as.integer(z_true==1)
    eps <- 1e-12
    # metrics
    roc_obj <- pROC::roc(response = y_true01, predictor = p_hat, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    y_pred01 <- as.integer(p_hat >= 0.5)
    TP <- sum(y_pred01==1 & y_true01==1)
    TN <- sum(y_pred01==0 & y_true01==0)
    FP <- sum(y_pred01==1 & y_true01==0)
    FN <- sum(y_pred01==0 & y_true01==1)
    ACC <- (TP+TN)/length(y_true01)
    Precision <- ifelse((TP+FP)>0, TP/(TP+FP), NA)
    Recall <- ifelse((TP+FN)>0, TP/(TP+FN), NA)
    F1 <- ifelse(is.na(Precision)|is.na(Recall)|(Precision+Recall)==0,
                 NA, 2*Precision*Recall/(Precision+Recall))
    LogLoss <- -mean( y_true01*log(pmin(pmax(p_hat,eps),1-eps)) +
                        (1-y_true01)*log(pmin(pmax(1-p_hat,eps),1-eps)) )
    out[r,] <- c(auc_val, ACC, Precision, Recall, F1, LogLoss)
  }
  # 汇总
  means <- colMeans(out, na.rm=TRUE)
  sds   <- apply(out, 2, sd, na.rm=TRUE)
  summary <- data.frame(metric = names(means), mean = means, sd = sds, row.names = NULL)
  list(raw = out, summary = summary)
}




set.seed(2026)
exp_fix  <- simulate_and_eval(reps=20, estimate_lambda = FALSE)
set.seed(2026)
exp_free <- simulate_and_eval(reps=20, estimate_lambda = TRUE)

print(exp_fix$summary)
print(exp_free$summary)

#> print(exp_fix$summary)
#metric      mean         sd
#1       AUC 0.7037677 0.04240444
#2       ACC 0.6583000 0.03974273
#3 Precision 0.6682489 0.03147902
#4    Recall 0.6292017 0.15359713
#5        F1 0.6645173 0.02865068
#6   LogLoss 0.6485786 0.09066817
#> print(exp_free$summary)
#metric      mean         sd
#1       AUC 0.7037212 0.04248708
#2       ACC 0.6583000 0.03970563
#3 Precision 0.6684098 0.03125658
#4    Recall 0.6289027 0.15330994
#5        F1 0.6644513 0.02808808
#6   LogLoss 0.6500592 0.09585404