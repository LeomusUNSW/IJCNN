## ============================================================
## Banknote MAR(70%) —— Fast ECM-AO / ECM-Logit / Logistic baseline
## 关键优化：
##  - CM1 使用 GMM 的闭式 M-step（不再对 θ 做 BFGS）
##  - 仅对缺失模型参数 α（可选 λ）做低维优化
##  - 仍保留：协方差正定化、权重裁剪、数值护栏
## ============================================================

set.seed(2025)
suppressPackageStartupMessages({
  library(mvtnorm)
  library(stats)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(RColorBrewer)
  library(tidyr)
  library(scales)
})

## ========= 0) 读取与预处理 =========
dat <- read.csv("magic07.csv")
X_raw <- as.matrix(dat[, c("fAlpha","fLength","fM3Long","fSize")])
X <- scale(X_raw)
z_obs  <- dat$class          # 观测标签(含 NA；1/2)
z_true <- dat$class_true     # 真实标签(1/2)
y_true01 <- as.integer(z_true == 1)

## ========= 1) 工具函数（含数值护栏） =========
clip01 <- function(p, eps=1e-16) pmin(pmax(p, eps), 1-eps)

safe_dmvnorm <- function(Y, mean, sigma) {
  out <- try(mvtnorm::dmvnorm(Y, mean=mean, sigma=sigma), silent=TRUE)
  if (inherits(out, "try-error") || any(!is.finite(out))) rep(1e-300, nrow(Y)) else out
}

sanitize_sym <- function(S, eps_diag = 1e-8) {
  d <- ncol(S); S[!is.finite(S)] <- 0; S <- (S + t(S)) / 2
  diag(S) <- pmax(diag(S), eps_diag); S
}

weighted_cov <- function(X, w = NULL, ridge = 1e-6) {
  X <- as.matrix(X)
  ok <- rowSums(!is.finite(X)) == 0
  X <- X[ok, , drop = FALSE]
  if (!is.null(w)) w <- pmax(as.numeric(w[ok]), 0)
  n <- nrow(X); d <- ncol(X)
  if (n < 2 || d < 1) return(diag(1, d))
  if (is.null(w)) {
    S <- stats::cov(X)
  } else {
    sw <- sum(w); if (!is.finite(sw) || sw <= 1e-12) return(diag(1, d))
    mu <- colSums(X * w) / sw
    Xc <- sweep(X, 2, mu, "-")
    S  <- (t(Xc) %*% (Xc * w)) / sw
  }
  S <- sanitize_sym(S) + diag(ridge, ncol(S)); S
}

ensure_pd <- function(S, eps = 1e-8, max_try = 5) {
  S <- sanitize_sym(S, eps)
  ev <- try(eigen(S, symmetric = TRUE), silent = TRUE)
  if (!inherits(ev, "try-error")) {
    lam <- pmax(ev$values, eps); U <- ev$vectors
    S <- U %*% (lam * t(U)); S <- sanitize_sym(S, eps)
  } else S <- S + diag(eps, ncol(S))
  add <- eps; k <- 0
  while (inherits(try(chol(S), silent=TRUE), "try-error") && k < max_try) {
    add <- add * 10; S <- S + diag(add, ncol(S)); k <- k + 1
  }
  S
}

ao_inv <- function(eta, lambda) {
  lam <- max(lambda, 1e-8)
  q <- 1 - (1 + lam * exp(eta))^(-1/lam)
  clip01(q, 1e-12)
}
logistic_inv <- function(eta) clip01(1/(1 + exp(-eta)), 1e-12)

## ========= 2) 共享：E步（计算 τ 与 margin） =========
E_step <- function(X, pi1, mu1, S1, mu0, S0, z_obs) {
  f1 <- safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8)
  f0 <- safe_dmvnorm(X, mu0, S0) * clip01(1-pi1,1e-8)
  den <- pmax(f1 + f0, 1e-300)
  tau1 <- clip01(f1 / den, 1e-12)
  # Clamp 已标注样本
  obs_id <- which(!is.na(z_obs))
  if (length(obs_id)) {
    lab1 <- (z_obs[obs_id] == 1)
    tau1[obs_id] <- as.numeric(lab1)
  }
  tau0 <- 1 - tau1
  delta <- abs(2*tau1 - 1)
  list(tau1=tau1, tau0=tau0, delta=delta)
}

## ========= 3) 共享：CM1（闭式 M-step 更新 π, μ, Σ） =========
M_step_params <- function(X, tau1, tau0, ridge=1e-4) {
  tau1 <- clip01(tau1, 1e-12); tau0 <- clip01(tau0, 1e-12)
  w1 <- sum(tau1); w0 <- sum(tau0); sw <- w1 + w0
  w1 <- max(w1, 1e-6); w0 <- max(w0, 1e-6)
  pi1 <- clip01(w1 / sw, 1e-6)
  mu1 <- colSums(X * tau1) / w1
  mu0 <- colSums(X * tau0) / w0
  S1  <- ensure_pd( weighted_cov(X, tau1, ridge) )
  S0  <- ensure_pd( weighted_cov(X, tau0, ridge) )
  list(pi1=pi1, mu1=mu1, S1=S1, mu0=mu0, S0=S0)
}

## ========= 4) 共享：CM2（更新缺失模型参数 α；可选 λ） =========
negloglik_miss_AO <- function(par, delta2, m, lambda) {
  q <- ao_inv(par[1] + par[2]*delta2, lambda)
  val <- -sum((1-m)*log(clip01(1-q)) + m*log(clip01(q)))
  if (!is.finite(val)) 1e12 else val
}
negloglik_miss_logit <- function(par, delta2, m) {
  q <- logistic_inv(par[1] + par[2]*delta2)
  val <- -sum((1-m)*log(clip01(1-q)) + m*log(clip01(q)))
  if (!is.finite(val)) 1e12 else val
}

## ========= 5) ECM-GMM（AO-link）—— 快速版 =========
ECM_GMM_MAR_AO_fast <- function(X, z_obs, lambda=0.7, estimate_lambda=FALSE,
                                max_iter=200, tol=1e-5, ridge=1e-4, verbose=FALSE) {
  n <- nrow(X); d <- ncol(X); m <- ifelse(is.na(z_obs), 1L, 0L)
  
  # 初始化（kmeans + 协方差正定化）
  km <- kmeans(X, centers=2, nstart=10)
  pi1 <- clip01(as.numeric(mean(km$cluster==1)), 1e-3)
  mu1 <- colMeans(X[km$cluster==1,,drop=FALSE]); mu0 <- colMeans(X[km$cluster==2,,drop=FALSE])
  S1  <- ensure_pd(weighted_cov(X[km$cluster==1,,drop=FALSE], NULL, ridge))
  S0  <- ensure_pd(weighted_cov(X[km$cluster==2,,drop=FALSE], NULL, ridge))
  
  # 初始 α：匹配缺失率（固定斜率 alpha1）
  delta0 <- E_step(X, pi1, mu1, S1, mu0, S0, z_obs)$delta
  alpha1 <- -4
  target <- mean(m)
  f_cal <- function(a0) mean(ao_inv(a0 + alpha1*(delta0^2), lambda)) - target
  alpha0 <- tryCatch(uniroot(f_cal, c(-60,60))$root, error=function(e) 0)
  
  ll_trace <- numeric(max_iter)
  for (t in 1:max_iter) {
    # E
    e <- E_step(X, pi1, mu1, S1, mu0, S0, z_obs)
    delta2 <- (e$delta)^2
    
    # CM1：闭式 M-step 更新 GMM 参数
    m1 <- M_step_params(X, e$tau1, e$tau0, ridge=ridge)
    pi1 <- m1$pi1; mu1 <- m1$mu1; S1 <- m1$S1; mu0 <- m1$mu0; S0 <- m1$S0
    
    # CM2：优化 α（2 维 BFGS，极快）
    opt_a <- optim(c(alpha0, alpha1), negloglik_miss_AO,
                   delta2=delta2, m=m, lambda=lambda, method="BFGS",
                   control=list(maxit=200, reltol=1e-8))
    alpha0 <- opt_a$par[1]; alpha1 <- opt_a$par[2]
    
    # 可选：profile λ
    if (estimate_lambda) {
      prof <- function(tval) {
        lmb <- exp(tval)
        optim(c(alpha0, alpha1), negloglik_miss_AO,
              delta2=delta2, m=m, lambda=lmb, method="BFGS",
              control=list(maxit=100, reltol=1e-8))$value
      }
      grid_t <- seq(log(0.1), log(5), length.out=11)
      vals <- sapply(grid_t, prof)
      t0 <- grid_t[which.min(vals)]
      opt_t <- optimize(prof, c(t0-0.5, t0+0.5))
      lambda <- max(exp(opt_t$minimum), 1e-3)
    }
    
    # 监控联合对数似然（观测部分）
    f1 <- safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8)
    f0 <- safe_dmvnorm(X, mu0, S0) * clip01(1-pi1,1e-8)
    den <- pmax(f1+f0, 1e-300)
    obs_id <- which(!is.na(z_obs)); mis_id <- which(is.na(z_obs))
    ll_gmm <- 0
    if (length(obs_id)) {
      ll_gmm <- ll_gmm + sum(log(ifelse(z_obs[obs_id]==1, f1[obs_id], f0[obs_id])))
    }
    if (length(mis_id)) {
      ll_gmm <- ll_gmm + sum(log(den[mis_id]))
    }
    q <- ao_inv(alpha0 + alpha1*delta2, lambda)
    ll_miss <- sum((1-m)*log(clip01(1-q)) + m*log(clip01(q)))
    ll <- ll_gmm + ll_miss
    ll_trace[t] <- ll
    if (verbose && t %% 5 == 0) message(sprintf("AO iter %3d | ll=%.4f | pi=%.3f", t, ll, pi1))
    if (t>1 && abs(ll_trace[t] - ll_trace[t-1]) < tol) { ll_trace <- ll_trace[1:t]; break }
  }
  list(pi1=pi1, mu1=mu1, S1=S1, mu0=mu0, S0=S0,
       alpha=c(alpha0=alpha0, alpha1=alpha1), lambda=lambda,
       tau = cbind(
         safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8) /
           pmax(safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8) +
                  safe_dmvnorm(X, mu0, S0) * clip01(1-pi1,1e-8), 1e-300),
         NA),
       loglik=ll_trace)
}

## ========= 6) ECM-GMM（logit-link）—— 快速版 =========
ECM_GMM_MAR_LOGIT_fast <- function(X, z_obs,
                                   max_iter=200, tol=1e-5, ridge=1e-4, verbose=FALSE) {
  n <- nrow(X); d <- ncol(X); m <- ifelse(is.na(z_obs), 1L, 0L)
  
  km <- kmeans(X, centers=2, nstart=10)
  pi1 <- clip01(as.numeric(mean(km$cluster==1)), 1e-3)
  mu1 <- colMeans(X[km$cluster==1,,drop=FALSE]); mu0 <- colMeans(X[km$cluster==2,,drop=FALSE])
  S1  <- ensure_pd(weighted_cov(X[km$cluster==1,,drop=FALSE], NULL, ridge))
  S0  <- ensure_pd(weighted_cov(X[km$cluster==2,,drop=FALSE], NULL, ridge))
  
  # 初始 α（固定 alpha1）
  delta0 <- E_step(X, pi1, mu1, S1, mu0, S0, z_obs)$delta
  alpha1 <- -4; target <- mean(m)
  f_cal <- function(a0) mean(logistic_inv(a0 + alpha1*(delta0^2))) - target
  alpha0 <- tryCatch(uniroot(f_cal, c(-60,60))$root, error=function(e) 0)
  
  ll_trace <- numeric(max_iter)
  for (t in 1:max_iter) {
    e <- E_step(X, pi1, mu1, S1, mu0, S0, z_obs)
    delta2 <- (e$delta)^2
    
    m1 <- M_step_params(X, e$tau1, e$tau0, ridge=ridge)
    pi1 <- m1$pi1; mu1 <- m1$mu1; S1 <- m1$S1; mu0 <- m1$mu0; S0 <- m1$S0
    
    opt_a <- optim(c(alpha0, alpha1), negloglik_miss_logit,
                   delta2=delta2, m=m, method="BFGS",
                   control=list(maxit=200, reltol=1e-8))
    alpha0 <- opt_a$par[1]; alpha1 <- opt_a$par[2]
    
    f1 <- safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8)
    f0 <- safe_dmvnorm(X, mu0, S0) * clip01(1-pi1,1e-8)
    den <- pmax(f1+f0, 1e-300)
    obs_id <- which(!is.na(z_obs)); mis_id <- which(is.na(z_obs))
    ll_gmm <- 0
    if (length(obs_id)) ll_gmm <- ll_gmm + sum(log(ifelse(z_obs[obs_id]==1, f1[obs_id], f0[obs_id])))
    if (length(mis_id)) ll_gmm <- ll_gmm + sum(log(den[mis_id]))
    q <- logistic_inv(alpha0 + alpha1*delta2)
    ll_miss <- sum((1-m)*log(clip01(1-q)) + m*log(clip01(q)))
    ll <- ll_gmm + ll_miss
    ll_trace[t] <- ll
    if (verbose && t %% 5 == 0) message(sprintf("LOGIT iter %3d | ll=%.4f | pi=%.3f", t, ll, pi1))
    if (t>1 && abs(ll_trace[t] - ll_trace[t-1]) < tol) { ll_trace <- ll_trace[1:t]; break }
  }
  list(pi1=pi1, mu1=mu1, S1=S1, mu0=mu0, S0=S0,
       alpha=c(alpha0=alpha0, alpha1=alpha1),
       tau = cbind(
         safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8) /
           pmax(safe_dmvnorm(X, mu1, S1) * clip01(pi1,1e-8) +
                  safe_dmvnorm(X, mu0, S0) * clip01(1-pi1,1e-8), 1e-300),
         NA),
       loglik=ll_trace)
}

## ========= 7) 训练（快速版 ECM + Logistic 基线） =========
fit_AO <- ECM_GMM_MAR_AO_fast(X, z_obs, lambda=0.7, estimate_lambda=FALSE,
                              max_iter=150, tol=1e-6, ridge=1e-4, verbose=FALSE)
p_AO <- fit_AO$tau[,1]

fit_LG <- ECM_GMM_MAR_LOGIT_fast(X, z_obs, max_iter=150, tol=1e-6, ridge=1e-4, verbose=FALSE)
p_LG <- fit_LG$tau[,1]

obs_id <- which(!is.na(z_obs))
df_train <- data.frame(y = as.integer(z_true[obs_id]==1), X[obs_id, , drop=FALSE])
glm_fit <- glm(y ~ ., data=df_train, family=binomial)
p_BL  <- predict(glm_fit, newdata=data.frame(X), type="response")

## ========= 8) 评估与可视化 =========
metrics_at_thr <- function(p_hat, y, thr=0.5) {
  y_pred <- as.integer(p_hat >= thr)
  TP <- sum(y_pred==1 & y==1); TN <- sum(y_pred==0 & y==0)
  FP <- sum(y_pred==1 & y==0); FN <- sum(y_pred==0 & y==1)
  ACC <- (TP+TN)/length(y)
  Precision <- ifelse(TP+FP>0, TP/(TP+FP), NA)
  Recall    <- ifelse(TP+FN>0, TP/(TP+FN), NA)
  F1 <- ifelse(is.na(Precision)|is.na(Recall)|(Precision+Recall)==0,
               NA, 2*Precision*Recall/(Precision+Recall))
  CM <- matrix(c(TN, FP, FN, TP), 2, 2, byrow=TRUE,
               dimnames=list("True"=c("Class0","Class1"),
                             "Pred"=c("Class0","Class1")))
  list(ACC=ACC, Precision=Precision, Recall=Recall, F1=F1, CM=CM)
}

eps <- 1e-12
roc_AO <- roc(y_true01, p_AO, quiet=TRUE)
roc_LG <- roc(y_true01, p_LG, quiet=TRUE)
roc_BL <- roc(y_true01, p_BL, quiet=TRUE)

auc_AO <- as.numeric(auc(roc_AO))
auc_LG <- as.numeric(auc(roc_LG))
auc_BL <- as.numeric(auc(roc_BL))

ll_AO <- -mean(y_true01*log(clip01(p_AO,eps)) + (1-y_true01)*log(clip01(1-p_AO,eps)))
ll_LG <- -mean(y_true01*log(clip01(p_LG,eps)) + (1-y_true01)*log(clip01(1-p_LG,eps)))
ll_BL <- -mean(y_true01*log(clip01(p_BL,eps)) + (1-y_true01)*log(clip01(1-p_BL,eps)))

br_AO <- mean((p_AO - y_true01)^2)
br_LG <- mean((p_LG - y_true01)^2)
br_BL <- mean((p_BL - y_true01)^2)

thr_AO  <- as.numeric(coords(roc_AO, "best", best.method="youden", ret="threshold"))
thr_LG  <- as.numeric(coords(roc_LG, "best", best.method="youden", ret="threshold"))
thr_BL  <- as.numeric(coords(roc_BL, "best", best.method="youden", ret="threshold"))

m_AO_fix <- metrics_at_thr(p_AO, y_true01, 0.5)
m_LG_fix <- metrics_at_thr(p_LG, y_true01, 0.5)
m_BL_fix <- metrics_at_thr(p_BL, y_true01, 0.5)

m_AO_opt <- metrics_at_thr(p_AO, y_true01, thr_AO)
m_LG_opt <- metrics_at_thr(p_LG, y_true01, thr_LG)
m_BL_opt <- metrics_at_thr(p_BL, y_true01, thr_BL)

cat("\n===== 混淆矩阵 (阈值=0.5) =====\n")
cat("\nECM-AO (0.5)\n");  print(m_AO_fix$CM)
cat("\nECM-Logit (0.5)\n");  print(m_LG_fix$CM)
cat("\nLogistic Baseline (0.5)\n");  print(m_BL_fix$CM)

cat("\n===== 混淆矩阵 (Youden 最优阈值) =====\n")
cat(sprintf("\nECM-AO (thr=%.3f)\n", thr_AO));  print(m_AO_opt$CM)
cat(sprintf("\nECM-Logit (thr=%.3f)\n", thr_LG));  print(m_LG_opt$CM)
cat(sprintf("\nLogistic Baseline (thr=%.3f)\n", thr_BL));  print(m_BL_opt$CM)

make_row <- function(name, auc, ll, br, m_fix, m_opt) {
  data.frame(
    Method = name,
    AUC = auc, LogLoss = ll, Brier = br,
    ACC_0.5 = m_fix$ACC, Precision_0.5 = m_fix$Precision, Recall_0.5 = m_fix$Recall, F1_0.5 = m_fix$F1,
    ACC_opt = m_opt$ACC, Precision_opt = m_opt$Precision, Recall_opt = m_opt$Recall, F1_opt = m_opt$F1
  )
}
result_table <- bind_rows(
  make_row("ECM (AO link)", auc_AO, ll_AO, br_AO, m_AO_fix, m_AO_opt),
  make_row("ECM (Logit link)", auc_LG, ll_LG, br_LG, m_LG_fix, m_LG_opt),
  make_row("Logistic (baseline)", auc_BL, ll_BL, br_BL, m_BL_fix, m_BL_opt)
)

# ROC 合图
# ROC 合图（统一配色：BL蓝 / AO橙 / Logit紫；全部实线）
lab_AO <- sprintf("ECM-AO (AUC=%.3f)", auc_AO)
lab_LG <- sprintf("ECM-Logit (AUC=%.3f)", auc_LG)
lab_BL <- sprintf("Logistic (AUC=%.3f)", auc_BL)

df_roc <- bind_rows(
  data.frame(fpr = 1-roc_AO$specificities, tpr = roc_AO$sensitivities, Method = lab_AO),
  data.frame(fpr = 1-roc_LG$specificities, tpr = roc_LG$sensitivities, Method = lab_LG),
  data.frame(fpr = 1-roc_BL$specificities, tpr = roc_BL$sensitivities, Method = lab_BL)
)

# 让图例顺序固定：baseline 在前
df_roc$Method <- factor(df_roc$Method, levels = c(lab_BL, lab_AO, lab_LG))

roc_colors <- setNames(
  c("#1F77B4", "#FF7F0E", "#9467BD"),
  c(lab_BL, lab_AO, lab_LG)
)

print(
  ggplot(df_roc, aes(fpr, tpr, color = Method)) +
    geom_path(size = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    coord_equal() +
    labs(
      title = "ROC: ECM-AO vs ECM-Logit vs Logistic (magic, MAR=70%)",
      x = "1 - Specificity", y = "Sensitivity"
    ) +
    scale_color_manual(values = roc_colors) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank())
)

## ===== 8.1) 阈值扫描：ECM-AO vs Logistic（默认 0.3–0.7） =====
# 阈值网格：默认 0.3–0.7，步长 0.05
thr_grid <- seq(0.3, 0.7, by = 0.05)

get_metrics_curve <- function(p_hat, y, label, thr_seq = thr_grid) {
  rows <- lapply(thr_seq, function(t) {
    m <- metrics_at_thr(p_hat, y, thr = t)
    c(threshold = t,
      Precision = m$Precision,
      Recall    = m$Recall,
      Accuracy  = m$ACC)
  })
  df <- as.data.frame(do.call(rbind, rows))
  df$Method <- label
  num_cols <- c("threshold","Precision","Recall","Accuracy")
  df[num_cols] <- lapply(df[num_cols], function(x) round(as.numeric(x), 4))
  df
}

curve_AO <- get_metrics_curve(p_AO, y_true01, "ECM-AO")
curve_BL <- get_metrics_curve(p_BL, y_true01, "Logistic")
# 如需加入 ECM-Logit：
# curve_LG <- get_metrics_curve(p_LG, y_true01, "ECM-Logit")

metrics_curve <- dplyr::bind_rows(curve_AO, curve_BL)  # 可再加 curve_LG
cat("\n===== Precision / Recall / Accuracy vs Threshold (0.3–0.7) =====\n")
print(metrics_curve)

write.csv(metrics_curve, "metrics_threshold_curve_AO_vs_Logistic.csv", row.names = FALSE)
cat('\n已保存到: "metrics_threshold_curve_AO_vs_Logistic.csv"\n')

metrics_long <- metrics_curve |>
  tidyr::pivot_longer(cols = c(Precision, Recall, Accuracy),
                      names_to = "Metric", values_to = "Value")


thr_colors <- c("Logistic" = "#1F77B4", "ECM-AO" = "#FF7F0E")

print(
  ggplot(metrics_long, aes(x = threshold, y = Value, color = Method)) +
    geom_line(size = 1.1) +
    geom_point(size = 1.8) +
    facet_wrap(~ Metric, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = thr_grid) +
    scale_color_manual(values = thr_colors) +
    labs(
      title = "ECM-AO vs Logistic (threshold 0.3–0.7)",
      x = "Classification threshold", y = "Metric value"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom", legend.title = element_blank())
)

## ========= 9) 合成实验：AUC vs Missing Rate（GMM under MAR） =========
set.seed(2025)

## ===== 9.1) GMM 数据生成 =====
n <- 500
pi <- c(0.5, 0.5)
mu <- list(c(-1,0), c(1,0))
Sigma <- list(diag(2), diag(2))

posterior_prob <- function(y, pi, mu, Sigma){
  f1 <- dmvnorm(y, mean=mu[[1]], sigma=Sigma[[1]]) * pi[1]
  f2 <- dmvnorm(y, mean=mu[[2]], sigma=Sigma[[2]]) * pi[2]
  f1 / (f1 + f2)
}
margin_conf <- function(y){ abs(2*posterior_prob(y,pi,mu,Sigma) - 1) }

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
data_full <- sample_data(n, pi, mu, Sigma)
data_full$mc <- apply(data_full[,c("x1","x2")],1,margin_conf)

## ===== 9.2) AO 缺失机制 =====
q_aranda_ordaz <- function(delta2, alpha0, alpha1, lambda){
  eta <- alpha0 + alpha1 * delta2
  q <- 1 - (1 + lambda * exp(eta))^(-1/lambda)
  pmin(pmax(q, 1e-8), 1-1e-8)
}
calibrate_alpha0 <- function(delta2, alpha1, lambda, target){
  f <- function(a0) mean(q_aranda_ordaz(delta2, a0, alpha1, lambda)) - target
  uniroot(f, c(-50,50))$root
}

## ===== 9.3) 主循环：不同缺失比例 =====
miss_levels <- seq(0.5, 0.95, by=0.05)
auc_results <- data.frame(MissingRate = miss_levels,
                          AUC_ECM_AO = NA, AUC_Logistic = NA)

lambda <- 0.5; alpha1 <- -6

for (i in seq_along(miss_levels)) {
  target_missing <- miss_levels[i]
  
  delta2 <- data_full$mc^2
  alpha0 <- calibrate_alpha0(delta2, alpha1, lambda, target_missing)
  q_miss <- q_aranda_ordaz(delta2, alpha0, alpha1, lambda)
  
  dat2 <- data_full
  dat2$label_true <- dat2$label
  dat2$label_missing <- rbinom(nrow(dat2), 1, q_miss)
  dat2$label_obs <- dat2$label_true
  dat2$label_obs[dat2$label_missing == 1] <- NA_integer_
  
  ## === ECM-AO ===
  fit_AO2 <- ECM_GMM_MAR_AO_fast(
    as.matrix(dat2[,c("x1","x2")]),
    dat2$label_obs,
    lambda=lambda,
    estimate_lambda=FALSE,
    max_iter=120,
    tol=1e-5
  )
  p_AO2 <- fit_AO2$tau[,1]
  auc_AO2 <- as.numeric(auc(dat2$label_true==1, p_AO2))
  
  ## === Logistic baseline ===
  obs_id2 <- which(!is.na(dat2$label_obs))
  df_train2 <- data.frame(
    y = as.integer(dat2$label_true[obs_id2]==1),
    dat2[obs_id2,c("x1","x2")]
  )
  glm_fit2 <- glm(y ~ ., data=df_train2, family=binomial)
  p_logit2 <- predict(glm_fit2, newdata=dat2, type="response")
  auc_logit2 <- as.numeric(auc(dat2$label_true==1, p_logit2))
  
  auc_results$AUC_ECM_AO[i] <- auc_AO2
  auc_results$AUC_Logistic[i] <- auc_logit2
  
  cat(sprintf("缺失率 %.0f%%: AUC_ECM_AO=%.3f, AUC_Logistic=%.3f\n",
              100*target_missing, auc_AO2, auc_logit2))
}

## ===== 9.4) 绘图：AUC vs Missing Rate =====
df_plot <- auc_results %>%
  pivot_longer(cols=c(AUC_ECM_AO, AUC_Logistic),
               names_to="Method", values_to="AUC") %>%
  mutate(Method = factor(Method,
                         levels=c("AUC_ECM_AO","AUC_Logistic"),
                         labels=c("ECM-AO","Logistic")))

print(
  ggplot(df_plot, aes(x=MissingRate, y=AUC, color=Method)) +
    geom_line(size=1.2) + geom_point(size=2) +
    scale_x_continuous(labels=scales::percent_format(accuracy=1)) +
    labs(title="AUC vs Missing Rate (GMM under MAR)",
         x="Missing proportion", y="AUC", color="Method") +
    theme_minimal(base_size=13) +
    theme(legend.position="bottom", legend.title = element_blank())
)
