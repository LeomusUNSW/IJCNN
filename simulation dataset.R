set.seed(2025)
library(mvtnorm)

# ---- Parameters ----
n <- 500   # sample size
d <- 2     # dimension
pi <- c(0.5, 0.5)  # mixing proportions
mu <- list(c(-1,0), c(1,0))  # means of two components
Sigma <- list(diag(2), diag(2)) # covariance matrices

# ---- Posterior probability function ----
posterior_prob <- function(y, pi, mu, Sigma){
  f1 <- dmvnorm(y, mean=mu[[1]], sigma=Sigma[[1]]) * pi[1]
  f2 <- dmvnorm(y, mean=mu[[2]], sigma=Sigma[[2]]) * pi[2]
  tau1 <- f1 / (f1 + f2)
  return(tau1)
}

# ---- Margin confidence function ----
margin_conf <- function(y) {
  tau1 <- posterior_prob(y, pi, mu, Sigma)
  abs(2*tau1 - 1)
}

# ---- Generate data with margin confidence in [0,0.6] ----
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
  return(as.data.frame(out))
}

data <- sample_data(n, pi, mu, Sigma)

# ---- Add margin confidence ----
data$mc <- apply(data[, c("x1","x2")], 1, margin_conf)

# ---- Aranda–Ordaz missingness ----
q_aranda_ordaz <- function(delta2, alpha0, alpha1, lambda){
  eta <- alpha0 + alpha1 * delta2
  q   <- 1 - (1 + lambda * exp(eta))^(-1/lambda)
  pmin(pmax(q, 1e-8), 1-1e-8)
}

calibrate_alpha0 <- function(delta2, alpha1, lambda, target){
  f <- function(a0){
    mean(q_aranda_ordaz(delta2, a0, alpha1, lambda)) - target
  }
  uniroot(f, c(-50, 50))$root
}

delta2  <- data$mc^2
lambda  <- 0.5
alpha1  <- -6
target_missing <- 0.7

alpha0 <- calibrate_alpha0(delta2, alpha1, lambda, target_missing)
q_miss <- q_aranda_ordaz(delta2, alpha0, alpha1, lambda)

# ---- Store true labels & apply missingness ----
data$label_true    <- data$label                  
data$label_missing <- rbinom(nrow(data), 1, q_miss) 
data$q_miss        <- q_miss                        
data$label[data$label_missing == 1] <- NA_integer_ 

# ---- Checks ----
mean(is.na(data$label))   # realized missing proportion
tapply(q_miss, cut(data$mc, breaks=c(0,0.2,0.4,0.6), include.lowest=TRUE), mean)

# ---- Packages ----
library(ggplot2)
library(cowplot)
library(scales)

# ---- Palette & base theme ----
pal <- c("Component 1" = "#2F5597",
         "Component 2" = "#E15759",
         "Missing"     = "#9E9E9E")

base_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3, colour = "#E6E6E6")
  )

# ---- Prepare 3-way label for plotting ----
plot_df <- data
plot_df$label3 <- ifelse(is.na(plot_df$label), "Missing",
                         ifelse(plot_df$label == 1, "Component 1", "Component 2"))
plot_df$label3 <- factor(plot_df$label3,
                         levels = c("Component 1","Component 2","Missing"))

# ---- Scatter (left) ----
plot_scatter <- ggplot(plot_df, aes(x = x1, y = x2, fill = label3)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.85,
             colour = "white", stroke = 0.2) +
  scale_fill_manual(values = pal, name = "Label") +
  labs(title = "GMM Sample under MAR", x = "x1", y = "x2") +
  base_theme +
  theme(legend.position = "none")

# ---- Pie (right, percentages centered) ----
pie_tab <- as.data.frame(table(plot_df$label3), stringsAsFactors = FALSE)
names(pie_tab) <- c("category", "n")
pie_tab$prop <- pie_tab$n / sum(pie_tab$n)

plot_pie <- ggplot(pie_tab, aes(x = "", y = prop, fill = category)) +
  geom_col(width = 1, colour = "white", linewidth = 0.6) +
  coord_polar(theta = "y", start = pi/2, direction = -1) +
  geom_text(aes(label = percent(prop, accuracy = 0.1)),
            position = position_stack(vjust = 0.5),
            colour = "white", size = 4.4, fontface = "bold") +
  scale_fill_manual(values = pal, name = "Label",
                    breaks = c("Component 1","Component 2","Missing")) +
  labs(title = "Label Composition under MAR") +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    legend.position = "none"
  )

# ---- Shared legend (bottom) ----
legend_shared <- get_legend(
  ggplot(plot_df, aes(x = x1, y = x2, fill = label3)) +
    geom_point(shape = 21) +
    scale_fill_manual(values = pal, name = "Label") +
    theme_minimal() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text  = element_text(size = 10))
)

# ---- Combine: scatter | pie + legend ----
row_main   <- plot_grid(plot_scatter, plot_pie,
                        ncol = 2, rel_widths = c(1.35, 1))
final_plot <- plot_grid(row_main, legend_shared,
                        ncol = 1, rel_heights = c(1, 0.12))

# ---- Show ----
print(final_plot)

