library(tidyverse)
library(moments)
library(alabama)
library(lmom)
library(robustbase)
library(VGAM)
library(evd)
library(quantreg)


# ================================
# 0. Load and preprocess data
# ================================

file_path <- "C:\\Data\\snp500_2008_2024.csv"
df_org <- read_csv(file_path)
df_org <- df_org %>% rename(Date = `...1`)


# Extract date and closing price
df <- df_org %>% select(Date, Close)

df <- df %>%
  mutate(Return = c(NA, diff(Close) / head(Close, -1))) %>%
  drop_na()

ggplot(df, aes(x = Date, y = Close)) +
  geom_line() +
  scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
  labs(title = "Close Price Over Time",
       x = "Date",
       y = "Close") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Load functions
# Define the objective function for optimization
objective_function <- function(beta, X_ordered, z_p) {
  beta_0 <- beta[1]
  beta_1 <- beta[2]
  beta_2 <- beta[3]
  beta_3 <- beta[4]
  residuals <- X_ordered - (beta_0 + beta_1 * z_p + beta_2 * z_p^2 + beta_3 * z_p^3)
  return(sum(abs(residuals)))
}

# Calculate residuals
calculate_residuals <- function(X_ordered, z_p, optimized_b) {
  b_0 <- optimized_b[1]
  b_1 <- optimized_b[2]
  b_2 <- optimized_b[3]
  b_3 <- optimized_b[4]
  
  residuals <- X_ordered - (b_0 + b_1 * z_p + b_2 * z_p^2 + b_3 * z_p^3)
  return(residuals)
}


# Function to calculate initial beta values based on skewness and kurtosis
# Only used New1 
calculate_initial_beta <- function(X_ordered, z_p, method = "New1") {
  gamma <- skewness(X_ordered)
  kappa <- kurtosis(X_ordered)
  if (method == "New1") {
    initial_beta <- c(
      -gamma / 6,
      1 - (1 / 8) * kappa + (5 / 36) * gamma^2,
      gamma / 6,
      (kappa / 24) - (gamma^2 / 18)
    )
  } else if (method == "New2") {
    initial_beta <- c(
      -gamma / 6,
      1 - (1 / 8) * kappa,
      gamma / 6,
      (kappa / 24)
    )
  }
  return(as.vector(initial_beta))
}


# Function to perform optimization
perform_optimization_new1 <- function(initial_beta, X_ordered, z_p) {
  result <- tryCatch({
    optim(
      par = initial_beta,
      fn = function(beta) objective_function(beta, X_ordered, z_p),
      method = "BFGS",
      control = list(maxit = 1000, trace = FALSE)
    )
  }, error = function(e) {
    return(NULL)
  })
  return(result)
}


# Function to estimate alpha_3(skewness), alpha_4(kurtosis) from optimized beta values
calculate_alphas <- function(optimized_b) {
  b_0 <- optimized_b[1]
  b_1 <- optimized_b[2]
  b_2 <- optimized_b[3]
  b_3 <- optimized_b[4]
  moment_2 <- b_1^2+2*b_2^2+15*b_3^2+6*b_1*b_3
  moment_3 <- 6 * b_1^2 * b_2 + 8 * b_2^3 + 72 * b_1 * b_2 * b_3 + 270 * b_2 * b_3^2 #3rd moment
  moment_4 <- 3 * (b_1^4 + 20 * b_1^3 * b_3 + 210 * b_1^2 * b_3^2 + 1260 * b_1 * b_3^3 + 3465 * b_3^4) +
    12 * b_2^2 * (5 * b_1^2 + 5 * b_2^2 + 78 * b_1 * b_3 + 375 * b_3^2) #4th moment
  return(list(alpha_3 = moment_3/((moment_2)^(3/2)), alpha_4 = moment_4/(moment_2^2))) # Return alpha values as a list
}

### Section 2. Models of Existing Methods ###
# Function to calculate L-moments
calculate_l_moments <- function(X) {
  lmr <- samlmu(X)
  tau_3 <- lmr[3]
  tau_4 <- lmr[4]
  return(list(alpha_3 = tau_3, alpha_4 = tau_4))
}

# Function to calculate trimmed L-moments
calculate_tl_moments <- function(X, trim = 1) {
  lmr <- samlmu(X, trim = trim)
  tau_3 <- lmr[3]
  tau_4 <- lmr[4]
  return(list(alpha_3 = tau_3, alpha_4 = tau_4))
}

# Function to calculate skewness and kurtosis using quantiles (Method 1 - Bowley skewness, Moors kurtosis)
calculate_quantile1_skewness_kurtosis <- function(X) {
  Q1 <- quantile(X, probs = 0.25)
  Q2 <- quantile(X, probs = 0.5)
  Q3 <- quantile(X, probs = 0.75)
  sk1 <- (Q3 + Q1 - 2 * Q2) / (Q3 - Q1)
  kr1 <- ((quantile(X, probs = 0.875) - quantile(X, probs = 0.625)) +
            (quantile(X, probs = 0.375) - quantile(X, probs = 0.125))) /
    (quantile(X, probs = 0.75) - quantile(X, probs = 0.25))
  return(list(alpha_3 = sk1, alpha_4 = kr1))
}

# Function to calculate skewness and kurtosis using quantiles (Method 2 - Groeneveld skewness, Hogg coefficient kurtosis)
calculate_quantile2_skewness_kurtosis <- function(X) {
  Q2 <- quantile(X, probs = 0.5) # Median
  mu <- mean(X)                 
  abs_dev <- mean(abs(X - Q2))
  
  sk2 <- (mu - Q2) / abs_dev
  
  calculate_U <- function(alpha) {
    y_values <- seq(1 - alpha, 1, length.out = 1000) 
    mean(sapply(y_values, function(y) quantile(X, probs = y))) 
  }
  calculate_L <- function(alpha) {
    y_values <- seq(0, alpha, length.out = 1000) 
    mean(sapply(y_values, function(y) quantile(X, probs = y)))
  }
  
  U_0.05 <- calculate_U(0.05)
  L_0.05 <- calculate_L(0.05)
  U_0.5 <- calculate_U(0.5)
  L_0.5 <- calculate_L(0.5)
  
  kr2 <- (U_0.05 - L_0.05) / (U_0.5 - L_0.5) - 2.59
  return(list(alpha_3 = sk2, alpha_4 = kr2))
}


# Function to calculate skewness and kurtosis using quantiles (Method 3 - Pearson skewness, centered Crow kurtosis)
calculate_quantile3_skewness_kurtosis <- function(X) {
  sk3 <- (mean(X) - quantile(X, probs = 0.5)) / sd(X)
  kr3 <- (quantile(X, probs = 0.975) + quantile(X, probs = 0.025)) /
    (quantile(X, probs = 0.75) - quantile(X, probs = 0.25)) - 2.91
  return(list(alpha_3 = sk3, alpha_4 = kr3))
}


# Function to calculate medcouple (MC) skewness
calculate_mc_skewness <- function(X) {
  sk4 <- mc(X)
  return(list(alpha_3 = sk4, alpha_4 = 0))
}

# Function to calculate modal skewness
calculate_modal_skewness <- function(X) {
  sorted_X <- sort(X)
  n <- length(X)
  half_range_size <- ceiling(n / 2) # Define half-range size
  min_range_width <- Inf # Initialize minimum range width
  mode_value <- NA # Initialize mode value
  
  # Find the mode by minimizing the range width in the half-range
  for (i in 1:(n - half_range_size + 1)) {
    range_width <- sorted_X[i + half_range_size - 1] - sorted_X[i]
    if (range_width < min_range_width) {
      min_range_width <- range_width
      mode_value <- mean(sorted_X[i:(i + half_range_size - 1)])
    }
  }
  count_less_than_mode <- sum(X < mode_value) # Count values below the mode
  sk5 <- 1 - (2 * count_less_than_mode / n)
  return(list(alpha_3 = sk5, alpha_4 = 0))
}

calculate_kde_moments <- function(X, method = "silverman") {
  n <- length(X)
  s2 <- mean((X - mean(X))^2)
  s3 <- mean((X - mean(X))^3)
  s4 <- mean((X - mean(X))^4)
  
  sd_X <- sd(X)
  iqr_X <- IQR(X)
  
  if (method == "silverman") {
    h <- 0.9 * min(sd_X, iqr_X / 1.34) * n^(-1/5)
  } else if (method == "scott") {
    h <- 1.06 * sd_X * n^(-1/5)
  } else {
    stop("Method must be either 'silverman' or 'scott'")
  }
  
  sigma_kde2 <- s2 + h^2
  skew_kde <- s3 / (sigma_kde2)^(3/2)
  kurt_kde <- (s4 + 6 * h^2 * s2 + 3 * h^4) / (sigma_kde2^2)
  
  return(list(alpha_3 = skew_kde, alpha_4 = kurt_kde))
}


# Function to run the model with optional reverse ordering of z_p
run_model <- function(X, method = "New1", reverse_z_p = FALSE) {
  mu_hat <- mean(X)
  sigma_hat <- sd(X)
  X.centered <- X - mu_hat
  X.std <- X.centered / sigma_hat
  X_ordered <- sort(X.std)
  n <- length(X.std)
  
  quantiles <- seq(0, 1, length.out = n + 2)[-c(1, n + 2)]
  z_p <- qnorm(quantiles, 0, 1)
  
  if (reverse_z_p) {
    z_p <- rev(z_p)
  }
  
  initial_beta <- calculate_initial_beta(X_ordered, z_p, method)
  result <- perform_optimization_new1(initial_beta, X_ordered, z_p)
  
  if (!is.null(result)) {
    optimized_b <- result$par
    alphas <- calculate_alphas(optimized_b)
    failure <- 0
  } else {
    optimized_b <- NA
    alphas <- list(alpha_3 = NA, alpha_4 = NA)
    failure <- 1
  }
  return(list(alphas = alphas, optimized_b = optimized_b, failure = failure))
}

# Function to run the smallresid model (choosing beta with smaller residual sum)
run_new5_model <- function(X) {
  pos_model <- run_model(X, "New1")
  neg_model <- run_model(X, "New1", reverse_z_p = TRUE)
  
  if (pos_model$failure == 1 || neg_model$failure == 1) {
    return(list(alphas = list(alpha_3 = NA, alpha_4 = NA), failure = 1))
  }
  
  X_ordered <- sort((X - mean(X)) / sd(X))
  n <- length(X)
  z_p <- qnorm(seq(1, n) / (n + 1), 0, 1)
  z_p_reversed <- rev(z_p)
  
  resid_pos <- sum(abs(calculate_residuals(X_ordered, z_p, pos_model$optimized_b)))
  resid_neg <- sum(abs(calculate_residuals(X_ordered, z_p_reversed, neg_model$optimized_b)))
  
  if (resid_pos < resid_neg) {
    best_b <- pos_model$optimized_b
  } else {
    best_b <- neg_model$optimized_b
  }
  
  alphas_new5 <- calculate_alphas(best_b)
  
  return(list(alphas = alphas_new5, failure = 0))
}



# ================================
# 1. Sliding window analysis
# ================================

sliding_window_analysis <- function(df, window_size) {
  n <- nrow(df)
  result <- tibble()
  
  for (i in 1:(n - window_size + 1)) {
    window_data <- df$Return[i:(i + window_size - 1)]
    current_date <- df$Date[i + window_size - 1]
    
    row <- tibble(Date = current_date)
    
    # --- Empirical moments ---
    row$Empirical_skew <- skewness(window_data)
    row$Empirical_kurt <- kurtosis(window_data)
    
    # --- New1 ---
    new1 <- run_model(window_data, method = "New1")
    row$New1_skew <- new1$alphas$alpha_3
    row$New1_kurt <- new1$alphas$alpha_4
    
    # --- New5 (small residual) ---
    new5 <- run_new5_model(window_data)
    row$New5_smallresid_skew <- new5$alphas$alpha_3
    row$New5_smallresid_kurt <- new5$alphas$alpha_4
    
    # --- L-moments ---
    lmom <- calculate_l_moments(window_data)
    row$lmom_skew <- lmom$alpha_3
    row$lmom_kurt <- lmom$alpha_4
    
    # --- TL (trimmed L-moments) ---
    tl <- calculate_tl_moments(window_data)
    row$TL_skew <- tl$alpha_3
    row$TL_kurt <- tl$alpha_4
    
    # --- Quantile 1 ~ 3 ---
    q1 <- calculate_quantile1_skewness_kurtosis(window_data)
    row$Quantile1_skew <- q1$alpha_3
    row$Quantile1_kurt <- q1$alpha_4
    
    q2 <- calculate_quantile2_skewness_kurtosis(window_data)
    row$Quantile2_skew <- q2$alpha_3
    row$Quantile2_kurt <- q2$alpha_4
    
    q3 <- calculate_quantile3_skewness_kurtosis(window_data)
    row$Quantile3_skew <- q3$alpha_3
    row$Quantile3_kurt <- q3$alpha_4
    
    # --- MC method ---
    mc <- calculate_mc_skewness(window_data)
    row$mc_skew <- mc$alpha_3
    
    # --- Modal method ---
    modal <- calculate_modal_skewness(window_data)
    row$modal_skew <- modal$alpha_3
    
    # --- KDE methods ---
    kde_sm <- calculate_kde_moments(window_data,method = "silverman")
    row$silverman_skew<-kde_sm$alpha_3
    row$silverman_kurt<-kde_sm$alpha_4
    
    kde_scott <- calculate_kde_moments(window_data,method = "scott")
    row$scott_skew<-kde_scott$alpha_3
    row$scott_kurt<-kde_scott$alpha_4

    result <- bind_rows(result, row)
  }
  
  return(result)
}

# ================================
# 2. Run analysis (100-day windows)
# ================================
result_100 <- sliding_window_analysis(df, window_size = 100)
