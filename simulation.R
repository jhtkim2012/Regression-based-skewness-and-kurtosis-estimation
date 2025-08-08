
library(moments)
library(alabama)
library(lmom)
library(robustbase)
library(VGAM)
library(evd)
library(quantreg)


### Section 1. Functions for the New Method ###
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
    (quantile(X, probs = 0.75) - quantile(X, probs = 0.25))
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


### Section 2. Models of New Methods ###

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


### Section 3. Generate Mixture Distribution Data and Compute Moments ###
generate_mixture_data <- function(n, weight1, dist1 = "gamma", dist2 = "gamma", 
                                  shape1 = 1, scale1 = 1, shape2 = 3, scale2 = 1, 
                                  mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1, df1 = NULL, df2 = NULL, xi1 = NULL, xi2 = NULL) {
  n1 <- round(n * weight1)  # Number of samples from the first distribution
  n2 <- n - n1  # Number of samples from the second distribution
  
  generate_data <- function(n, dist, shape, scale, mu, sigma, df, xi) {
    if (dist == "gamma") {
      return(rgamma(n, shape = shape, scale = scale))
    } else if (dist == "lognormal") {
      return(rlnorm(n, meanlog = mu, sdlog = sigma))
    } else if (dist == "t") {
      if (is.null(df)) stop("Degrees of freedom (df) must be specified for t distribution.")
      
      scale_t <- sigma / sqrt(df / (df - 2))
      
      return(rt(n, df = df) * scale_t + mu)
    } else if (dist == "gpd") {
      return(evd::rgpd(n, loc = 0, scale = scale, shape = xi))
    } else {
      stop("Unsupported distribution type")
    }
  }
  
  X1 <- generate_data(n1, dist1, shape1, scale1, mu1, sigma1, df1, xi1)
  X2 <- generate_data(n2, dist2, shape2, scale2, mu2, sigma2, df2, xi2)
  
  X_mixture <- c(X1, X2)
  return(list(X_mixture = X_mixture))
}

# Function to compute moments based on the given distribution parameters
calculate_moments <- function(dist, shape, scale, mu, sigma, df, xi) {
  if (dist == "gamma") {
    E <- shape * scale
    Var <- shape * scale^2
    Skew <- 2 / sqrt(shape)
    Kurt <- 6 / shape + 3
  } else if (dist == "lognormal") {
    E <- exp(mu + (sigma^2) / 2)
    Var <- (exp(sigma^2) - 1) * exp(2 * mu + sigma^2)
    Skew <- (exp(sigma^2) + 2) * sqrt(exp(sigma^2) - 1)
    Kurt <- exp(4 * sigma^2) + 2 * exp(3 * sigma^2) + 3 * exp(2 * sigma^2) - 3
  } else if (dist == "t") {
    
    scale_t <- sigma / sqrt(df / (df - 2))  
    
    E <- ifelse(df > 1, mu, NA)
    Var <- ifelse(df > 2, scale_t^2 * df / (df - 2), NA)
    Skew <- 0
    Kurt <- ifelse(df > 4, 6 / (df - 4) + 3, Inf)
  } else if (dist == "gpd") {
    E <- ifelse(xi < 1, scale / (1 - xi), Inf)
    Var <- ifelse(xi < 0.5, (scale^2) / ((1 - xi)^2 * (1 - 2 * xi)), Inf)
    Skew <- ifelse(xi < 1/3, 2 * (1 + xi) * sqrt(1 - 2 * xi) / (1 - 3 * xi), Inf)
    Kurt <- ifelse(xi < 1/4, 3 * (1 - 2 * xi) * (2 * xi^2 + xi + 3) / ((1 - 3 * xi) * (1 - 4 * xi)), Inf)
  } else {
    stop("Unsupported distribution type")
  }
  return(list(E = E, Var = Var, Skew = Skew, Kurt = Kurt))
}


# Function to compute moments for a mixture distribution
calculate_mixture_moments <- function(weight1, dist1 = "gamma", dist2 = "gamma", 
                                      shape1 = 1, scale1 = 1, shape2 = 3, scale2 = 1, 
                                      mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1, df1 = NULL, df2 = NULL, xi1 = NULL, xi2 = NULL) {
  
  
  moments1 <- calculate_moments(dist1, shape1, scale1, mu1, sigma1, df1, xi1)
  moments2 <- calculate_moments(dist2, shape2, scale2, mu2, sigma2, df2, xi2)
  
  E_mix <- weight1 * moments1$E + (1 - weight1) * moments2$E
  Var_mix <- weight1 * (moments1$Var + moments1$E^2) + (1 - weight1) * (moments2$Var + moments2$E^2) - E_mix^2
  
  E3_mix <- weight1 * ((moments1$E - E_mix)^3 + 3 * (moments1$E - E_mix) * moments1$Var + moments1$Skew * (moments1$Var^(3/2))) + 
    (1 - weight1) * ((moments2$E - E_mix)^3 + 3 * (moments2$E - E_mix) * moments2$Var + moments2$Skew * (moments2$Var^(3/2)))
  Skew_mix <- E3_mix / (Var_mix^(3/2))
  
  E4_mix <- weight1 * ((moments1$E - E_mix)^4 + 
                         6 * (moments1$E - E_mix)^2 * moments1$Var + 
                         4 * (moments1$E - E_mix) * moments1$Skew * (moments1$Var^(3/2)) + 
                         moments1$Kurt * moments1$Var^2) + 
    (1 - weight1) * ((moments2$E - E_mix)^4 + 
                       6 * (moments2$E - E_mix)^2 * moments2$Var + 
                       4 * (moments2$E - E_mix) * moments2$Skew * (moments2$Var^(3/2)) + 
                       moments2$Kurt * moments2$Var^2)
  
  Kurt_mix <- E4_mix / (Var_mix^2)
  
  return(list(E_mix = E_mix, Var_mix = Var_mix, Skew_mix = Skew_mix, Kurt_mix = Kurt_mix))
}


### Section 3. Simulation ###
simulation_results <- function(n, dist_type = "gamma", shape_gamma = 2, scale_gamma = 1, mu_lognormal = 1, sigma_lognormal = 0.5, df_t = NULL, loc_t = 0, scale_t = 1, sigma_gpd = 1, xi_gpd = 0.5) {
  nsim <- 1000 # Number of simulations
  
  # Compute theoretical skewness and kurtosis for different distributions
  if (dist_type == "gamma") {
    true_skewness <- 2 / sqrt(shape_gamma)
    true_kurtosis <- 6 / shape_gamma + 3
  } else if (dist_type == "lognormal") {
    true_skewness <- (exp(sigma_lognormal^2) + 2) * sqrt(exp(sigma_lognormal^2) - 1)
    true_kurtosis <- exp(4 * sigma_lognormal^2) + 2 * exp(3 * sigma_lognormal^2) + 3 * exp(2 * sigma_lognormal^2) - 3
  } else if (dist_type == "t") {
    if (is.null(df_t)) stop("Degrees of freedom (df_t) must be specified for t distribution.")
    
    true_skewness <- 0  
    if (df_t > 4) {
      true_kurtosis <- 6 / (df_t - 4) + 3
    } else {
      true_kurtosis <- Inf  # Kurtosis defined only if df > 4
    }
  } else if (dist_type == "gpd") {
    if (xi_gpd >= 0) {
      true_skewness <- 2 * (1 + xi_gpd) * sqrt(1 - 2 * xi_gpd) / (1 - 3 * xi_gpd)
      true_kurtosis <- 3*(1 - 2 * xi_gpd) * (2 * xi_gpd^2 + xi_gpd + 3) / ((1 - 3 * xi_gpd) * (1 - 4 * xi_gpd))
    } else {
      stop("Invalid shape parameter for GPD. Shape must be < 0.5 to define skewness and kurtosis.")
    }
  } else {
    stop("Unsupported distribution type. Use 'gamma', 'lognormal', 't', or 'gpd'.")
  }
  
  true_values <- c(true_skewness, true_kurtosis) # Store true skewness and kurtosis
  
  methods <- c("Empirical", "New1", "New5_smallresid", "lmom", "TL", "Quantile1",
               "Quantile2", "Quantile3", "mc", "modal", "KDE_Silverman", "KDE_Scott")
  stat_names <- c("Skewness", "Kurtosis")
  
  results <- list()
  estimates_storage <- list()
  
  # Create data frames to store bias, variance, and sqrt(MSE)
  for (stat in stat_names) {
    results[[stat]] <- data.frame(
      Method = methods,
      True_Value = rep(0, length(methods)),
      Bias = rep(0, length(methods)),
      Variance = rep(0, length(methods)),
      root_MSE = rep(0, length(methods)),
      Overestimate = rep(0, length(methods)),
      Underestimate = rep(0, length(methods)),
      number_of_failure = rep(0, length(methods))
    )
    results[[stat]]$True_Value <- ifelse(stat == "Skewness", true_skewness, true_kurtosis)
    estimates_storage[[stat]] <- matrix(NA, nrow = nsim, ncol = length(methods))
  }
  
  set.seed(1)
  isim <- 0
  total_attempts <- 0
  failures <- rep(0, length(methods))
  
  while (isim < nsim) {
    total_attempts <- total_attempts + 1
    
    # Generate random samples based on selected distribution
    if (dist_type == "gamma") {
      X <- rgamma(n, shape = shape_gamma, scale = scale_gamma)
    } else if (dist_type == "lognormal") {
      X <- rlnorm(n, meanlog = mu_lognormal, sdlog = sigma_lognormal)
    } else if (dist_type == "t") {
      X <- rt(n, df = df_t) * scale_t + loc_t
    } else if (dist_type == "gpd") {
      X <- evd::rgpd(n, loc = 0, scale = sigma_gpd, shape = xi_gpd)
    }
    
    # Apply different estimation methods
    methods_res <- list(
      Empirical = list(alphas = list(alpha_3 = skewness(X), alpha_4 = kurtosis(X)), failure = 0),
      New1 = run_model(X, "New1"),
      New5_smallresid = run_new5_model(X),
      lmom = list(alphas = calculate_l_moments(X), failure = 0),
      TL = list(alphas = calculate_tl_moments(X), failure = 0),
      Quantile1 = list(alphas = calculate_quantile1_skewness_kurtosis(X), failure = 0),
      Quantile2 = list(alphas = calculate_quantile2_skewness_kurtosis(X), failure = 0),
      Quantile3 = list(alphas = calculate_quantile3_skewness_kurtosis(X), failure = 0),
      mc = list(alphas = calculate_mc_skewness(X), failure = 0),
      modal = list(alphas = calculate_modal_skewness(X), failure = 0),
      KDE_Silverman = list(alphas = calculate_kde_moments(X, method = "silverman"), failure = 0),
      KDE_Scott = list(alphas = calculate_kde_moments(X, method = "scott"), failure = 0)
    )
    
    if (methods_res$New1$failure == 1) {
      failures[2] <- failures[2] + methods_res$New1$failure
      next
    }
    isim <- isim + 1
    #print(isim)
    
    # Store results for each method
    for (method_index in seq_along(methods)) {
      method <- methods[method_index]
      estimates <- methods_res[[method]]$alphas
      for (i in 1:2) {
        true_value <- true_values[i]
        estimate <- ifelse(i == 1, estimates$alpha_3, estimates$alpha_4)
        if (!is.na(estimate)) {
          results[[stat_names[i]]]$Bias[method_index] <- 
            results[[stat_names[i]]]$Bias[method_index] + (estimate - true_value)
          results[[stat_names[i]]]$root_MSE[method_index] <- 
            results[[stat_names[i]]]$root_MSE[method_index] + (estimate - true_value)^2
          estimates_storage[[stat_names[i]]][isim, method_index] <- estimate
          if (estimate > true_value) {
            results[[stat_names[i]]]$Overestimate[method_index] <- 
              results[[stat_names[i]]]$Overestimate[method_index] + 1
          } else {
            results[[stat_names[i]]]$Underestimate[method_index] <- 
              results[[stat_names[i]]]$Underestimate[method_index] + 1
          }
        }
      }
      failures[method_index] <- failures[method_index] + methods_res[[method]]$failure
    }
  }
  
  # Finalize bias, variance, and sqrt(MSE) calculations
  for (stat in stat_names) {
    results[[stat]]$Bias <- results[[stat]]$Bias / nsim
    results[[stat]]$Variance <- apply(estimates_storage[[stat]], 2, var, na.rm = TRUE)
    results[[stat]]$root_MSE <- sqrt(results[[stat]]$Bias^2 + results[[stat]]$Variance)
    results[[stat]]$number_of_failure <- unname(failures)
  }
  
  # Print summary results
  if (dist_type == "gamma") {
    cat("Gamma", shape_gamma, scale_gamma, "\n")
  } else if (dist_type == "lognormal") {
    cat("LN", mu_lognormal, sigma_lognormal, "\n")
  } else if (dist_type == "t") {
    cat("t", df_t, loc_t, scale_t, "\n")
  } else if (dist_type == "gpd") {
    cat("GPD", sigma_gpd, xi_gpd, "\n")
  }
  
  cat("n=",n,"\n")
  cat("Total attempts:", total_attempts, "\n")
  cat("Skewness Results:\n")
  print(results$Skewness)
  cat("Kurtosis Results:\n")
  kurtosis_results <- results$Kurtosis[!results$Kurtosis$Method %in% c("mc", "modal"), ]
  print(kurtosis_results)
  
  skewness_results <- results$Skewness
  kurtosis_results <- results$Kurtosis[!results$Kurtosis$Method %in% c("mc", "modal"), ]
  
}


### Section 4. mixture simulation ###
simulation_mixture_results <- function(n, dist1, dist2, shape1=1, scale1=1, shape2=1, scale2=1, 
                                       mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1, 
                                       df1 = NULL, df2 = NULL, xi1 = NULL, xi2 = NULL, weight1=0.5) {
  set.seed(1)
  nsim <- 1000 # Number of simulations
  
  # Compute true skewness, kurtosis
  true_moments <- calculate_mixture_moments(weight1, dist1, dist2, shape1, scale1, shape2, scale2, mu1, sigma1, mu2, sigma2, df1, df2, xi1, xi2)
  true_skewness <- true_moments$Skew_mix
  true_kurtosis <- true_moments$Kurt_mix
  true_values <- c(true_skewness, true_kurtosis)  # Store true skewness and kurtosis
  
  # Define estimation methods
  methods <- c("Empirical", "New1", "New5_smallresid", "lmom", "TL", "Quantile1",
               "Quantile2", "Quantile3", "mc", "modal", "KDE_Silverman", "KDE_Scott")
  stat_names <- c("Skewness", "Kurtosis")
  
  results <- list()
  estimates_storage <- list()
  
  # Create data frames to store bias, variance, and sqrt(MSE)
  for (stat in stat_names) {
    results[[stat]] <- data.frame(
      Method = methods,
      True_Value = rep(0, length(methods)),
      Bias = rep(0, length(methods)),
      Variance = rep(0, length(methods)),
      root_MSE = rep(0, length(methods)),
      Overestimate = rep(0, length(methods)),
      Underestimate = rep(0, length(methods)),
      number_of_failure = rep(0, length(methods))
    )
    results[[stat]]$True_Value <- ifelse(stat == "Skewness", true_skewness, true_kurtosis)
    estimates_storage[[stat]] <- matrix(NA, nrow = nsim, ncol = length(methods))
  }
  
  isim <- 0
  total_attempts <- 0
  failures <- rep(0, length(methods))
  
  while (isim < nsim) {
    total_attempts <- total_attempts + 1
    
    # Generate data according to weight1
    data <- generate_mixture_data(n, weight1, dist1, dist2, shape1, scale1, shape2, scale2, mu1, sigma1, mu2, sigma2, df1, df2, xi1, xi2)
    X <- data$X_mixture
    
    # Apply different estimation methods
    methods_res <- list(
      Empirical = list(alphas = list(alpha_3 = skewness(X), alpha_4 = kurtosis(X)), failure = 0),
      New1 = run_model(X, "New1"),
      New5_smallresid = run_new5_model(X),
      lmom = list(alphas = calculate_l_moments(X), failure = 0),
      TL = list(alphas = calculate_tl_moments(X), failure = 0),
      Quantile1 = list(alphas = calculate_quantile1_skewness_kurtosis(X), failure = 0),
      Quantile2 = list(alphas = calculate_quantile2_skewness_kurtosis(X), failure = 0),
      Quantile3 = list(alphas = calculate_quantile3_skewness_kurtosis(X), failure = 0),
      mc = list(alphas = calculate_mc_skewness(X), failure = 0),
      modal = list(alphas = calculate_modal_skewness(X), failure = 0),
      KDE_Silverman = list(alphas = calculate_kde_moments(X, method = "silverman"), failure = 0),
      KDE_Scott = list(alphas = calculate_kde_moments(X, method = "scott"), failure = 0)
    )
    
    if (methods_res$New5_smallresid$failure == 1) {
      failures[2] <- failures[2] + methods_res$New1$failure
      next
    }
    isim <- isim + 1
    #print(isim)
    
    # Store results for each method
    for (method_index in seq_along(methods)) {
      method <- methods[method_index]
      estimates <- methods_res[[method]]$alphas
      for (i in 1:2) {
        true_value <- true_values[i]
        estimate <- ifelse(i == 1, estimates$alpha_3, estimates$alpha_4)
        if (!is.na(estimate)) {
          bias_value <- estimate - true_value
          results[[stat_names[i]]]$Bias[method_index] <- 
            results[[stat_names[i]]]$Bias[method_index] + (estimate - true_value)
          results[[stat_names[i]]]$root_MSE[method_index] <- 
            results[[stat_names[i]]]$root_MSE[method_index] + (estimate - true_value)^2
          estimates_storage[[stat_names[i]]][isim, method_index] <- estimate
          if (bias_value > 0) {
            results[[stat_names[i]]]$Overestimate[method_index] <- 
              results[[stat_names[i]]]$Overestimate[method_index] + 1
          } else if (bias_value < 0) { # Clearly distinguish Overestimate and Underestimate
            results[[stat_names[i]]]$Underestimate[method_index] <- 
              results[[stat_names[i]]]$Underestimate[method_index] + 1
          }
          
        }
      }
      failures[method_index] <- failures[method_index] + methods_res[[method]]$failure
    }
  }
  
  # Finalize bias, variance, and sqrt(MSE) calculations
  for (stat in stat_names) {
    results[[stat]]$Bias <- results[[stat]]$Bias / nsim
    results[[stat]]$Variance <- apply(estimates_storage[[stat]], 2, var, na.rm = TRUE)
    results[[stat]]$root_MSE <- sqrt(results[[stat]]$Bias^2 + results[[stat]]$Variance)
    results[[stat]]$number_of_failure <- unname(failures)
  }
  
  
  cat("n=",n,"\n")
  cat("Total attempts:", total_attempts, "\n")
  cat("Skewness Results:\n")
  print(results$Skewness)
  cat("Kurtosis Results:\n")
  kurtosis_results <- results$Kurtosis[!results$Kurtosis$Method %in% c("mc", "modal"), ]
  print(kurtosis_results)

}


#============================================================

# Set up file for saving results
output_file <- "simulation_kde_results_250804.txt"
sink(output_file, split = TRUE)  # Output to both file and console simultaneously


sample_sizes <- c(50, 100, 200, 500, 1000, 2000)


###1 Component distribution###
# t-distribution (df = 5, mu = 0.0127*12, sigma = 0.0348*12)
loc_t <- 0.0127*12
scale_t <- 0.0348*12 / sqrt(5 / (5 - 2))  # Scaling for t-distribution with df=5
for (n in sample_sizes) {
  cat("\n# t-distribution (df = 5), n =", n, "\n")
  simulation_results(n = n, dist_type = "t", df_t = 5, loc_t = loc_t, scale_t = scale_t)
}

# Gamma(1, 1)
for (n in sample_sizes) {
  cat("\n# Gamma(1, 1), n =", n, "\n")
  simulation_results(n = n, dist_type = "gamma", shape_gamma = 1, scale_gamma = 1)
}

#Lognormal(1, sigma=0.5)
for (n in sample_sizes) {
  cat("\n# Lognormal (mu = 1, sigma = 0.5), n =", n, "\n")
  simulation_results(n = n, dist_type = "lognormal", mu_lognormal = 1, sigma_lognormal = 0.5)
}

# GPD with xi = 0.1, sigma = 2
for (n in sample_sizes) {
  cat("\n# GPD (xi = 0.1, sigma = 2), n =", n, "\n")
  simulation_results(n = n, dist_type = "gpd", sigma_gpd = 2, xi_gpd = 0.1)
}


######################
###mixture###

# t-distribution (df = 5) mixture with t-distribution (df = 5)
df_t <- 5

# Dist1 (t-dist): mu = 0.0127*12, sigma = 0.0348*12
mu_t1 <- 0.0127*12
sigma_t1 <- 0.0348*12
scale_t1 <- sigma_t1 / sqrt(df_t / (df_t - 2))

# Dist2 (t-dist): mu = -0.0161*12, sigma = 0.0748*12
mu_t2 <- -0.0161*12
sigma_t2 <- 0.0748*12
scale_t2 <- sigma_t2 / sqrt(df_t / (df_t - 2))

for (n in sample_sizes) {
  cat("\n# Mixture: t-distribution (df = 5, 0.0127*12, 0.0348*12) & t-distribution (df = 5, -0.0161*12, 0.0748*12), n =", n, "weight = 0.95","\n")
  simulation_mixture_results(n = n, 
                             dist1 = "t", df1 = df_t, mu1 = mu_t1, sigma1 = sigma_t1,  
                             dist2 = "t", df2 = df_t, mu2 = mu_t2, sigma2 = sigma_t2, 
                             weight1=0.95)
}


# Gamma(1, 1) mixture with Gamma(2,1)
for (n in sample_sizes) {
  cat("\n# Mixture: Gamma(1, 1) & Gamma(2,1), n =", n, "weight = 0.95","\n")
  simulation_mixture_results(n = n, dist1 = "gamma", shape1 = 1, scale1 = 1, 
                             dist2 = "gamma", shape2 = 2, scale2 = 1, weight1=0.95)
}

# GPD mixture with Gamma(2,1)
for (n in sample_sizes) {
  cat("\n# Mixture: GPD (xi = 0.1, sigma = 2) & Gamma(2,1), n =", n, "weight = 0.95","\n")
  simulation_mixture_results(n = n, dist1 = "gpd", sigma1 = 2, xi1 = 0.1, 
                             dist2 = "gamma", shape2 = 2, scale2 = 1, weight1=0.95)
}

# Lognormal distribution mixture with Gamma(2,1)
for (n in sample_sizes) {
  cat("\n# Mixture: Lognormal (mu = 1, sd = 0.5) & Gamma(2,1), n =", n, "weight = 0.95", "\n")
  simulation_mixture_results(n = n, dist1 = "lognormal", mu1 = 1, sigma1 = 0.5, 
                             dist2 = "gamma", shape2 = 2, scale2 = 1, weight1=0.95)
}


# Saving completed
sink()

library(dplyr)

# Skewness Method name change
all_results <- all_results %>%
  mutate(Method = case_when(
    moments == "Skewness" & Method == "New5_smallresid" ~ "New",
    moments == "Skewness" & Method == "Quantile1" ~ "Bowley",
    moments == "Skewness" & Method == "Quantile2" ~ "Groeneveld",
    moments == "Skewness" & Method == "Quantile3" ~ "Pearson",
    moments == "Skewness" & Method == "lmom" ~ "Lmom",
    moments == "Skewness" & Method == "modal" ~ "Modal",
    TRUE ~ Method
  ))

# Kurtosis Method  name change
all_results <- all_results %>%
  mutate(Method = case_when(
    moments == "Kurtosis" & Method == "New5_smallresid" ~ "New",
    moments == "Kurtosis" & Method == "Quantile1" ~ "Moors",
    moments == "Kurtosis" & Method == "Quantile2" ~ "Hogg",
    moments == "Kurtosis" & Method == "Quantile3" ~ "Crow",
    moments == "Kurtosis" & Method == "lmom" ~ "Lmom",
    TRUE ~ Method
  ))

#  name change for two-point mixture distns
all_mix_results <- all_mix_results %>%
  mutate(Method = case_when(
    moments == "Skewness" & Method == "New5_smallresid" ~ "New",
    moments == "Skewness" & Method == "Quantile1" ~ "Bowley",
    moments == "Skewness" & Method == "Quantile2" ~ "Groeneveld",
    moments == "Skewness" & Method == "Quantile3" ~ "Pearson",
    moments == "Skewness" & Method == "lmom" ~ "Lmom",
    moments == "Skewness" & Method == "modal" ~ "Modal",
    moments == "Kurtosis" & Method == "New5_smallresid" ~ "New",
    moments == "Kurtosis" & Method == "Quantile1" ~ "Moors",
    moments == "Kurtosis" & Method == "Quantile2" ~ "Hogg",
    moments == "Kurtosis" & Method == "Quantile3" ~ "Crow",
    moments == "Kurtosis" & Method == "lmom" ~ "Lmom",
    TRUE ~ Method
  ))

write_csv(all_results, "all_results.csv")
write_csv(all_mix_results, "all_mix_results.csv")

