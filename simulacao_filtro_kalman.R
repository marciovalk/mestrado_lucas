generate_ar1_with_na <- function(n, phi, sigma, k, j) {
  # Argumentos:
  # n: Tamanho do processo
  # phi: AR(1) Coeficiente (entre -1 e 1)
  # sigma: SD do erro branco
  # k: Index de início dos NAs
  # j: Qtd consecutiva de NAs
  
  # Corrigindo inputs inválidos
  if (k + j - 1 > n) {
    stop("The range k to k+j-1 exceeds the length of the time series n.")
  }
  
  if (phi <= -1 || phi >= 1) {
    stop("phi must be between -1 and 1 for stationarity.")
  }
  
  # Gera o processo
  ar1 <- numeric(n)
  noise <- rnorm(n, mean = 0, sd = sigma)
  
  ar1[1] <- noise[1]  
  for (t in 2:n) {
    ar1[t] <- phi * ar1[t - 1] + noise[t]
  }
  
  # Introduzindo os NAs
  ar1[k:(k + j - 1)] <- NA
  
  return(ar1)
}

# Imputação com filtro de Kalman
impute_with_kalman <- function(ts) {
  imputed_ts <- na_kalman(ts)
  return(imputed_ts)
}

# Imputação por Interpolação Linear
impute_with_linear_interpolation <- function(ts) {
  imputed_ts <- na_interpolation(ts, option = "linear")
  return(imputed_ts)
}

# Imputação por Média Móvel
impute_with_moving_average <- function(ts, window_size = 5) {
  imputed_ts <- na_ma(ts, k = window_size, weighting = "simple")
  return(imputed_ts)
}

# Imputação com a Média Global
impute_with_global_mean <- function(ts) {
  ts[is.na(ts)] <- mean(ts, na.rm = TRUE)
  return(ts)
}

# Imputação por Hot Deck
impute_with_hot_deck <- function(ts) {
  library(VIM)
  # Converte a série temporal para data frame
  ts_df <- data.frame(value = ts)
  imputed_df <- hotdeck(ts_df)
  imputed_ts <- imputed_df$value
  return(imputed_ts)
}


extract_non_missing_segments <- function(ts) {
  non_na_indices <- which(!is.na(ts))
  new_ts <- ts[non_na_indices]
  return(new_ts)
}

ols_estimate_ar1 <- function(ts) {
  n <- length(ts)
  y <- ts[2:n]
  x <- ts[1:(n-1)]
  
  phi_hat <- sum(x * y) / sum(x^2)
  return(phi_hat)
}

calculate_mse <- function(estimates, true_phi) {
  mse <- mean((estimates - true_phi)^2)
  return(mse)
}

monte_carlo_simulation <- function(n_simulations, n, phi, sigma, k, j) {
  phi_estimates_kalman <- numeric(n_simulations)
  phi_estimates_linear <- numeric(n_simulations)
  phi_estimates_ma <- numeric(n_simulations)
  phi_estimates_mean <- numeric(n_simulations)
  phi_estimates_hot_deck <- numeric(n_simulations)
  phi_estimates_non_missing <- numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    ts_with_na <- generate_ar1_with_na(n, phi, sigma, k, j)
    
    # Imputação com Kalman
    ts_kalman <- impute_with_kalman(ts_with_na)
    phi_estimates_kalman[i] <- ols_estimate_ar1(ts_kalman)
    
    # Imputação por interpolação linear
    ts_linear <- impute_with_linear_interpolation(ts_with_na)
    phi_estimates_linear[i] <- ols_estimate_ar1(ts_linear)
    
    # Imputação por média móvel
    ts_ma <- impute_with_moving_average(ts_with_na)
    phi_estimates_ma[i] <- ols_estimate_ar1(ts_ma)
    
    # Imputação por média global
    ts_mean <- impute_with_global_mean(ts_with_na)
    phi_estimates_mean[i] <- ols_estimate_ar1(ts_mean)
    
    # Imputação por Hot Deck
    ts_hot_deck <- impute_with_hot_deck(ts_with_na)
    phi_estimates_hot_deck[i] <- ols_estimate_ar1(ts_hot_deck)
    
    # Estimativa sem dados faltantes (segmentos não faltantes)
    ts_non_missing <- extract_non_missing_segments(ts_with_na)
    phi_estimates_non_missing[i] <- ols_estimate_ar1(ts_non_missing)
  }
  
  # Calcula o EQM para cada abordagem
  mse_kalman <- calculate_mse(phi_estimates_kalman, phi)
  mse_linear <- calculate_mse(phi_estimates_linear, phi)
  mse_ma <- calculate_mse(phi_estimates_ma, phi)
  mse_mean <- calculate_mse(phi_estimates_mean, phi)
  mse_hot_deck <- calculate_mse(phi_estimates_hot_deck, phi)
  mse_non_missing <- calculate_mse(phi_estimates_non_missing, phi)
  
  return(list(
    phi_estimates_kalman = phi_estimates_kalman,
    phi_estimates_linear = phi_estimates_linear,
    phi_estimates_ma = phi_estimates_ma,
    phi_estimates_mean = phi_estimates_mean,
    phi_estimates_hot_deck = phi_estimates_hot_deck,
    phi_estimates_non_missing = phi_estimates_non_missing,
    mse_kalman = mse_kalman,
    mse_linear = mse_linear,
    mse_ma = mse_ma,
    mse_mean = mse_mean,
    mse_hot_deck = mse_hot_deck,
    mse_non_missing = mse_non_missing
  ))
}

n_simulations <- 1000
n <- 300
phi <- 0.5
sigma <- 1
k <- 20
j <- 120

results <- monte_carlo_simulation(n_simulations, n, phi, sigma, k, j)

cat("Média da estimativa de phi (Kalman):", mean(results$phi_estimates_kalman), "\n")
cat("Média da estimativa de phi (Interpolação Linear):", mean(results$phi_estimates_linear), "\n")
cat("Média da estimativa de phi (Média Móvel):", mean(results$phi_estimates_ma), "\n")
cat("Média da estimativa de phi (Média Global):", mean(results$phi_estimates_mean), "\n")
cat("Média da estimativa de phi (Hot Deck):", mean(results$phi_estimates_hot_deck), "\n")
cat("Média da estimativa de phi (Sem NAs):", mean(results$phi_estimates_non_missing), "\n")

cat("EQM (Kalman):", results$mse_kalman, "\n")
cat("EQM (Interpolação Linear):", results$mse_linear, "\n")
cat("EQM (Média Móvel):", results$mse_ma, "\n")
cat("EQM (Média Global):", results$mse_mean, "\n")
cat("EQM (Hot Deck):", results$mse_hot_deck, "\n")
cat("EQM (Sem NAs):", results$mse_non_missing, "\n")
