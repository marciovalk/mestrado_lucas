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

impute_with_kalman <- function(ts) {
  # Argumentos:
  # ts: Série temporal com valores nulos
  
  # Imputação de valores com filtro de Kalman
  imputed_ts <- na_kalman(ts)
  
  return(imputed_ts)
}


extract_non_missing_segments <- function(ts) {
  # Argumentos:
  # ts: Série temporal com valores nulos
  
  # Identificando os termos nulos
  non_na_indices <- which(!is.na(ts))
  
  # Extrai nulos da série
  new_ts <- ts[non_na_indices]
  
  return(new_ts)
}


# Função para calcular o estimador de Mínimos Quadrados (OLS) para AR(1)
ols_estimate_ar1 <- function(ts) {
  n <- length(ts)
  y <- ts[2:n]
  x <- ts[1:(n-1)]
  
  # Estimador OLS para AR(1) é dada por phi_hat = (sum(x * y)) / sum(x^2)
  phi_hat <- sum(x * y) / sum(x^2)
  
  return(phi_hat)
}

# Função para calcular o Erro Quadrático Médio (EQM)
calculate_mse <- function(estimates, true_phi) {
  # Calcula o EQM
  mse <- mean((estimates - true_phi)^2)
  return(mse)
}

# Função para realizar a simulação de Monte Carlo
monte_carlo_simulation <- function(n_simulations, n, phi, sigma, k, j) {
  phi_estimates_imputed <- numeric(n_simulations)
  phi_estimates_non_missing <- numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    # Gerar série temporal AR(1) com NAs
    ts_with_na <- generate_ar1_with_na(n, phi, sigma, k, j)
    
    # Imputação com Kalman
    ts_imputed <- impute_with_kalman(ts_with_na)
    
    # Estimador OLS para a série imputada
    phi_estimates_imputed[i] <- ols_estimate_ar1(ts_imputed)
    
    # Extrair os segmentos não faltantes
    ts_non_missing <- extract_non_missing_segments(ts_with_na)
    
    # Estimador OLS para a série com segmentos não faltantes
    phi_estimates_non_missing[i] <- ols_estimate_ar1(ts_non_missing)
  }
  
  # Calcular o EQM para cada abordagem
  mse_imputed <- calculate_mse(phi_estimates_imputed, phi)
  mse_non_missing <- calculate_mse(phi_estimates_non_missing, phi)
  
  # Retorna as estimativas de phi e os EQMs
  return(list(
    phi_estimates_imputed = phi_estimates_imputed,
    phi_estimates_non_missing = phi_estimates_non_missing,
    mse_imputed = mse_imputed,
    mse_non_missing = mse_non_missing
  ))
}

# Configurações da simulação
n_simulations <- 1000  # Número de simulações de Monte Carlo
n <- 350  # Tamanho das séries temporais
phi <- 0.5  # Coeficiente AR(1)
sigma <- 1  # Desvio padrão do erro branco
k <- 30  # Índice de início dos NAs
j <- 90  # Quantidade consecutiva de NAs

# Executar a simulação de Monte Carlo
results <- monte_carlo_simulation(n_simulations, n, phi, sigma, k, j)

# Exibir as médias das estimativas de phi e os EQMs
cat("Média da estimativa de phi (com imputação Kalman): ", mean(results$phi_estimates_imputed), "\n")
cat("Média da estimativa de phi (sem NAs - segmentos não faltantes): ", mean(results$phi_estimates_non_missing), "\n")
cat("EQM para a série imputada com Kalman: ", results$mse_imputed, "\n")
cat("EQM para a série sem NAs (segmentos não faltantes): ", results$mse_non_missing, "\n")







