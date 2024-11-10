# Carregar a biblioteca necessária
library(stats)
library(ggplot2)  # Para gráficos mais avançados

# Definir a função para criar a série temporal AR(1) combinada
fn_cria_ar_colado <- function(n, tau, phi) {
  # Gerar duas séries AR(1) com o mesmo phi
  X_t <- arima.sim(model = list(ar = phi), n = n)
  Y_t <- arima.sim(model = list(ar = phi), n = n)
  
  # Selecionar segmentos de cada série
  X_segmento <- X_t[1:tau]
  Y_segmento <- Y_t[1:(n - tau)]
  
  # Combinar os segmentos para criar Z_t
  Z_t <- c(X_segmento, Y_segmento)
  
  return(Z_t)
}

# Estimador para phi usando o método de MLE
estimate_phi_combined <- function(Z_t, tau) {
  # Dividir Z_t em segmentos X_t e Y_t
  X_segment <- Z_t[1:tau]                 # X_t
  Y_segment <- Z_t[(tau + 1):length(Z_t)] # Y_{t - \tau}
  
  # Calcular o numerador e o denominador do estimador
  numerator <- sum(head(X_segment, -1) * tail(X_segment, -1)) + 
    sum(head(Y_segment, -1) * tail(Y_segment, -1))
  denominator <- sum(head(X_segment, -1)^2) + sum(head(Y_segment, -1)^2)
  
  # Calcular phi_hat
  phi_hat <- numerator / denominator
  
  return(phi_hat)
}

# Função de simulação de Monte Carlo
monte_carlo_simulation <- function(n, tau, phi, n_simulations) {
  phi_hats <- numeric(n_simulations)
  
  for (i in 1:n_simulations) {
    Z_t <- fn_cria_ar_colado(n, tau, phi)
    phi_hats[i] <- estimate_phi_combined(Z_t, tau)
  }
  
  # Cálculo das estatísticas
  phi_mean <- mean(phi_hats)
  phi_bias <- phi_mean - phi        # Viés em relação ao verdadeiro phi
  phi_variance <- var(phi_hats)
  
  return(list(
    phi_estimates = phi_hats,
    phi_mean = phi_mean,
    phi_bias = phi_bias,
    phi_variance = phi_variance
  ))
}

# Função para gerar o histograma e o gráfico de densidade
generate_phi_plot <- function(phi_estimates, n, tau_factor) {
  # Definir o nome da combinação de n e tau
  plot_title <- paste("Histograma e Densidade para n =", n, "e tau =", tau_factor * 100, "%")
  
  # Criar o histograma e gráfico de densidade
  p <- ggplot(data.frame(phi_hat = phi_estimates), aes(x = phi_hat)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    geom_density(color = "red", size = 1) +
    labs(title = plot_title, x = expression(hat(phi)), y = "Densidade") +
    theme_minimal()
  
  # Exibir o gráfico
  print(p)
}

# Parâmetros para a simulação
phi <- 0.5           # Coeficiente AR(1) verdadeiro para ambas as séries
n_simulations <- 1000  # Número de simulações de Monte Carlo

# Valores de n e tau
n_values <- c(50, 100, 250, 500, 1000, 5000, 10000)
tau_factors <- c(0.25, 0.5, 0.75)

# Tabela para armazenar os resultados
results <- data.frame(
  n = integer(),
  tau_factor = numeric(),
  phi_mean = numeric(),
  phi_bias = numeric(),
  phi_variance = numeric()
)

# Executar a simulação para diferentes valores de n e tau
for (n in n_values) {
  for (tau_factor in tau_factors) {
    tau <- round(tau_factor * n)  # Calcular tau com base no fator
    
    # Executar a simulação de Monte Carlo
    simulation_results <- monte_carlo_simulation(n, tau, phi, n_simulations)
    
    # Armazenar os resultados
    results <- rbind(results, data.frame(
      n = n,
      tau_factor = tau_factor,
      phi_mean = simulation_results$phi_mean,
      phi_bias = simulation_results$phi_bias,
      phi_variance = simulation_results$phi_variance
    ))
    
    # Gerar e exibir o gráfico de histograma e densidade
    generate_phi_plot(simulation_results$phi_estimates, n, tau_factor)
  }
}

# Exibir a tabela de resultados
print(results)
