fn_cria_serie_acoplada <- function(n, phi, tau, mu) {
  # Verificar se n > tau
  if(n <= tau){
    return(NA)  # Retorna NA caso n <= tau
  }
  
  # Inicializar o vetor X_t
  X_t <- numeric(n)
  
  # Gerar o processo AR(1) para X_t
  X_t[1] <- rnorm(1)  # Valor inicial (aleatório)
  
  for (t in 2:n) {
    X_t[t] <- mu + phi * X_t[t - 1] + rnorm(1)  # Processo AR(1)
  }
  
  # Inicializar o vetor Y_t
  Y_t <- numeric(n)
  
  # Gerar o processo AR(1) para Y_t
  Y_t[1] <- rnorm(1)  # Valor inicial (aleatório)
  
  for (t in 2:n) {
    Y_t[t] <- mu + phi * Y_t[t - 1] + rnorm(1)  # Processo AR(1)
  }
  
  # Inicializar o vetor Z_t
  Z_t <- numeric(n)
  
  # Definir Z_t conforme a condição
  for (t in 1:n) {
    if (t <= tau) {
      Z_t[t] <- X_t[t]  # Para t <= τ, Z_t = X_t
    } else {
      Z_t[t] <- Y_t[t - tau]  # Para t > τ, Z_t = Y_{t - τ}
    }
  }
  
  # Média da série temporal
  media_zt = mean(Z_t)
  
  return(media_zt)
}

# Função de simulação de Monte Carlo para a média e variância amostral
simulacao_montecarlo <- function(n, phi, tau, mu, n_simulacoes) {
  # Vetor para armazenar as médias de Z_t em cada simulação
  medias_z <- numeric(n_simulacoes)
  
  # Realizar n_simulacoes e calcular a média de Z_t para cada uma
  for (i in 1:n_simulacoes) {
    # Chamar a função para gerar a série acoplada e calcular a média de Z_t
    medias_z[i] <- fn_cria_serie_acoplada(n, phi, tau, mu)
  }
  
  # Calcular a média das médias amostrais
  media_amostral <- mean(medias_z, na.rm = TRUE)  # Remove NAs, caso n <= tau em algum caso
  
  # Calcular a variância das médias amostrais
  variancia_amostral <- var(medias_z, na.rm = TRUE)  # Remove NAs ao calcular a variância
  
  # Retornar a média das médias, a variância das médias e as médias simuladas
  return(list(media = media_amostral, variancia = variancia_amostral, medias = medias_z))
}

# Definir os parâmetros
n_vals <- c(50, 100, 500, 1000, 5000)
phi <- 0.5
mu <- 0
n_simulacoes <- 100  # Número de simulações


combinacoes <- expand.grid(n = n_vals, tau = c(0.25, 0.5, 0.75))

# Calcular os valores de tau para cada n
combinacoes <- combinacoes %>%
  mutate(tau = round(tau * n))  # Multiplicar n por 0.25, 0.5, e 0.75

resultados <- combinacoes %>%
  rowwise() %>%
  mutate(
    resultados_simulacao = list(simulacao_montecarlo(n, phi, tau, mu, n_simulacoes))
  ) %>%
  # Usando unnest() para desmembrar a lista em colunas separadas
  unnest(cols = c(resultados_simulacao)) %>% View()
  rename(
    media_simulada = media,
    variancia_simulada = variancia,
    medias_simuladas = medias
  )
