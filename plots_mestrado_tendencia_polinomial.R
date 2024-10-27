# Script para gerar gráficos da seção 2.2

library(dplyr)
library(ggplot2)


# 1. Gera série temporal com tendência polinomial 

gerar_serie_temporal <- function(n, intercepto = 10, coeficiente = 1, var = 1) {
  # n: número de pontos na série temporal
  # intercepto: valor inicial da série
  # coeficiente: coeficiente da tendência (slope)
  # ruído: desvio padrão do ruído aleatório
  
  # Gerar o tempo
  tempo <- 1:n
  
  # Gerar a tendência linear
  tendencia <- intercepto + coeficiente * tempo
  
  # Adicionar ruído aleatório
  ruido <- rnorm(n, mean = 0, sd = var)
  
  # Criar a série temporal
  serie_temporal <- tendencia + ruido
  
  return(serie_temporal)
}


# Função para gerar série temporal com tendência polinomial e variáveis tau e k
gerar_serie_temporal_modificado <- function(n, tau, k, intercepto = 10, coeficiente = 1, var = 1) {
  # Gerar a série temporal original
  serie_temporal <- gerar_serie_temporal(n, intercepto, coeficiente, var)
  
  # Criar a nova série
  nova_serie <- c(serie_temporal[1:tau], rep(NA, k))
  
  # Adicionar os valores a partir de X_{tau + k + 1}
  if ((tau + k + 1) <= n) {
    nova_serie <- c(nova_serie, serie_temporal[(tau + k + 1):n])
  }
  
  return(nova_serie)
}

# Gerando a série 
# Para reprodutibilidade

set.seed(301263) 
n <- 100
tau <- 30
k <- 25

nova_serie <- gerar_serie_temporal_modificado(n, tau, k, intercepto = 10, coeficiente = 0.5, var = 1)

# Plotar a nova série usando ggplot2

dados <- data.frame(Tempo = 1:length(nova_serie), Valor = nova_serie)

plot_1 <- 
  ggplot(dados, aes(x = Tempo, y = Valor)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_point(data = dados %>% filter(is.na(Valor)), aes(y = Valor), color = "red", size = 3) +
  labs(title = "",
       x = "Tempo",
       y = "Valor") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(plot = plot_1, filename = 'tendencia_dados_faltantes.png')


# Série modificada para acoplamento 


gerar_serie_temporal_acoplamento <- function(n, tau, k, intercepto = 10, coeficiente = 1, var = 1) {
  # Gerar a série temporal original
  serie_temporal <- gerar_serie_temporal(n, intercepto, coeficiente, var)
  
  # Criar a nova série
  nova_serie <- c(serie_temporal[1:tau])
  
  # Adicionar os valores a partir de X_{tau + k + 1}
  if ((tau + k + 1) <= n) {
    kappa <- serie_temporal[tau + k + 1] - serie_temporal[tau]
    nova_serie <- c(nova_serie, serie_temporal[(tau + k + 1):n] - kappa)
  }
  
  return(nova_serie)
}


set.seed(301263) 
n <- 100
tau <- 30
k <- 25

nova_serie <- gerar_serie_temporal_acoplamento(n, tau, k, intercepto = 10, coeficiente = 0.5, var = 1)

# Preparar os dados para o ggplot2
dados <- data.frame(Tempo = 1:length(nova_serie), Valor = nova_serie)

# Plotar a nova série usando ggplot2
plot_2 <- 
  ggplot(dados, aes(x = Tempo, y = Valor)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  geom_point(data = dados %>% filter(is.na(Valor)), aes(y = Valor), color = "red", size = 3) +
  labs(title = "",
       x = "Tempo",
       y = "Valor") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))


ggsave(plot = plot_2, filename = 'tendencia_dados_faltantes_acoplada.png')
