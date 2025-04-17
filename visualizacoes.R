gerar_arma_modificada <- function(n, p, q, x_1, tau) {
  # Carregar os pacotes necessários
  if (!requireNamespace("forecast", quietly = TRUE)) {
    install.packages("forecast")
  }
  library(forecast)
  
  if (!requireNamespace("imputeTS", quietly = TRUE)) {
    install.packages("imputeTS")
  }
  library(imputeTS)
  
  # Gerar coeficientes aleatórios para AR e MA
  ar_coeffs <- runif(p, -0.9, 0.9)  
  ma_coeffs <- runif(q, -0.9, 0.9)  
  
  # Criar o processo ARMA
  ts_arma <- arima.sim(n = n, model = list(ar = ar_coeffs, ma = ma_coeffs))
  
  # Criar uma cópia para modificação
  ts_arma_mod <- ts_arma
  
  # Criar a lacuna de valores NA entre x_1 e x_1 + tau
  if (x_1 + tau <= n) {  
    ts_arma_mod[(x_1 + 1):(x_1 + tau)] <- NA
  } else {
    warning("x_1 + tau excede o tamanho da série. Ajustando para o máximo possível.")
    ts_arma_mod[(x_1 + 1):n] <- NA
  }
  
  # Imputação dos dados ausentes por outros métodos
  ts_media <- na_replace(ts_arma_mod, mean(ts_arma_mod, na.rm = TRUE))  
  ts_mediana <- na_replace(ts_arma_mod, median(ts_arma_mod, na.rm = TRUE))  
  ts_kalman <- na_kalman(ts_arma_mod)  
  
  # Hot Deck Aleatório: substituir NA por valores aleatórios da série original
  valores_disponiveis <- ts_arma[!is.na(ts_arma_mod)]  
  ts_hotdeck <- ts_arma_mod
  ts_hotdeck[is.na(ts_hotdeck)] <- sample(valores_disponiveis, sum(is.na(ts_hotdeck)), replace = TRUE)
  
  # Forward Fill: preenche com o último valor observado
  ts_forward <- na_locf(ts_arma_mod)
  
  # Backward Fill: preenche com o próximo valor observado
  ts_backward <- na_locf(ts_arma_mod, option = "nocb")  
  
  # Imputação por Filtro Gaussiano
  imputar_filtro_gaussiano <- function(ts, kernel_size = 3, sigma = 1) {
    kernel <- dnorm(seq(-kernel_size, kernel_size), mean = 0, sd = sigma)  
    kernel <- kernel / sum(kernel)  
    
    ts_gaussiano <- ts  
    nas <- which(is.na(ts))  
    
    for (i in nas) {
      vizinhos <- max(1, i - kernel_size):min(length(ts), i + kernel_size)
      valores_vizinhos <- ts[vizinhos]
      pesos <- kernel[(kernel_size + 1) - (i - vizinhos)]
      
      valores_validos <- !is.na(valores_vizinhos)
      
      if (sum(valores_validos) > 0) {
        ts_gaussiano[i] <- sum(valores_vizinhos[valores_validos] * pesos[valores_validos]) / sum(pesos[valores_validos])
      } else {
        ts_gaussiano[i] <- mean(ts, na.rm = TRUE)  
      }
    }
    
    return(ts_gaussiano)
  }
  
  ts_gaussiano <- imputar_filtro_gaussiano(ts_arma_mod, kernel_size = 3, sigma = 1)
  
  # Cole a série modificada, juntando as partes e excluindo os NAs
  ts_mod_colada <- ts_arma_mod[!is.na(ts_arma_mod)]
  
  # Criar um dataframe com as séries imputadas (com colunas de tamanho n)
  df_series <- data.frame(
    Tempo = 1:n, 
    Original = ts_arma, 
    Modificada = ts_arma_mod,
    Media = ts_media, 
    Mediana = ts_mediana,
    HotDeck = ts_hotdeck, 
    Kalman = ts_kalman,
    ForwardFill = ts_forward,
    BackwardFill = ts_backward,
    Gaussiano = ts_gaussiano
  )
  
  # Retorna uma lista contendo o dataframe completo e a série modificada "colada"
  return(list(df_series = df_series, Modificada_Colada = ts_mod_colada))
}

# Exemplo de uso
set.seed(123)  
resultado <- gerar_arma_modificada(n = 500, p = 1, q = 1, x_1 = 150, tau = 60)
series_df <- resultado$df_series
ts_mod_colada <- resultado$Modificada_Colada

# Plotar as séries
par(mfrow = c(3, 3)) 

plot(series_df$Tempo, series_df$Original, type = "l", col = "blue", lwd = 2,
     main = "Série Original", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$Modificada, type = "l", col = "red", lwd = 2,
     main = "Série Modificada (com NAs)", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$Media, type = "l", col = "green", lwd = 2,
     main = "Imputação pela Média", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$Mediana, type = "l", col = "purple", lwd = 2,
     main = "Imputação pela Mediana", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$HotDeck, type = "l", col = "orange", lwd = 2,
     main = "Imputação Hot Deck", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$Kalman, type = "l", col = "cyan", lwd = 2,
     main = "Imputação Kalman", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$ForwardFill, type = "l", col = "brown", lwd = 2,
     main = "Imputação Forward Fill", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$BackwardFill, type = "l", col = "pink", lwd = 2,
     main = "Imputação Backward Fill", ylab = "Valor", xlab = "Tempo")

plot(series_df$Tempo, series_df$Gaussiano, type = "l", col = "darkgreen", lwd = 2,
     main = "Imputação Filtro Gaussiano", ylab = "Valor", xlab = "Tempo")

par(mfrow = c(1, 1))

