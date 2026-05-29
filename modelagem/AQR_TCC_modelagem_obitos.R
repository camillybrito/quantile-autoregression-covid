# Modelo autorregressivo quantílico --------------------------------------------
library(tidyverse)
library(qpcR) # Hannan Quinn Criterion
require(qpcR)
library(stats)
library(quantreg) # regressao quantilica
library(caret)
library(tidyquant)
library(TSA)
library(ICglm) # Hannan Quinn Criterion
library(scales)
library(zoo)
library(geobr)

# para número de mortes

# https://acervolima.com/regressao-de-quantis-na-programacao-r/
# https://eranraviv.com/quantile-autoregression-in-r/
# https://www.youtube.com/watch?v=cm29i-n53SM
# https://pt.stackoverflow.com/questions/71993/decomposi%C3%A7%C3%A3o-de-s%C3%A9ries-temporais-di%C3%A1rias

dados_mensais

dados_mensais_ts <- ts(dados_mensais$obitosMensais, start = c(2020, 02, 01), 
                     frequency = 12)
is.ts(dados_mensais_ts)

options(digits = 4)
library(qpcR) # for the Hannan Quinn Criterion
library(ICglm) # para o HQIC
require(qpcR)
y = dados_mensais_ts
plot(y, main = "Evolução das mortes por covid ao longo do tempo", ylab = "Quantidade", 
     lwd = 2,  type = "b", cex.main = 1.5)

lagorder = function(x,maxdef){ 
  # maxdef é o número max de defasagens utilizadas
  ordmax = maxdef ; HQ = NULL ; AK = NULL ; SC = NULL ; lagmat = NULL
  l1 = length(x)
  for (i in 1:ordmax){
    lagmat = cbind(lagmat[-i,],x[(1):(l1-i)]) # lagged matrix
    armod <- lm(x[(i+1):l1]~lagmat)
    HQ[i] = HQIC(armod)
    AK[i] = AIC(armod)
    SC[i] = BIC(armod)
  }
  return(c(which.min(HQ), which.min(AK),which.min(SC) ))
}

#  função chamada lagorder que é usada para determinar a ordem adequada de um modelo 
#  autoregressivo (AR) com base em diferentes critérios de informação, como o Critério de 
#  Informação de Hannan-Quinn (HQIC), o Critério de Informação de Akaike (AIC) e o Critério de 
#  Informação Bayesiano (BIC).


lagorder(y, maxdef = 12)

TT = length(y)
lm0 = lm(y[-1]~y[-TT]) ; summary(lm0)

tauseq = seq(.05,.95,.05)
qslope = NULL ; qr0 = list()
pvalor.0 = NULL; pvalor.1 = NULL
for (i in 1:length(tauseq)) {
  mod <- rq(y[-1] ~ y[-TT], tau = tauseq[i])
  qr0[[i]] <- mod
  resumo <- summary(mod, se = "nid")  # <-- ESSENCIAL
  pvalor.0[i] <- resumo$coefficients[1, 4]  # alpha0
  pvalor.1[i] <- resumo$coefficients[2, 4]  # alpha1
  qslope[i] <- mod$coef[2]
}

qr0 # para ver os coeficientes do modelo
pvalor.0  # p-valores de alpha0
pvalor.1  # p-valores de alpha1




# teste -  para cada quantil

mod = rq(y[-1]~y[-TT], tau = 0.95) # ajusta um modelo QAR(1) no quantil "tau = ***"   
res = mod$residuals  

Box.test(res, lag = 1, type = "Ljung-Box", fitdf = 0)   # testa se ainda existe autocorrelação
# Interpretação do teste Ljung-Box
# Hipóteses:
#   H0: resíduos não autocorrelacionados (modelo bem especificado)
#   H1: resíduos autocorrelacionados
# Resultado:
#   p-valor pequeno (< 0.05) modelo inadequado (falta estrutura)
#   p-valor grande modelo ok


# para p-valor muito pequeno use <0.0001




layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
layout.show(3)
plot(y, ylab = "Número de mortes", main  = "Evolução das mortes por covid ao longo do tempo no Brasil")
plot(qslope~tauseq, ty = "b", xlab = "Quantil", ylab = "Coeficiente AR")
plot(y[-1]~y[-TT], col = 2 , ylab = "Número de mortes", xlab = "Número de mortes (defasagens)")
for (i in 1:length(tauseq)) {
  lines(y[-TT], qr0[[i]]$fitted.values, 
        lwd = 1, col = rgb(0, 0, 0, 0.2))
}


# Se as linhas:
#   são paralelas → efeito constante
#   se abrem / cruzam → heterocedasticidade ou não linearidade



head(y[-1]); tail(y[-1]);
head(y[-TT]); tail(y[-TT]);


# tentativa de gráfico (p-valor INTERCEPTO x quantil)

gr <- data.frame(
  P_valor = pvalor.0,
  Quantil = tauseq
)

plot(gr$Quantil, gr$P_valor, main="P-valor do Intercepto para cada quantil analisado",
     xlab="Quantil", 
     ylab="P-valor")
abline(h = 0.05, col = "red", lty = 1)

text(x = min(gr$Quantil), y = 0.05, labels = "0.05", col = "red", pos = 3, offset = 0.5)



# tentativa de gráfico (p-valor y[-TT] x quantil)

gr <- data.frame(
  P_valor = pvalor.1,
  Quantil = tauseq
)

plot(gr$Quantil, gr$P_valor, main="P-valor para o parâmetro Autorregressivo de ordem 1
     para cada quantil analisado",
     xlab="Quantil", 
     ylab="P-valor")
abline(h = 0.05, col = "red", lty = 1)

text(x = min(gr$Quantil), y = 0.05, labels = "0.05", col = "red", pos = 3, offset = 0.5)



#---------------------------------------------
# Medidas de erro
#---------------------------------------------

# Vetores
mse  <- numeric(length(tauseq))
rmse <- numeric(length(tauseq))
mae  <- numeric(length(tauseq))
mape <- numeric(length(tauseq))

for (i in 1:length(tauseq)) {
  
  # Valores ajustados
  yhat <- qr0[[i]]$fitted.values
  
  # Valores observados
  y_real <- y[-1]
  
  # Erros
  erro <- y_real - yhat
  
  # MSE
  mse[i] <- mean(erro^2, na.rm = TRUE)
  
  # RMSE
  rmse[i] <- sqrt(mean(erro^2, na.rm = TRUE))
  
  # MAE
  mae[i] <- mean(abs(erro), na.rm = TRUE)
  
  # Índices válidos (remove zeros)
  idx_validos <- y_real != 0
  
  # MAPE corrigido
  mape[i] <- mean(
    abs(
      erro[idx_validos] /
        y_real[idx_validos]
    ),
    na.rm = TRUE
  ) * 100
}

# Tabela final
medidas_erro <- data.frame(
  Tau  = tauseq,
  MSE  = round(mse, 2),
  RMSE = round(rmse, 2),
  MAE  = round(mae, 2),
  MAPE = round(mape, 2)
)

print(medidas_erro)


# Quantis com melhor desempenho preditivo
quantis_centrais <- c(0.25, 0.30, 0.35, 0.40, 0.45, 0.50)

# Vetores para armazenar resultados
mape_centrais <- numeric(length(quantis_centrais))
mse_centrais  <- numeric(length(quantis_centrais))
rmse_centrais <- numeric(length(quantis_centrais))
mae_centrais  <- numeric(length(quantis_centrais))

for(i in seq_along(quantis_centrais)) {
  
  # Índice do quantil correspondente
  idx <- which(abs(tauseq - quantis_centrais[i]) < 1e-8)
  
  # Valores ajustados
  yhat <- qr0[[idx]]$fitted.values
  
  # Valores observados
  y_real <- y[-1]
  
  # Erro
  erro <- y_real - yhat
  
  # Índices válidos (remove zeros)
  idx_validos <- y_real != 0
  
  # MAPE ignorando zeros
  mape_centrais[i] <- mean(
    abs(
      erro[idx_validos] /
        y_real[idx_validos]
    ),
    na.rm = TRUE
  ) * 100
  
  # MSE
  mse_centrais[i] <- mean(erro^2, na.rm = TRUE)
  
  # RMSE
  rmse_centrais[i] <- sqrt(mean(erro^2, na.rm = TRUE))
  
  # MAE
  mae_centrais[i] <- mean(abs(erro), na.rm = TRUE)
}

# Resultado final
resultado_mape <- data.frame(
  Tau  = quantis_centrais,
  MSE  = round(mse_centrais, 2),
  RMSE = round(rmse_centrais, 2),
  MAE  = round(mae_centrais, 2),
  MAPE = round(mape_centrais, 2)
)

print(resultado_mape)
