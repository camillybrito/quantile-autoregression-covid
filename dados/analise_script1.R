# PACOTES UTILIZADOS ---------------------------------------------------------------
# install.packages("COVID19")
# install.packages("tidyverse")
# #instalando o pacote covid19 que contém dados sobre a pandemia
# install.packages("RSQLite")
# install.packages("stats")
# install.packages("dplyr")
# install.packages("quantreg")
# install.packages("ggplot2")
# install.packages("caret")
# install.packages("sf")

#https://covid19datahub.io/articles/doc/data.html
#Variables explanation
#Cran: https://cran.r-project.org/web/packages/COVID19/COVID19.pdf

library(COVID19) # carregando pacote
library(tidyverse)
library(qpcR) # for the Hannan Quinn Criterion
require(qpcR)
library(stats)
library(quantreg) #para regressao quantilica
library(caret)
library(tidyquant)
library(TSA)
library(ICglm) # para o HQIC
library(scales)
library(zoo)
library(geobr)
library(sf)

# ACESSANDO DADOS---------------------------------------------------------------
brcovid <- covid19(
  country = "BRA", 
  level = 2, 
  start = '2020-03-11', 
  end = '2020-12-31'
) %>% 
  arrange(date)  # ordena por data


# DADOS NACIONAIS (Nível 1)
covid_nacional <- covid19(
  country = "BRA",
  level = 1,
  start = "2020-02-01",
  end = "2022-12-31"
) %>% 
  arrange(date) %>%
  mutate(casos_diarios = confirmed - lag(confirmed))%>%
  mutate(
    mortes_diarias = deaths - lag(deaths)  # cálculo das mortes por dia
  )


# 
# brcovidnovo <- data.frame(
#   Data = brcovid$date,
#   Mortes = deaths_diariobr,
#   Confirmados = conf_diariobr,
#   Vacinados = vaccin_diariobr, 
#   # Estado = brcovid$administrative_area_level_2,
#   Lat = brcovid$latitude,
#   Long = brcovid$longitude,
#   País = brcovid$administrative_area_level_1
# )
# 
# brcovidnovo <- brcovidnovo %>% mutate_all(funs(replace(., is.na(.), 0)))
# # ou utilizar (criando outro banco):
# # bancobrasil <- replace(x = brcovidnovo, list = is.na(brcovidnovo), values = 0) 
# ## bancobrasil tem zero no lugar de NA

# 
# # Dados por estado
# covid_uf <- covid19("BRA", level = 2, start = "2020-02-01", end = "2022-12-31") %>%
#   group_by(administrative_area_level_2) %>%
#   summarise(
#     casos_confirmados = max(confirmed, na.rm = TRUE),
#     mortes_confirmadas = max(deaths, na.rm = TRUE)
#   ) %>%
#   rename(uf = administrative_area_level_2)

# GRÁFICO MAPA DE FREQUENCIA----------------------------------------------------
# Mapa dos estados brasileiros
mapa_uf <- read_state(year = 2020)

# Verificando nomes de estado
unique(mapa_uf$name_state)
unique(covid_uf$uf)

# Novo join por nome do estado
mapa_casos <- left_join(mapa_uf, covid_uf, by = c("name_state" = "uf"))

ggplot(mapa_casos) +
  geom_sf(aes(fill = casos_confirmados), color = "white") +
  scale_fill_gradient(
    low = "#cce5ff",  # azul claro
    high = "#08306b", # azul escuro
    labels = addUnits
  ) +
  labs(
    #title = "Casos Confirmados de COVID-19 por Estado (2020–2022)",
    fill = "Casos"
  ) +
  theme_minimal()

ggplot(mapa_casos) +
  geom_sf(aes(fill = mortes_confirmadas), color = "white") +
  scale_fill_gradient(low = "#ffd6d6", high = "#990000", labels = addUnits) +
  labs(
    #title = "Mortes Confirmadas por COVID-19 por Estado (2020–2022)", 
    fill = "Mortes") +
  theme_minimal()


# GRÁFICOS SERIES CASOS E MORTES----------------------------------------------

# Gráfico com destaque para picos (média móvel de 7 dias)
ggplot(covid_nac, aes(x = date)) +
  geom_col(aes(y = mortes_diarias), fill = "gray40", alpha = 0.6) +
  geom_line(aes(y = media7), color = "red", size = 1) +
  scale_y_continuous(labels = addUnits)+
  labs(
    #title = "Mortes Diárias por COVID-19 no Brasil",
    subtitle = "Com média móvel de 7 dias",
    x = "Data",
    y = "Número de Mortes"
  ) +
  theme_minimal()


ggplot(covid_nacional, aes(x = date, y = casos_diarios)) +
  geom_col(fill = "tomato") +
  scale_y_continuous(labels = addUnits) +
  labs(
    #title = "Casos Diários de COVID-19 no Brasil (Fev 2020 - Dez 2022)",
    x = "Data",
    y = "Número de Casos Diários"
  ) +
  theme_minimal()


ggplot(covid_nacional, aes(x = date, y = mortes_diarias)) +
  geom_col(fill = "gray20") +
  scale_y_continuous(labels = addUnits) +
  labs(
    #title = "Mortes Diárias por COVID-19 no Brasil (Fev 2020 - Dez 2022)",
    x = "Data",
    y = "Número de Mortes"
  ) +
  theme_minimal()



# ANÁLISE EXPLORATÓRIA DE SERIES TEMPORAIS----------------------------------------  

#modificando temporariamente os dados
brcovidnovo <- brcovidnovo %>% select(Data, Mortes, Confirmados)
str(brcovidnovo)

# transformar em um objeto de série temporal:
brcovidnovo_ts <- ts(brcovidnovo$Mortes, start = c(2020, 3, 11), 
                     frequency = 365.25)
is.ts(brcovidnovo_ts)
# decomposição aditiva
decomposicao <- decompose(brcovidnovo_ts)

# Visualizar os componentes
plot(decomposicao)

# componentes separadamente
# tendencia <- decomposicao$trend
# sazonalidade <- decomposicao$seasonal
# residuos <- decomposicao$residuals

# # plotando os componentes
# par(mfrow=c(3,1))
# plot(tendencia, main="Tendência", col="blue")
# plot(sazonalidade, main="Sazonalidade", col="red")
# plot(residuos, main="Resíduos", col="green")


