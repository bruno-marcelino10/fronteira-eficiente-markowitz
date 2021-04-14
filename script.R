##### Análise de Carteiras - Fronteira Eficiente de Markowitz ##############
##### 08/02/2021 ###########################################################
##### Bruno Marcelino ######################################################
##### Análise aplicada a diversas empresas entre 2010 e 2019 ###############

# Importação de bibliotecas... 

library("tidyverse")
library("tidyquant")
library("Quandl")
library("timetk")

#### ---- Importação e Tratamento dos Dados ----

## ESCOLHA OS ATIVOS E O PERÍODO

periodo <- as.Date("2009-12-31") 

ativo_1 <- "ABEV3.SA"
ativo_2 <- "ALPA4.SA"
ativo_3 <- "ITUB4.SA"
ativo_4 <- "CYRE3.SA"
ativo_5 <- "VALE3.SA"
ativo_6 <- "BTOW3.SA"
ativo_7 <- "WEGE3.SA"
ativo_8 <- "ITSA4.SA"
ativo_9 <- "BBAS3.SA"
ativo_10 <- "BBDC3.SA"
ativo_11 <- "GGBR4.SA"
ativo_12 <- "CSNA3.SA"
ativo_13 <- "SAPR4.SA"

# Ativos que serão importados (acrescente quantos desejar)

ativos <- c(ativo_1,
            ativo_2,
            ativo_3,
            ativo_4,
            ativo_5,
            ativo_6,
            ativo_7,
            ativo_8,
            ativo_9,
            ativo_10,
            ativo_11,
            ativo_12,
            ativo_13
)

# ---- Importando cotações ----

cotacoes <- ativos %>% 
    tq_get(get = "stock.prices", from = periodo, periodicity = "monthly") %>% 
    mutate(adjusted = na.locf(adjusted)) %>% 
    group_by(symbol)

# ---- Calculando os retornos para cada ativo ----

df_retornos <- cotacoes %>% 
    tq_transmute(select = adjusted,
                 mutate_fun = periodReturn,
                 period = "monthly",
                 type = "arithmetic",
                 col_rename = "retornos")

# ---- Importando a taxa SELIC da API do BACEN (de 2010 até hoje) ----

Quandl.api_key('TC1ow5j6G7s4SFHTzgDz') # Acrescentando a chave da API - Comando necessário pra acessar o Quandl
selic_mensal <-  as.data.frame(Quandl("BCB/4390", type = "xts")/100) # Importando a serie do selic do Bacen
selic_nova <- selic_mensal %>% 
    rownames_to_column() %>% 
    tibble() %>% 
    rename("date" = rowname, "retorno" = V1)

# Manipulando dados SELIC
meses_do_periodo <- (12*11)+1
df_selic <- slice(selic_nova, length(selic_nova$retorno)-meses_do_periodo+1:length(selic_nova$retorno))

media_selic <- df_selic %>% summarise(mean(retorno))
media_selic <- as.numeric(media_selic)


# ---- Separando df amostral (80% dos dados) e de avaliação (20% dos dados)

periodo_amostral <- round(0.8*length(df_retornos$date), 0)

data_amostra <- periodo %m+% months(periodo_amostral)

df_amostral <- df_retornos %>%  # 106 meses
    filter(date <= data_amostra)  # 28 meses

df_avaliacao <- df_retornos %>%  
    filter(date >= data_amostra)

# ----  Estatísticas Básicas do Portfólio ----

media_retornos <- df_amostral %>% summarise(mean(retornos)) # retorno médio no período de amostra
media_retornos <- as.numeric(media_retornos$`mean(retornos)`)

cov_retornos <- df_amostral %>% # covariância dos retornos do período de amostra
    spread(key = symbol, value = retornos) %>% 
    select(-date) %>% 
    cov()

dp_retornos <- cov_retornos %>% # desvio-padrão de cada ativo no período de amostra
    diag() %>% 
    sqrt()

# ----- Estimando o Portfólio Tangente -----

cov_inversa_retornos <- solve(cov_retornos) # calcula a inversa da matriz de covariâncias
excesso_de_retorno <- media_retornos - media_selic # calcula o vetor de prêmios de risco
portfólio_tangente <- cov_inversa_retornos %*% excesso_de_retorno # multiplica a matriz da inversa da covariância pela matriz dos prêmios de risco
pesos_portfólio_tangente <- portfólio_tangente/sum(portfólio_tangente)

# ----- Estimando o Portfólio de Variância Mínima -----

vetor_de_uns <- rep(1,length(ativos)) # cria um vetor de 1's
portfólio_variancia_minima <- cov_inversa_retornos %*% vetor_de_uns # multiplica a matriz inversa da covariância pelo vetor de uns
pesos_portfólio_variancia_minima <- portfólio_variancia_minima/sum(portfólio_variancia_minima) 


# ---- Fazendo uma tabela comparativa entre PVM e PT ---- 
estatisticas_descritivas_do_portfolio <- function(pesos){ # função que calcula as estatísticas descritivas, dado um portfólio já calculado
    retorno <- sum(pesos*media_retornos)
    dp <- sqrt((t(pesos) %*% (cov_retornos %*% pesos)))
    sharpe <- (retorno - media_selic)/dp
    estatisticas <- cbind(retorno,sharpe,dp)
    return(estatisticas)
}

# Manipulando os dados dos pesos do PVM e PT
pesos_portfólio_variancia_minima <- t(pesos_portfólio_variancia_minima)
pesos_portfólio_variancia_minima <- c(as.numeric(pesos_portfólio_variancia_minima))

pesos_portfólio_tangente <- t(pesos_portfólio_tangente)
pesos_portfólio_tangente <- c(as.numeric(pesos_portfólio_tangente))

# Criando vetores com o retorno e o risco do PVM e PT
PVM <- as.numeric(estatisticas_descritivas_do_portfolio(pesos_portfólio_variancia_minima))
PT <- as.numeric(estatisticas_descritivas_do_portfolio(pesos_portfólio_tangente))

# ----- Criando a Fronteira Eficiente -----

numero_portfolios_eficientes_testados <- 10000
portfolios_eficientes <- cbind("PT" = pesos_portfólio_tangente, "PVM" = pesos_portfólio_variancia_minima)
combinacoes <- seq(from = 0, by = 1/numero_portfolios_eficientes_testados)[-1]

# Matrizes que irão receber as estatísticas dos i portfólios eficientes que serão criados
matriz_pesos_FE <- matrix(nrow = numero_portfolios_eficientes_testados, ncol = length(ativos))
matriz_retornos_FE <- vector("numeric", length = numero_portfolios_eficientes_testados)
matriz_riscos_FE <- vector('numeric', length = numero_portfolios_eficientes_testados)
matriz_sharpe_FE <- vector('numeric', length = numero_portfolios_eficientes_testados)

for (i in seq_along(matriz_retornos_FE)){
    
    peso_PVM <- runif(n = 1, min = -0.5, max = 1) # proporção a ser aplicada no PVM (número aleatório entre 0 e 1)
    peso_PT <- 1-peso_PVM   
    pesos_FE <- (pesos_portfólio_variancia_minima*peso_PVM) + (pesos_portfólio_tangente*peso_PT) # portfólios eficientes são uma combinação linear do PVM e do PT
    pesos_FE <- t(pesos_FE)
    pesos_FE <- c(as.numeric(pesos_FE))
    matriz_pesos_FE[i,] <- pesos_FE # gera uma matriz com i portfólios eficientes possíveis
    
    retorno_FE <- sum(pesos_FE*media_retornos) 
    matriz_retornos_FE[i] <- retorno_FE # calcula o retorno de cada portfólio eficiente i
    
    var_portfolio_FE <- t(pesos_FE) %*% (cov_retornos %*% pesos_FE) 
    dp_portfolio_FE <- sqrt(var_portfolio_FE)
    matriz_riscos_FE[i] <- dp_portfolio_FE # calcula o risco de cada portfólio eficiente i
    
    sharpe_portfolio_FE <- (retorno_FE - media_selic)/dp_portfolio_FE
    matriz_sharpe_FE[i] <- sharpe_portfolio_FE #calcula o sharpe de cada portfólio eficiente i
}

# Unindo os resultados obtidos para cada portfólio eficiente i em uma tabela
fronteira_eficiente <- tibble("Retornos" = matriz_retornos_FE,
                              "Sharpes" = matriz_sharpe_FE,
                              "Desvio-padrão" = matriz_riscos_FE
)

matriz_pesos_FE <- matriz_pesos_FE %>% tk_tbl()

names(matriz_pesos_FE) <- c(ativos)

fronteira_eficiente <- fronteira_eficiente %>% bind_cols(matriz_pesos_FE)

# ---- Gerando Gráfico da Fronteira Eficiente, destacando o desempenho das carteiras não-diversificadas ----

plot(fronteira_eficiente$`Desvio-padrão`,
     fronteira_eficiente$Retornos,
     xlab = "Risco",
     ylab = "Retorno",
     main = "Fronteira Eficiente de Markowitz",
     xlim = c(0, 0.15),
     ylim = c(0, 0.031)
)

vet_risco <- c(PVM[3], PT[3])
vet_retorno <- c(PVM[1], PT[1])

risco_retorno_portfolios <- data.frame(Risco = c(vet_risco, dp_retornos), Retorno = c(vet_retorno, media_retornos))
risco_retorno_portfolios <- t(risco_retorno_portfolios)
colnames(risco_retorno_portfolios) <- c("Variância Mínima", "Portfólio Tangente", ativos)

points(t(risco_retorno_portfolios), pch = 19, col = 1:length(ativos)+2)
text(t(risco_retorno_portfolios), colnames(risco_retorno_portfolios), col= 1:length(ativos)+2, cex = 1.3, pos = 4) 
grid()