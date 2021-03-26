#Homework by Jonibek Isomiddinov
install.packages("rio")
library(tidyverse)
library(purrr)
library(urca)
library(rugarch)
library(forecast)
library(rio)
library(vars)
library(tseries)
library(tsibble)
new_data <- import("PS.xls")
glimpse(new_data)
df <- ts(new_data$Euro)
st <- log(df) 
r <- diff(st)
acf(st, lag.max = 50)#not decreasing - clearly non-stationary
acf(r, lag.max = 50) #looks more like stationary
pacf(st, lag.max = 50)
pacf(r, lag.max = 50) 

# adf test: nodiff
adf0.1 <- ur.df(r, type = "none", selectlags = c("BIC")) # urca package
summary(adf0.1)
adf0.2 <- ur.df(r, type = "drift", selectlags = c("BIC"))
summary(adf0.2)
adf0.3 <- ur.df(r, type = "trend", selectlags = c("BIC"))
summary(adf0.3)

#Null hypothesis in ur.df is that there is unit root, however for all specifications
#we reject it at any reasonable significance level -> returns are stationary

PP.test(r, lshort = TRUE)
PP.test(r, lshort = FALSE)

#Null of Phillips-Perron Test that r has a unit root is rejected at 5% significance
#level

#returns are stationary, now let's look at log exchange rates
adf.test(st) #null of unit root is not rejected

PP.test(st, lshort = TRUE)# non-stationary
PP.test(st, lshort = FALSE)# non-stationary

# adf test: nodiff
adf0.1 <- ur.df(st, type = "none", selectlags = c("BIC")) # urca package
summary(adf0.1)
adf0.2 <- ur.df(st, type = "drift", selectlags = c("BIC"))
summary(adf0.2)
adf0.3 <- ur.df(st, type = "trend", selectlags = c("BIC"))
summary(adf0.3)
#for each of the specifications we do not reject the null of unit root, our 
#test statistic is not lower that critical values, i would argue for different specifications
#if there were contradictions among them but they are all clearly stating same result

# adf test: 1 diff
adf1.1 <- ur.df(diff(st, 1), type = "none", selectlags = c("BIC"))
summary(adf1.1)
adf1.2 <- ur.df(diff(st, 1), type = "drift", selectlags = c("BIC"))
summary(adf1.2)
adf1.3 <- ur.df(diff(st, 1), type = "trend", selectlags = c("BIC"))
summary(adf1.3)
#very good results for first differences, now the null of unit root is rejected for sure
#our process is I(1)
adf.test(diff(st)) #just to check all is ok
PP.test(diff(st)) #just to check all is ok 

#our log exchange rates are I(1)

# II Ljung Box Test
Box.test(r, lag = 20, fitdf = 0, type = "Lj")
Box.test(r, lag = 40, fitdf = 0, type = "Lj")
Box.test(r, lag = 60, fitdf = 0, type = "Lj")
#seems that null of independece is rejected for all specifications at 5% sign level
#not independent series

#III Box-Jenkins
mod_1 <- auto.arima(r)
summary(mod_1) #best model chosen by Box-Jenkins methodology is ARIMA (0,0,1)
AIC(mod_1) #AIC of chosen model
BIC(mod_1) #BIC of chosen model

# IV any heteroskedasticity?
acf(mod_1$residuals^2)
pacf(mod_1$residuals^2) 
#we have some arch effect
Box.test(mod_1$residuals^2, type = "Ljung-Box") #null of independence is rejected


# This script will estimate several garch models and find the best using the BIC
# criteria. 

#First I initialize parameters for model, then run through a loop to build 
#models with different specifications, finally from whole bunch of models
#the ones are chosen with least bic and aic

## MAIN OPTIONS

max_lag_AR <- 1 
max_lag_MA <- 1 
max_lag_ARCH <- 2 
max_lag_GARCH <- 1 
dist_to_use <- c('norm', 'std') # see rugarch::ugarchspec help for more
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch::rugarchspec help for more

## END OPTIONS

do_arch_test <- function(x, max_lag = 5) {
  require(FinTS)
  require(tidyverse)
  
do_single_arch <- function(x, used_lag)  {
    test_out <- FinTS::ArchTest(x, lags = used_lag)
    
    res_out <- tibble(Lag = used_lag,
                      `LMStatistic` = test_out$statistic, 
                      `pvalue` = test_out$p.value)}

out <- find_best_arch_model(x = r, 
                            type_models = models_to_estimate,
                            dist_to_use = dist_to_use,
                            max_lag_AR = max_lag_AR,
                            max_lag_MA = max_lag_MA,
                            max_lag_ARCH = max_lag_ARCH,
                            max_lag_GARCH = max_lag_GARCH)


find_best_arch_model <- function(x, 
                                 type_models, 
                                 dist_to_use,
                                 max_lag_AR,
                                 max_lag_MA,
                                 max_lag_ARCH,
                                 max_lag_GARCH) {
  
  require(tidyr)
  
  df_grid <- expand_grid(type_models = type_models,
                         dist_to_use = dist_to_use,
                         arma_lag = 0:max_lag_AR,
                         ma_lag = 0:max_lag_MA,
                         arch_lag = 1:max_lag_ARCH,
                         garch_lag = 1:max_lag_GARCH)
  
  
  l_out <- pmap(.l = list(x = rep(list(x), nrow(df_grid)), 
                          type_model = df_grid$type_models,
                          type_dist = df_grid$dist_to_use,
                          lag_ar = df_grid$arma_lag,
                          lag_ma = df_grid$ma_lag,
                          lag_arch = df_grid$arch_lag,
                          lag_garch  = df_grid$garch_lag),
                do_single_garch)
  
  tab_out <- bind_rows(l_out)
  
  # find by AIC
  idx <- which.min(tab_out$AIC)
  best_aic <- tab_out[idx, ]
  
  # find by BIC
  idx <- which.min(tab_out$BIC)
  best_bic <- tab_out[idx, ]
  
  l_out <- list(best_aic = best_aic,
                best_bic = best_bic,
                tab_out = tab_out)
  
  return(l_out)
}
do_single_garch <- function(x, 
                            type_model, 
                            type_dist, 
                            lag_ar, 
                            lag_ma, 
                            lag_arch, 
                            lag_garch) {
  require(rugarch)
  
  
  spec = ugarchspec(variance.model = list(model =  type_model, 
                                          garchOrder = c(lag_arch, lag_garch)),
                    mean.model = list(armaOrder = c(lag_ar, lag_ma)),
                    distribution = type_dist)
  
  message('Estimating ARMA(',lag_ar, ',', lag_ma,')-',
          type_model, '(', lag_arch, ',', lag_garch, ')', 
          ' dist = ', type_dist,
          appendLF = FALSE)
  
  try({
    my_rugarch <- list()
    my_rugarch <- ugarchfit(spec = spec, data = x)
  })
  
  if (!is.null(coef(my_rugarch))) {
    message('\tDone')
    
    AIC <- rugarch::infocriteria(my_rugarch)[1]
    BIC <- rugarch::infocriteria(my_rugarch)[2]
  } else {
    message('\tEstimation failed..')
    
    AIC <- NA
    BIC <- NA
  }
  
  est_tab <- tibble(lag_ar, 
                    lag_ma,
                    lag_arch,
                    lag_garch,
                    AIC =  AIC,
                    BIC = BIC,
                    type_model = type_model,
                    type_dist,
                    model_name = paste0('ARMA(', lag_ar, ',', lag_ma, ')+',
                                        type_model, '(', lag_arch, ',', lag_garch, ') ',
                                        type_dist) ) 
  
  return(est_tab)
}


out$best_aic
out$best_bic

#we can see that the best model is ARIMA (0,0,1) + Garch(1,1)

# ARIMA + GARCH just fit
garch.spec <- ugarchspec(variance.model = list(model = "sGARCH", # can be different functional form
                                               garchOrder = c(1, 1)), # can be different GARCH order
                         mean.model = list(armaOrder = c(0, 1)), # can be different ARMA order
                         distribution.model = "norm") # can be different distribution 

garch.model <- ugarchfit(spec = garch.spec, data = r)
garch.model



#!!!!! for b repeating all this for the S&P series
dff <- ts(new_data$SP500)
SP <- log(dff) 
SPr <- diff(SP)
acf(SP, lag.max = 50)#not decreasing - clearly non-stationary
acf(SPr, lag.max = 50) #looks more like stationary
pacf(SP, lag.max = 50)
pacf(SPr, lag.max = 50) 

# adf test: nodiff
adf0.1 <- ur.df(SPr, type = "none", selectlags = c("BIC")) # urca package
summary(adf0.1)
adf0.2 <- ur.df(SPr, type = "drift", selectlags = c("BIC"))
summary(adf0.2)
adf0.3 <- ur.df(SPr, type = "trend", selectlags = c("BIC"))
summary(adf0.3)

#Null hypothesis in ur.df is that there is unit root, however for all specifications
#we reject it at any reasonable significance level -> returns are stationary

PP.test(SPr, lshort = TRUE)
PP.test(SPr, lshort = FALSE)

#Null of Phillips-Perron Test that SPr has a unit root is rejected at 5% significance
#level

#returns are stationary, now let's look at log exchange rates
adf.test(SP) #null of unit root is not rejected

PP.test(SP, lshort = TRUE)# non-stationary
PP.test(SP, lshort = FALSE)# non-stationary

# adf test: nodiff
adf0.1 <- ur.df(SP, type = "none", selectlags = c("BIC")) # urca package
summary(adf0.1)
adf0.2 <- ur.df(SP, type = "drift", selectlags = c("BIC"))
summary(adf0.2)
adf0.3 <- ur.df(SP, type = "trend", selectlags = c("BIC"))
summary(adf0.3)
#for each of the specifications we do not reject the null of unit root (5% sign), our 
#test statistic is not lower that critical values, i would argue for different specifications
#if there were contradictions among them but they are all clearly stating same result

# adf test: 1 diff
adf1.1 <- ur.df(diff(SP, 1), type = "none", selectlags = c("BIC"))
summary(adf1.1)
adf1.2 <- ur.df(diff(SP, 1), type = "drift", selectlags = c("BIC"))
summary(adf1.2)
adf1.3 <- ur.df(diff(SP, 1), type = "trend", selectlags = c("BIC"))
summary(adf1.3)
#very good results for first differences, now the null of unit root is rejected for sure
#our process is I(1)
adf.test(diff(SP)) #just to check all is ok
PP.test(diff(SP)) #just to check all is ok 

#our log exchange rates are I(1)

# II Ljung Box Test
Box.test(SPr, lag = 20, fitdf = 0, type = "Lj")
Box.test(SPr, lag = 40, fitdf = 0, type = "Lj")
Box.test(SPr, lag = 60, fitdf = 0, type = "Lj")
#seems that null of independece is rejected for all specifications at 5% sign level
#not independent series

#III Box-Jenkins
mod_2 <- auto.arima(SPr)
summary(mod_2) #best model chosen by Box-Jenkins methodology is ARIMA (0,0,2)
AIC(mod_2) #AIC of chosen model
BIC(mod_2) #BIC of chosen model

# IV any heteroskedasticity?
acf(mod_2$residuals^2)
pacf(mod_2$residuals^2) 
#we have some arch effect
Box.test(mod_2$residuals^2, type = "Ljung-Box") #null of independence is rejected


# This script will estimate several garch models and find the best using the BIC
# criteria. 

#First I initialize parameters for model, then run through a loop to build 
#models with different specifications, finally from whole bunch of models
#the ones are chosen with least bic and aic

## MAIN OPTIONS

max_lag_AR <- 2
max_lag_MA <- 2 
max_lag_ARCH <- 2 
max_lag_GARCH <- 2 
dist_to_use <- c('norm', 'std') # see rugarch::ugarchspec help for more
models_to_estimate <- c('sGARCH', 'eGARCH', 'gjrGARCH') # see rugarch::rugarchspec help for more

## END OPTIONS
out <- find_best_arch_model(x = SPr, 
                              type_models = models_to_estimate,
                              dist_to_use = dist_to_use,
                              max_lag_AR = max_lag_AR,
                              max_lag_MA = max_lag_MA,
                              max_lag_ARCH = max_lag_ARCH,
                              max_lag_GARCH = max_lag_GARCH)


out$best_aic
out$best_bic
  
#we can see that the best model is ARIMA (1,0,1) + Garch(2,1) by aic
#and ARIMA(0,0,1) + Garch(2,1) by bic
  
# ARIMA + GARCH just fit
garch.spec <- ugarchspec(variance.model = list(model = "sGARCH", # can be different functional form
                                                 garchOrder = c(2, 1)), # can be different GARCH order
                           mean.model = list(armaOrder = c(1, 1)), # can be different ARMA order
                           distribution.model = "norm") # can be different distribution 
  
garch.model <- ugarchfit(spec = garch.spec, data = SPr)
garch.model


## C testing for breaks!
# LOG exchange rates --> st and log returns on SP500 ---> SPr
library(strucchange)
d_levels <- ts.intersect(y = st, y1 = stats::lag(st))

# BREAKS IN LEVELS
# consider AR(1) model
# compute a series of F statistics
qlr_levels <- Fstats(y ~ y1, from = 0.15, to = 0.85, data = d_levels)

# plot F-stats
plot(qlr_levels)
lines(breakpoints(qlr_levels))

# widen the window
qlr_levels1 <- Fstats(y ~ y1, from = 0.1, to = 0.9, data = d_levels)
plot(qlr_levels1)
breakpoints(qlr_levels1)
lines(breakpoints(qlr_levels1))

#break at 473
# for arima models - compute F-sats by hand
# consult https://stats.stackexchange.com/questions/104324/strucchange-package-on-arima-model for more info
arima_order <- c(0,0,1)
n <- length(st)

# fit unrestricted model
fit <- Arima(st, order = arima_order, include.mean = FALSE)
rss_all <- sum(residuals(fit)^2)
sigma2 <- fit$sigma2

# now fit separate models for each period
stats <- rep(NA, n)
for (i in seq.int(20, n-20)) {
  fit1 <- arima(st[seq(1,i)], order = arima_order, include.mean = FALSE)
  fit2 <- arima(st[seq(i+1,n)], order = arima_order, include.mean = FALSE)
  rss_sum <- sum(c(residuals(fit1), residuals(fit2))^2)
  stats[i] <- (rss_all - rss_sum)/sigma2
}

# plot F-sats and 95% quantile of chi-square distribution
plot(stats)
abline(h = qchisq(0.05, df = length(coef(fit)), lower.tail = FALSE), lty = 2, col = "red")

# find breakpoint as time that corresponds to the min p_value
which.min(1 - pchisq(stats, df = length(coef(fit))))

#here we see that after arima fitting the break is at 1442

# BREAKS IN VARIANCE
# back to AR(1) model
model <- lm(y ~ y1, data = d_levels)
summary(model)

unif_dum = rep(0, 3517)   #for this case 3517 for next with SPr 3516 i guess bcz of different orders of arima
#in case of error please put 3516 or vice versa in next step
unif_dum[473] = 1 # breakpoint number based on AR(1) model

# bind 2 time series: squared model residuals and dummy
d_var <- ts.intersect(y = model$residuals^2, y1 = ts(unif_dum))

# QLR test for new regression: model residuals ~ dummy
qlr_var <- Fstats(y ~ y1, from = 0.1, to = 0.9, data = d_var)
plot(qlr_var)
breakpoints(qlr_var)
lines(breakpoints(qlr_var))

#1423 break in variance very close to previous found

#now all the same but for SPr
d_levels <- ts.intersect(y = SPr, y1 = stats::lag(SPr))

# BREAKS IN LEVELS
# consider AR(1) model
# compute a series of F statistics
qlr_levels <- Fstats(y ~ y1, from = 0.15, to = 0.85, data = d_levels)

# plot F-stats
plot(qlr_levels)
lines(breakpoints(qlr_levels))

# widen the window
qlr_levels1 <- Fstats(y ~ y1, from = 0.1, to = 0.9, data = d_levels)
plot(qlr_levels1)
breakpoints(qlr_levels1)
lines(breakpoints(qlr_levels1))

#break at 1153
# for arima models - compute F-sats by hand
# consult https://stats.stackexchange.com/questions/104324/strucchange-package-on-arima-model for more info
arima_order <- c(0,0,2)
n <- length(SPr)

# fit unrestricted model
fit <- Arima(SPr, order = arima_order, include.mean = FALSE)
rss_all <- sum(residuals(fit)^2)
sigma2 <- fit$sigma2

# now fit separate models for each period
stats <- rep(NA, n)
for (i in seq.int(20, n-20)) {
  fit1 <- arima(SPr[seq(1,i)], order = arima_order, include.mean = FALSE)
  fit2 <- arima(SPr[seq(i+1,n)], order = arima_order, include.mean = FALSE)
  rss_sum <- sum(c(residuals(fit1), residuals(fit2))^2)
  stats[i] <- (rss_all - rss_sum)/sigma2
}

# plot F-sats and 95% quantile of chi-square distribution
plot(stats)
abline(h = qchisq(0.05, df = length(coef(fit)), lower.tail = FALSE), lty = 2, col = "red")

# find breakpoint as time that corresponds to the min p_value
which.min(1 - pchisq(stats, df = length(coef(fit))))

#here we see that after arima fitting the break is at 1149

# BREAKS IN VARIANCE
# back to AR(1) model
model <- lm(y ~ y1, data = d_levels)
summary(model)

unif_dum = rep(0, 3516)
unif_dum[473] = 1 # breakpoint number based on AR(1) model

# bind 2 time series: squared model residuals and dummy
d_var <- ts.intersect(y = model$residuals^2, y1 = ts(unif_dum))

# QLR test for new regression: model residuals ~ dummy
qlr_var <- Fstats(y ~ y1, from = 0.1, to = 0.9, data = d_var)
plot(qlr_var)
breakpoints(qlr_var)
lines(breakpoints(qlr_var))
#break at 2678

# NOW D!

d_levels <- ts.intersect(SP, st)
VARselect(d_levels, lag.max=7, type = "const")
VARselect(d_levels, lag.max = 7, type = "const")$selection
# eatimateed 3 lags by AIC

# estimate
model <- VAR(d_levels, p = 3, type = "const")
summary(model)

# granger causality (H_0 - doesn't granger cause)
causality(model, 'SP')
causality(model, 'st')

#we see that null is not rejected for both sides, neither of them Granger causes each other
#As both of the series are I(1) separately, there could be cointegration
#K = 3 as P = 3 by AIC

coint_h11 <- ca.jo(d_levels, type = "eigen", K = 3, spec = "transitory")
summary(coint_h11) #r=0 => no cointegration
# Determining the cointegration rank 
# r* is the number of cointegration vectors
# H_0: r = r* < k
# H_a: r = r* + 1
# testing proceeds sequentially for r* = 1,2,.. 
# with the first non-rejection used as an estimator for r
# since p=2, thus K=2
#we see that test statistic is lower than critical value, so no cointegration at all
#moreover both eigenvalues are ~~0
coint_h12 <- ca.jo(d_levels, type = "eigen", K = 3, spec = "longrun")
summary(coint_h12) #r=0 => no cointegration

#i received that there is not cointegration,
#but if i received that there is 1 cointegrating equation i would do as below code
#however in my case there would be enough just to compute var for both returns
#why returns as var requires both series stationary
d_levels <- ts.intersect(SPr, r)
VARselect(d_levels, lag.max=7, type = "const")
VAR(d_levels, lag = 2, type = "const")

#### THIS IS WHAT I WOULD DO IN CASE I HAD COINTEGRATION (lets say 1 cointegrating vector)
## eigenvectors and cointegration relations
coint_h12@V

################################################################################
################################## VECM ########################################
################################################################################

# error correction term
ec <- d_levels[, 'SP'] + coint_h12@V[2,1]*d_levels[, 'st']
autoplot(ec) + ggtitle("Error correction term")
tsdisplay(ec)

# VECM estimation
usconsumption_var2 <- VAR(d_levels, p = 3, type = "const")
usconsumption_vecm <- cajorls(coint_h12, r = 1) 
usconsumption_var <- vec2var(coint_h12, r = 1) 