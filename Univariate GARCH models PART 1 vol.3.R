#============================================================================
#                            UNIVARIATE GARCH MODELS PART 1 OF 2
#============================================================================

# load libraries
# install.packages(c("rugarch","FinTS","KernSmooth"), repos="http://cran.us.r-project.org")
#
# Install FinTS from the downloaded zipped gz file:
# https://cran.r-project.org/src/contrib/Archive/FinTS/
# install.packages("C:/Users/Dean/Downloads/FinTS_0.4-5.tar.gz", repos = NULL, type = "source")

library(PerformanceAnalytics)
library(quantmod)
library(rugarch)
library(car)
library(FinTS)



options(digits=4)

# Download data
symbol.vec = c("MSFT", "^GSPC")
getSymbols(symbol.vec, from ="2000-01-03", to = "2012-04-03")
colnames(MSFT)
start(MSFT)
end(MSFT)

# Extract adjusted closing prices
MSFT = MSFT[, "MSFT.Adjusted", drop=F]
GSPC = GSPC[, "GSPC.Adjusted", drop=F]

# Plot prices
plot(MSFT)
plot(GSPC)

# Calculate log-returns for GARCH analysis
MSFT.ret = CalculateReturns(MSFT, method="log")
GSPC.ret = CalculateReturns(GSPC, method="log")


# Remove first NA observation
MSFT.ret = MSFT.ret[-1,]
GSPC.ret = GSPC.ret[-1,]
colnames(MSFT.ret) ="MSFT"
colnames(GSPC.ret) = "GSPC"

# Create combined data series
MSFT.GSPC.ret = cbind(MSFT.ret,GSPC.ret)

# Plot returns
plot(MSFT.ret)
plot(GSPC.ret)

# Plot returns with squared and absolute returns
dataToPlot = cbind(MSFT.ret, MSFT.ret^2, abs(MSFT.ret))
colnames(dataToPlot) = c("Returns", "Returns^2", "abs(Returns)")
plot.zoo(dataToPlot, main="MSFT Daily Returns", col="blue")

dataToPlot = cbind(GSPC.ret, GSPC.ret^2, abs(GSPC.ret))
colnames(dataToPlot) = c("Returns", "Returns^2", "abs(Returns)")
plot.zoo(dataToPlot, main="GSPC Daily Returns", col="blue")

# Plot autocorrelations of returns, returns^2 and abs(returns)
par(mfrow=c(3,2))
  acf(MSFT.ret, main="MSFT Returns")
  acf(GSPC.ret, main="GSPC Returns")
  acf(MSFT.ret^2, main="MSFT Returns^2")
  acf(GSPC.ret^2, main="GSPC Returns^2")
  acf(abs(MSFT.ret), main="MSFT abs(Returns)")
  acf(abs(GSPC.ret), main="GSPC abs(Returns)")
par(mfrow=c(1,1))

# Compute summary statistics
table.Stats(MSFT.GSPC.ret)


# SIMULATE ARCH(1) PROCESS
# Use rugarch function ugarchsim
?ugarchspec
?ugarchpath

# Specify arch(1) model
arch1.spec = ugarchspec(variance.model = list(garchOrder=c(1,0)),
                        mean.model = list(armaOrder=c(0,0)),
                        fixed.pars=list(mu = 0, omega=0.1, alpha1=0.8))
class(arch1.spec)
arch1.spec

set.seed(123)
arch1.sim = ugarchpath(arch1.spec, n.sim=1000)
# result is an S4 object
class(arch1.sim)
# [1] "uGARCHpath"
# attr(,"package")
# [1] "rugarch"
slotNames(arch1.sim)
# [1] "path"  "model" "seed"
names(arch1.sim@path)
# [1] "sigmaSim"  "seriesSim" "residSim"

# Use the plot method to plot simulated series and conditional volatilities
par(mfrow=c(2,1))
plot(arch1.sim, which=2)
plot(arch1.sim, which=1)
par(mfrow=c(1,1))

par(mfrow=c(3,1))
  acf(arch1.sim@path$seriesSim, main="Returns")
  acf(arch1.sim@path$seriesSim^2, main="Returns^2")
  acf(abs(arch1.sim@path$seriesSim), main="abs(Returns)")
par(mfrow=c(1,1))

# Use qqPlot() function from car package
qqPlot(arch1.sim@path$seriesSim, ylab="ARCH(1) Returns")


# SIMULATE GARCH(1,1) PROCESS
# specify GARCH(1,1) model
garch11.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                        mean.model = list(armaOrder=c(0,0)),
                        fixed.pars=list(mu = 0, omega=0.1, alpha1=0.1,
                                        beta1 = 0.7))
set.seed(123)
garch11.sim = ugarchpath(garch11.spec, n.sim=1000)

# use the plot method to plot simulated series and conditional volatilities
par(mfrow=c(2,1))
  plot(garch11.sim, which=2)
  plot(garch11.sim, which=1)
par(mfrow=c(1,1))

par(mfrow=c(3,1))
  acf(garch11.sim@path$seriesSim, main="Returns")
  acf(garch11.sim@path$seriesSim^2, main="Returns^2")
  acf(abs(garch11.sim@path$seriesSim), main="abs(Returns)")
par(mfrow=c(1,1))

# use qqPlot() function from car package
qqPlot(garch11.sim@path$seriesSim, ylab="GARCH(1,1) Returns")



# TESTING FOR ARCH/GARCH EFFECTS IN RETURNS

# Use Box.test from stats package
Box.test(coredata(MSFT.ret^2), type="Ljung-Box", lag = 12)
Box.test(coredata(GSPC.ret^2), type="Ljung-Box", lag = 12)

# use ArchTest() function from FinTS package for Engle's LM test
ArchTest(MSFT.ret)
ArchTest(GSPC.ret)


# ESTIMATE GARCH(1,1)

# Specify GARCH(1,1) model with only constant in mean equation
garch11.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                          mean.model = list(armaOrder=c(0,0)))
MSFT.garch11.fit = ugarchfit(spec=garch11.spec, data=MSFT.ret,
                             solver.control=list(trace = 1))

# uGARCHfit Objects
class(MSFT.garch11.fit)
slotNames(MSFT.garch11.fit)
names(MSFT.garch11.fit@fit)
names(MSFT.garch11.fit@model)

# Show garch fit
MSFT.garch11.fit

# Use extractor functions
#
# Function        Description
# coef()          Extract estimated coefficients
# infocriteria()  Calculate information criteria for fit
# likelihood()    Extract likelihood
# nyblom()        Calculate Hansen-Nyblom coefficient stability test
# signbias()      Calculate Engle-Ng sign bias test
# newsimpact()    Calculate news impact curve
# as.data.frame() Extract data, fitted data, residuals and conditional vol
# sigma()         Extract conditional volatility estimates
# residuals()     Extract residuals
# fitted()        Extract fitted values
# getspec()       Extract model specification
# gof()           Compute goodness-of-fit statistics
# uncmean()       Extract unconditional mean
# uncvariance()   Extract unconditional variance
# plot()          Produce various plots
# persistence()   Calculate persistence of fitted model
# halflife()      Calculate half-life of the fitted model

# Estimated coefficients
coef(MSFT.garch11.fit)
# Unconditional mean in mean equation
uncmean(MSFT.garch11.fit)
# Unconditional variance: omega/(alpha1 + beta1)
uncvariance(MSFT.garch11.fit)
# Persistence: alpha1 + beta1
persistence(MSFT.garch11.fit)
# Half-life: this is the the number of days it takes for
# half of the expected reversion back towards the unconditional variance
halflife(MSFT.garch11.fit)

# Residuals: e(t)
plot.ts(residuals(MSFT.garch11.fit), ylab="e(t)", col="blue")
abline(h=0)

# Sigma(t) = conditional volatility
plot.ts(sigma(MSFT.garch11.fit), ylab="sigma(t)", col="blue")

# uGARCHfit plot  methods:
# 1:  Series with 2 Conditional SD Superimposed
# 2:  Series with 2.5% VaR Limits (with unconditional mean)
# 3:  Conditional SD
# 4:  ACF of Observations
# 5:  ACF of Squared Observations
# 6:  ACF of Absolute Observations
# 7:  Cross Correlation
# 8:  Empirical Density of Standardized Residuals
# 9:  QQ-Plot of Standardized Residuals
# 10: ACF of Standardized Residuals
# 11: ACF of Squared Standardized Residuals
# 12: News-Impact Curve

#plot(MSFT.garch11.fit)
plot(MSFT.garch11.fit, which=1)
plot(MSFT.garch11.fit, which="all")
plot(MSFT.garch11.fit, which=9)


# SIMULATE FROM FITTED MODEL
MSFT.garch11.sim = ugarchsim(MSFT.garch11.fit,
                             n.sim=nrow(MSFT.ret),
                             rseed=123,
                             startMethod="unconditional")
class(MSFT.garch11.sim)
slotNames(MSFT.garch11.sim)

# Plot actual returns and simulated returns
par(mfrow=c(2,1))
  plot(MSFT.ret, main="Actual MSFT returns")
  plot(as.xts(MSFT.garch11.sim@simulation$seriesSim,
              order.by=index(MSFT.ret)),
       main="Simulated GARCH(1,1) Returns")
par(mfrow=c(1,1))




# MODEL SELECTION ON GARCH(P,Q) MODELS
arch.order = 1:5
arch.names = paste("arch", arch.order, sep="")

# Fit all arch models with p <= 5
infocriteria_all=NULL
arch.list = list()
for (p in arch.order) {
  arch.spec = ugarchspec(variance.model = list(garchOrder=c(p,0)),
                         mean.model = list(armaOrder=c(0,0)))
  arch.fit = ugarchfit(spec=arch.spec, data=MSFT.ret,
                       solver.control=list(trace = 0))
  arch.list[[p]] = arch.fit
  infocriteria_all = cbind(infocriteria_all, infocriteria(arch.fit))
}
names(arch.list) = arch.names

# Refit GARCH(1,1) to cleaned data
garch11.fit = ugarchfit(spec=garch11.spec, data=MSFT.ret,
                        solver.control=list(trace = 0))
infocriteria_all = cbind(infocriteria_all, infocriteria(garch11.fit))
colnames(infocriteria_all) = c(arch.names,"garch11")
infocriteria_all


# Compare arch(5) to garch(1,1)
par(mfrow=c(2,1))
plot.ts(sigma(garch11.fit), main="GARCH(1,1) conditional vol",
        ylab="vol", col="blue")
plot.ts(sigma(arch.list$arch5), main="ARCH(5) conditional vol",
        ylab="vol", col="blue")
par(mfrow=c(1,1))


# FORECASTING
# Forecast from GARCH(1,1)
MSFT.garch11.fcst = ugarchforecast(MSFT.garch11.fit, n.ahead=100)
class(MSFT.garch11.fcst)
slotNames(MSFT.garch11.fcst)
names(MSFT.garch11.fcst@forecast)

# Forecast object method functions
# Function      Description
# as.array      Extracts the forecast array
# as.data.frame Extracts the forecasts
# as.list       Extracts the forecast list will all rollframes
# plot          Forecasts plots
# fpm           Forecast performance measures
# show          Forecast summary

MSFT.garch11.fcst
#plot(MSFT.garch11.fcst)

# Plot selections:
# 0: Exit
# 1: Time Series Prediction (unconditional)
# 2: Time Series Prediction (rolling)
# 3: Conditional SD Prediction

par(mfrow=c(2,1))
  plot(MSFT.garch11.fcst, which=1)
  plot(MSFT.garch11.fcst, which=3)
par(mfrow=c(1,1))


# Forecast from ARCH(5)
MSFT.arch5.fcst = ugarchforecast(arch.list$arch5, n.ahead=100)
par(mfrow=c(2,1))
  plot(MSFT.arch5.fcst, which=1)
  plot(MSFT.arch5.fcst, which=3)
par(mfrow=c(1,1))





# Computing VaR
# look at example to see how to implement test
?VaRTest
# we use the original estimates without the cleaning
plot(MSFT.garch11.fit, which=2)

# BOOTSTRAP FORECAST DENSITIES
MSFT.garch11.boot = ugarchboot(MSFT.garch11.fit, method="Partial",
                               n.ahead=100, n.bootpred=2000)
class(MSFT.garch11.boot)
MSFT.garch11.boot
plot(MSFT.garch11.boot, which=3)


# ROLLING ESTIMATION OF GARCH(1,1)
MSFT.garch11.roll = ugarchroll(garch11.spec, MSFT.ret, n.ahead=1,
                               forecast.length = 1000,
                               refit.every=20, refit.window="moving")
class(MSFT.garch11.roll)
#plot(MSFT.garch11.roll)
# VaR plot
plot(MSFT.garch11.roll, which=4)
# Coef plot`
plot(MSFT.garch11.roll, which=5)
# show backtesting report
?report 
report(MSFT.garch11.roll, type="VaR")
report(MSFT.garch11.roll, type="fpm")
