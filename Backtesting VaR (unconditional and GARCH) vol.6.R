#============================================================================
#                     BACKTESTING VaR MODELS
#============================================================================
# Load libraries
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
MSFT.GSPC.ret = merge(MSFT.ret,GSPC.ret)

# Plot returns
plot(MSFT.ret)
plot(GSPC.ret)


#
# Backtesting unconditional VaR models
# normal VaR, HS and modified HS
#

# Set up estimation window and testing window
n.obs = nrow(MSFT.ret)
w.e = 1000		 # estimation window
w.t = n.obs - w.e  # test window
alpha = 0.99

# Loop over testing sample, compute VaR and record hit rates
backTestVaR <- function(x, p = 0.95) {
  normal.VaR = as.numeric(VaR(x, p=p, method="gaussian"))
  historical.VaR = as.numeric(VaR(x, p=p, method="historical"))
  modified.VaR = as.numeric(VaR(x, p=p, method="modified"))
  ans = c(normal.VaR, historical.VaR, modified.VaR)
  names(ans) = c("Normal", "HS", "Modified")
  return(ans)
}

# Rolling 1-step ahead estimates of VaR
# we use 'rollapply', which is a zoo a generic zoo function for applying
# a function to rolling margins of an array.
VaR.results = rollapply(as.zoo(MSFT.ret), width=w.e,
                        FUN = backTestVaR, p=0.99, by.column = FALSE,
                        align = "right")
VaR.results = lag(VaR.results, k=-1)
chart.TimeSeries(merge(MSFT.ret, VaR.results), legend.loc="topright")


violations.mat = matrix(0, 3, 5)
rownames(violations.mat) = c("Normal", "HS", "Modified")
colnames(violations.mat) = c("En1", "n1", "1-alpha", "Percent", "VR")
violations.mat[, "En1"] = (1-alpha)*w.t
violations.mat[, "1-alpha"] = 1 - alpha

# Show Normal VaR violations
# we use the 'index' zoo function which is a generic function
# for extracting the index of an object and replacing it.
normalVaR.violations = as.zoo(MSFT.ret[index(VaR.results), ]) < VaR.results[, "Normal"]
violation.dates = index(normalVaR.violations[which(normalVaR.violations)])

# Plot violations
plot(as.zoo(MSFT.ret[index(VaR.results),]), col="blue", ylab="Return")
abline(h=0)
lines(VaR.results[, "Normal"], col="black", lwd=2)
lines(as.zoo(MSFT.ret[violation.dates,]), type="p", pch=16, col="red", lwd=2)

for(i in colnames(VaR.results)) {
  VaR.violations = as.zoo(MSFT.ret[index(VaR.results), ]) < VaR.results[, i]
  violations.mat[i, "n1"] = sum(VaR.violations)
  violations.mat[i, "Percent"] = sum(VaR.violations)/w.t
  violations.mat[i, "VR"] = violations.mat[i, "n1"]/violations.mat[i, "En1"]
}
violations.mat

?VaRTest
VaR.test = VaRTest(1-alpha,
                   actual=coredata(MSFT.ret[index(VaR.results),]),
                   VaR=coredata(VaR.results[,"Normal"]))
names(VaR.test)
# LR test for correct number of exceedances
VaR.test[1:7]

# LR tests for independence of exceedances
VaR.test[8:12]


# Backtest VaR but re-fit every 20 obvs (e.g. every month)
VaR.results.20 = rollapply(as.zoo(MSFT.ret), width=w.e, by = 20,
                        FUN = backTestVaR, p=0.99, by.column = FALSE,
                        align = "right")
chart.TimeSeries(merge(MSFT.ret, VaR.results, fill=na.locf), legend.loc="topright")




# Rolling GARCH(1,1) with VaR violations



MSFT.garch11.roll = ugarchroll(garch11.spec, MSFT.ret, n.ahead=1,
                               forecast.length = w.t,
                               refit.every=20, refit.window="moving",
                               calculate.VaR=TRUE, VaR.alpha=0.01)
class(MSFT.garch11.roll)

#plot(MSFT.garch11.roll)
# VaR plot
plot(MSFT.garch11.roll, which=4)
# Coef plot`
plot(MSFT.garch11.roll, which=5)
# show backtesting report
?report
report.msft = report(MSFT.garch11.roll, type="VaR")
report(MSFT.garch11.roll, type="fpm")

#Alternatively you can see the same thing in the following way:
xx = VaRTest(alpha=0.01, actual=MSFT.garch11.roll@forecast$VaR[,2],
             VaR=MSFT.garch11.roll@forecast$VaR[,1])
xx 

# As you can see GARCH estimates are not much better than the previous
# three unconditional methods
# HOMEWORK: do the same thing, but with different GARCH models
