#============================================================================
#                            UNIVARIATE GARCH MODELS PART 2 OF 2
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
MSFT.GSPC.ret = cbind(MSFT.ret,GSPC.ret)

# Plot returns
plot(MSFT.ret)
plot(GSPC.ret)


# GARCH(1,1) MODEL
garch11.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                          mean.model = list(armaOrder=c(0,0)))
MSFT.garch11.fit = ugarchfit(spec=garch11.spec, data=MSFT.ret,
                             solver.control=list(trace = 1))
GSPC.garch11.fit = ugarchfit(spec=garch11.spec, data=GSPC.ret,
                             solver.control=list(trace = 1))
MSFT.garch11.fit
GSPC.garch11.fit

# Engle-Ng sign bias test
signbias(MSFT.garch11.fit)
signbias(GSPC.garch11.fit)


# ASYMMETRIC GARCH MODELS
# Nelson's egarch model
egarch11.spec = ugarchspec(variance.model=list(model="eGARCH",
                                               garchOrder=c(1,1)),
                           mean.model=list(armaOrder=c(0,0)))
MSFT.egarch11.fit = ugarchfit(egarch11.spec, MSFT.ret)
MSFT.egarch11.fit

# GJR garch model
gjrgarch11.spec = ugarchspec(variance.model=list(model="gjrGARCH",
                                                 garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(0,0)))
MSFT.gjrgarch11.fit = ugarchfit(gjrgarch11.spec, MSFT.ret)
MSFT.gjrgarch11.fit

# Aparch models
aparch11.1.spec = ugarchspec(variance.model=list(model="apARCH",
                                                 garchOrder=c(1,1)),
                             mean.model=list(armaOrder=c(0,0)),
                             fixed.pars=list(delta=1))

MSFT.aparch11.1.fit = ugarchfit(aparch11.1.spec, MSFT.ret)
MSFT.aparch11.1.fit


nic.garch11 = newsimpact(MSFT.garch11.fit)
nic.egarch11 = newsimpact(MSFT.egarch11.fit)
nic.gjrgarch11 = newsimpact(MSFT.gjrgarch11.fit)
nic.aparch11.1 = newsimpact(MSFT.aparch11.1.fit)

# Compare information criteria
model.list = list(garch11 = MSFT.garch11.fit,
                  egarch11 = MSFT.egarch11.fit,
                  gjrgarch11 = MSFT.gjrgarch11.fit,
                  aparch11.1 = MSFT.aparch11.1.fit)
info.mat = sapply(model.list, infocriteria)
rownames(info.mat) = rownames(infocriteria(MSFT.garch11.fit))
info.mat

# Show news impact curve from estimated garch(1,1) and egarch(1,1)
par(mfrow=c(2,2))
plot(nic.garch11$zx, type="l", lwd=2, col="blue", main="GARCH(1,1)",
     nic.garch11$zy, ylab=nic.garch11$yexpr, xlab=nic.garch11$xexpr)
plot(nic.egarch11$zx, type="l", lwd=2, col="blue", main="EGARCH(1,1)",
     nic.egarch11$zy, ylab=nic.egarch11$yexpr, xlab=nic.egarch11$xexpr)
plot(nic.gjrgarch11$zx, type="l", lwd=2, col="blue", main="TGARCH(1,1)",
     nic.gjrgarch11$zy, ylab=nic.gjrgarch11$yexpr, xlab=nic.gjrgarch11$xexpr)
plot(nic.aparch11.1$zx, type="l", lwd=2, col="blue", main="APARCH(1,1,1)",
     nic.aparch11.1$zy, ylab=nic.aparch11.1$yexpr, xlab=nic.aparch11.1$xexpr)
par(mfrow=c(1,1))



# GARCH WITH NON-NORMAL ERRORS
# recall normal GARCH(1,1)
# examine standardized residuals
MSFT.garch11.zt = residuals(MSFT.garch11.fit)/sigma(MSFT.garch11.fit)
MSFT.garch11.zt = xts(MSFT.garch11.zt, order.by=index(MSFT.ret))
qqPlot(coredata(MSFT.garch11.zt))
plot(MSFT.garch11.fit, which=8)
plot(MSFT.garch11.fit, which=9)


# GARCH(1,1)  with Student-t errors
garch11.t.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "std")
MSFT.garch11.t.fit = ugarchfit(spec=garch11.t.spec, data=MSFT.ret)
MSFT.garch11.t.fit
plot(MSFT.garch11.t.fit, which=9)

aparch11.1.t.spec = ugarchspec(variance.model = list(model="apARCH",
                                                     garchOrder=c(1,1)),
                               mean.model = list(armaOrder=c(0,0)),
                               distribution.model = "std",
                               fixed.pars=list(delta=1))
MSFT.aparch11.1.t.fit = ugarchfit(spec=aparch11.1.t.spec, data=MSFT.ret)
MSFT.aparch11.1.t.fit

# Fit skewed t
garch11.st.spec = ugarchspec(variance.model = list(garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0)),
                            distribution.model = "sstd")
MSFT.garch11.st.fit = ugarchfit(spec=garch11.st.spec, data=MSFT.ret)
MSFT.garch11.st.fit
plot(MSFT.garch11.st.fit, which=9)


# PLOT FORECASTS FROM COMPETING MODELS
MSFT.garch11.fcst = ugarchforecast(MSFT.garch11.fit, n.ahead=250)
MSFT.garch11.t.fcst = ugarchforecast(MSFT.garch11.t.fit, n.ahead=250)
MSFT.aparch11.1.fcst = ugarchforecast(MSFT.aparch11.1.fit, n.ahead=250)
MSFT.aparch11.1.t.fcst = ugarchforecast(MSFT.aparch11.1.t.fit, n.ahead=250)

# Extract volatility forecasts
MSFT.garch11.sigma = as.data.frame(MSFT.garch11.fcst@forecast$sigmaFor)
MSFT.garch11.t.sigma = as.data.frame(MSFT.garch11.t.fcst@forecast$sigmaFor)
MSFT.aparch11.1.sigma = as.data.frame(MSFT.aparch11.1.fcst@forecast$sigmaFor)
MSFT.aparch11.1.t.sigma = as.data.frame(MSFT.aparch11.1.t.fcst@forecast$sigmaFor)

ymax = max(MSFT.garch11.sigma,MSFT.garch11.t.sigma,MSFT.aparch11.1.sigma, MSFT.aparch11.1.t.sigma)
ymin = min(MSFT.garch11.sigma,MSFT.garch11.t.sigma,MSFT.aparch11.1.sigma, MSFT.aparch11.1.t.sigma)

plot.ts(MSFT.garch11.sigma, main="Volatility Forecasts",
        ylim=c(ymin,ymax), col="black",
        lwd=2, ylab="sigma(t+h|t)", xlab="h")
lines(MSFT.garch11.t.sigma, col="blue", lwd=2)
lines(MSFT.aparch11.1.sigma, col="green", lwd=2)
lines(MSFT.aparch11.1.t.sigma, col="red", lwd=2)
legend(x="topleft", legend=c("GARCH-n", "GARCH-t", "APARCH-n", "APARCH-t"),
       col=c("black", "blue","green","red"), lwd=2, lty = "solid")

# Evaluate rolling forecasts
# re-fit models leaving 100 out-of-sample observations for forecast evaluation statistics
MSFT.garch11.fit = ugarchfit(spec=garch11.spec, data=MSFT.ret,
                             out.sample=100)
MSFT.garch11.t.fit = ugarchfit(spec=garch11.t.spec, data=MSFT.ret,
                               out.sample=100)
MSFT.aparch11.1.fit = ugarchfit(aparch11.1.spec, MSFT.ret,
                                out.sample=100)
MSFT.aparch11.1.t.fit = ugarchfit(spec=aparch11.1.t.spec, data=MSFT.ret,
                                  out.sample=100)

# Compare persistence and unconditional variance
c.mat = matrix(0, 4, 2)
colnames(c.mat) = c("Persistence", "E[sig(t)]")
rownames(c.mat) = c("GARCH-n", "GARCH-t", "APARCH-n","APARCH-t")
c.mat["GARCH-n","Persistence"] = persistence(MSFT.garch11.fit)
c.mat["GARCH-t","Persistence"] = persistence(MSFT.garch11.t.fit)
c.mat["APARCH-n","Persistence"] = persistence(MSFT.aparch11.1.fit)
c.mat["APARCH-t","Persistence"] = persistence(MSFT.aparch11.1.t.fit)

c.mat["GARCH-n","E[sig(t)]"] = sqrt(uncvariance(MSFT.garch11.fit))
c.mat["GARCH-t","E[sig(t)]"] = sqrt(uncvariance(MSFT.garch11.t.fit))
c.mat["APARCH-n","E[sig(t)]"] = sqrt(uncvariance(MSFT.aparch11.1.fit))
c.mat["APARCH-t","E[sig(t)]"] = sqrt(uncvariance(MSFT.aparch11.1.t.fit))

c.mat

# Compute 100 1-step ahead rolling forecasts
MSFT.garch11.fcst = ugarchforecast(MSFT.garch11.fit, n.roll=100, n.ahead=1)
MSFT.garch11.t.fcst = ugarchforecast(MSFT.garch11.t.fit, n.roll=100, n.ahead=1)
MSFT.aparch11.1.fcst = ugarchforecast(MSFT.aparch11.1.fit, n.roll=100, n.ahead=1)
MSFT.aparch11.1.t.fcst = ugarchforecast(MSFT.aparch11.1.t.fit, n.roll=100, n.ahead=1)

# compute forecast evaluation statistics
fcst.list = list(garch11=MSFT.garch11.fcst,
                 garch11.t=MSFT.garch11.t.fcst,
                 aparch11.1=MSFT.aparch11.1.fcst,
                 aparch11.t.1=MSFT.aparch11.1.t.fcst)
fpm.mat = sapply(fcst.list, fpm)
fpm.mat 

