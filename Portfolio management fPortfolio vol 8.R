#============================================================================
#   PORTFOLIO MANAGEMENT: EXAMPLES USING THE fPortfolio PACKAGE
#============================================================================
#
# ATTENTION: IF YOU RUN THIS PROGRAM IN RSTUDIO, YOU WILL HAVE TO DECREASE THE MARGINS OF THE 
# PLOTS TO MAKE THEM VISIBLE (for example par(mar=c(1,1,1,1)))
# IN THIS CASE SUGGEST TO YOU, TO RUN THE CODE DIRECTLY IN R


# Chapter 16
install.packages("fPortfolio", repos="http://cran.us.r-project.org")

library(fPortfolio)
mvSpec <- portfolioSpec()
print(mvSpec)


# Chapter 17
# 17.1

colnames(LPP2005.RET)
lppData <- 100 * LPP2005.RET[, 1:6]
ewSpec <- portfolioSpec()
nAssets <- ncol(lppData)
setWeights(ewSpec) <- rep(1/nAssets, times = nAssets)

ewPortfolio <- feasiblePortfolio(
data = lppData,
spec = ewSpec,
constraints = "LongOnly")
print(ewPortfolio)


col <- divPalette(ncol(lppData), "RdBu")
weightsPie(ewPortfolio, radius = 0.7, col = col)
mtext(text = "Equally Weighted MV Portfolio", side = 3, line = 1.5, font = 2, cex = 0.7, adj = 0)
weightedReturnsPie(ewPortfolio, radius = 0.7, col = col)
mtext(text = "Equally Weighted MV Portfolio", side = 3, line = 1.5,font = 2, cex = 0.7, adj = 0)
covRiskBudgetsPie(ewPortfolio, radius = 0.7, col = col)
mtext(text = "Equally Weighted MV Portfolio", side = 3, line = 1.5,font = 2, cex = 0.7, adj = 0)

# 17.2

minriskSpec <- portfolioSpec()
targetReturn <- getTargetReturn(ewPortfolio@portfolio)["mean"]
setTargetReturn(minriskSpec) <- targetReturn

minriskPortfolio <- efficientPortfolio(
data = lppData,
spec = minriskSpec,
constraints = "LongOnly")
print(minriskPortfolio)

#17.3

globminSpec <- portfolioSpec()
globminPortfolio <- minvariancePortfolio(
data = lppData,
spec = globminSpec,
constraints = "LongOnly")
print(globminPortfolio)


col <- seqPalette(ncol(lppData), "YlGn")
weightsPie(globminPortfolio, box = FALSE, col = col)
mtext(text = "Global Minimum Variance MV Portfolio", side = 3, line = 1.5, font = 2, cex = 0.7, adj = 0)
weightedReturnsPie(globminPortfolio, box = FALSE, col = col)
mtext(text = "Global Minimum Variance MV Portfolio", side = 3, line = 1.5, font = 2, cex = 0.7, adj = 0)
covRiskBudgetsPie(globminPortfolio, box = FALSE, col = col)
mtext(text = "Global Minimum Variance MV Portfolio", side = 3, line = 1.5, font = 2, cex = 0.7, adj = 0)

#17.4

tgSpec <- portfolioSpec()
setRiskFreeRate(tgSpec) <- 0

tgPortfolio <- tangencyPortfolio(
data = lppData,
spec = tgSpec,
constraints = "LongOnly")
print(tgPortfolio)

col <- seqPalette(ncol(lppData), "BuPu")
weightsPie(tgPortfolio, box = FALSE, col = col)
mtext(text = "Tangency MV Portfolio", side = 3, line = 1.5, font = 2, cex = 0.7, adj = 0)
weightedReturnsPie(tgPortfolio, box = FALSE, col = col)
mtext(text = "Tangency MV Portfolio", side = 3, line = 1.5,font = 2, cex = 0.7, adj = 0)
covRiskBudgetsPie(tgPortfolio, box = FALSE, col = col)
mtext(text = "Tangency MV Portfolio", side = 3, line = 1.5,font = 2, cex = 0.7, adj = 0)


# 17.5

col = rampPalette(ncol(lppData), "purple2green")
weights <- 100 * as.vector(getWeights(tgPortfolio))
names <- as.vector(getUnits(tgPortfolio))
barplot(height = weights, names.arg = names,horiz = TRUE, las = 1, col = col)
title(main = "Weights of Long-Only Tangency Portfolio",xlab = "Weights %")


# 18.1

lppData <- 100*LPP2005.RET[, 1:6]
colnames(lppData)
lppSpec <- portfolioSpec()
setNFrontierPoints(lppSpec) <- 5
longFrontier <- portfolioFrontier(lppData, lppSpec)
print(longFrontier)

longFrontier <- portfolioFrontier(lppData)
plot(longFrontier)

args(tailoredFrontierPlot)
#tailoredFrontierPlot()

setNFrontierPoints(lppSpec) <- 25
longFrontier <- portfolioFrontier(lppData, lppSpec)
tailoredFrontierPlot(object = longFrontier, mText = "MV Portfolio - LongOnly Constraints",risk = "Cov")

weightsPlot(longFrontier)
text <- "Mean-Variance Portfolio - Long Only Constraints"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
weightedReturnsPlot(longFrontier)
covRiskBudgetsPlot(longFrontier)


#18.2

set.seed(1953)
frontierPlot(object = longFrontier, pch = 19, xlim = c(0.05, 0.85), cex = 0.5)

par(mfrow = c(1, 1))
set.seed(1953)
frontierPlot(object = longFrontier, pch = 19, xlim = c(0.05,0.85), cex = 0.5)
monteCarloPoints(object = longFrontier, mcSteps = 1000, pch = 19,cex = 0.5)
twoAssetsLines(object = longFrontier, col = "orange", lwd = 2)
frontier <- frontierPoints(object = longFrontier)
lines(frontier, col = "red", lwd = 2)


#18.3
shortSpec <- portfolioSpec()
setNFrontierPoints(shortSpec) <- 5
setSolver(shortSpec) <- "solveRshortExact"
shortFrontier <- portfolioFrontier(
data = lppData,
spec = shortSpec,
constraints = "Short")
print(shortFrontier)

setNFrontierPoints(shortSpec) <- 20
shortFrontier <- portfolioFrontier(data = lppData, spec = shortSpec,constraints = "Short")
tailoredFrontierPlot(object = shortFrontier, mText = "MV Portfolio - Short Constraints", risk = "Cov")
weightsPlot(shortFrontier)
text <- "MV Portfolio - Short Constrained Portfolio"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
weightedReturnsPlot(shortFrontier)
covRiskBudgetsPlot(shortFrontier)


#18.4

boxSpec <- portfolioSpec()
setNFrontierPoints(boxSpec) <- 15
boxConstraints <- c(
"minW[1:6]=0.01",
"maxW[1:6]=0.5")
boxFrontier <- portfolioFrontier(
data = lppData,
spec = boxSpec,
constraints = boxConstraints)
print(boxFrontier)

setNFrontierPoints(boxSpec) <- 25
boxFrontier <- portfolioFrontier(data = lppData, spec = boxSpec,constraints = boxConstraints)
tailoredFrontierPlot(object = boxFrontier, mText = "MV Portfolio - Box Constraints", risk = "Cov")
weightsPlot(boxFrontier)
text <- "MV Portfolio - Box Constrained Portfolio"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
weightedReturnsPlot(boxFrontier)
covRiskBudgetsPlot(boxFrontier)


#18.5

groupSpec <- portfolioSpec()
setNFrontierPoints(groupSpec) <- 7
groupConstraints <- c("minsumW[c(1,4)]=0.3","maxsumW[c(2,5)]=0.5")
groupFrontier <- portfolioFrontier(
data = lppData,
spec = groupSpec,
constraints = groupConstraints)
print(groupFrontier)


groupSpec <- portfolioSpec()
setNFrontierPoints(groupSpec) <- 25
groupFrontier <- portfolioFrontier(data = lppData, spec = groupSpec, constraints = groupConstraints)
tailoredFrontierPlot(object = groupFrontier, mText = "MV Portfolio - Group Constraints", risk = "Cov")

weightsPlot(groupFrontier)
text <- "MV Portfolio - Group Constrained Portfolio"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
weightedReturnsPlot(groupFrontier)
covRiskBudgetsPlot(groupFrontier)

#18.6

boxgroupSpec <- portfolioSpec()
setNFrontierPoints(boxgroupSpec) <- 15
boxgroupConstraints <- c(boxConstraints, groupConstraints)
boxgroupFrontier <- portfolioFrontier(
data = lppData,
spec = boxgroupSpec,
constraints = boxgroupConstraints)
print(boxgroupFrontier)


boxgroupSpec <- portfolioSpec()
setNFrontierPoints(boxgroupSpec) <- 25
boxgroupFrontier <- portfolioFrontier(
data = lppData,
spec = boxgroupSpec,
constraints = boxgroupConstraints)
tailoredFrontierPlot(
object = boxgroupFrontier,
mText = "MV Portfolio - Box/Group Constraints",
risk = "Cov")
weightsPlot(boxgroupFrontier)
text <- "MV Portfolio - Box/Group Constrained Portfolio"
mtext(text, side = 3, line = 3, font = 2, cex = 0.9)
weightedReturnsPlot(boxgroupFrontier)
covRiskBudgetsPlot(boxgroupFrontier)


#18.7

frontierPlot(longFrontier, auto = TRUE)
frontierPlot(longFrontier, return = "mean", risk = "Cov", auto = FALSE)
frontierPlot(longFrontier, return = "mean", risk = "CVaR",auto = FALSE)
frontierPlot(longFrontier, return = "mean", risk = "VaR",auto = FALSE)

