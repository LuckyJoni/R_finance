#============================================================================
#   PORTFOLIO MANAGEMENT: DOW JONES CASE STUDY
#============================================================================

# Chapter 19

install.packages("corpcor", , repos="http://cran.us.r-project.org")

# no more available fro R version > 3
# install.packages("fEcofin", , repos="http://cran.us.r-project.org")
# library(fEcofin)


library(fPortfolio)
library(corpcor)

DowJones=read.table(file = "C:/DowJones.txt", header=TRUE)
djiData = as.timeSeries(DowJones30)
djiData.ret <- 100 * returns(djiData)
colnames(djiData)

c(start(djiData), end(djiData))


par(mfrow = c(1, 1), ask = TRUE)
for (i in 1:3) plot(djiData.ret[, (10 * i - 9):(10 * i)])
for (i in 1:3) plot(djiData[, (10 * i - 9):(10 * i)])
assetsCorImagePlot(djiData.ret)
plot(assetsSelect(djiData.ret))
assetsCorEigenPlot(djiData.ret)

frontier <- portfolioFrontier(djiData.ret)
tailoredFrontierPlot(frontier)
weightsPlot(frontier)

selection <- assetsSelect(djiData.ret, method = "kmeans")
cluster <- selection$cluster
cluster[cluster == 1]
cluster[cluster == 2]
cluster[cluster == 3]
cluster[cluster == 4]
cluster[cluster == 5]

constraints <- c(
'maxsumW[c("AXP","T","C","KO","XOM","GE","JPM","JNJ","MCD","MRK","MO","PG","SBC","DIS")] = 0.30',
'maxsumW[c("HD","WMT")] = 0.30',
'maxsumW[c("HWP","IBM")] = 0.30',
'maxsumW[c("INTC","MSFT")] = 0.30',
'maxsumW[c("AA","BA","CAT","DD","EK","GM","HON","IP","MMM","UTX")] = 0.30')


djiSpec <- portfolioSpec()
setNFrontierPoints(djiSpec) <- 25
setEstimator(djiSpec) <- "shrinkEstimator"
djiFrontier <- portfolioFrontier(djiData.ret, djiSpec)
col = seqPalette(30, "YlOrRd")
weightsPlot(djiFrontier, col = col)
