library(crypto)
# Use Market capitalization data for Lykke from coinmarketcap.com
#  for the 2016 fiscal year
data.lyk <- read.csv("D:/Dean/LEZIONI/Basic R for finance ICEF/Lectures_notes_4_lectures/4_fourth_lecture/Credit_risk2/Lykke_market_cap.csv", sep=";")
# select data
data.select <-data.lyk$Market.cap

library(rootSolve)

# INPUTS of the MERTON’s MODEL
  # I can consider two options for the Average Equity in millions CHF
  # in December 2016:
  # a) Use Lykke Report Business 2016
  # equity=1
  # b) Use average market capitalization + USDCHF rate in December 2016
  equity=mean(data.select/1000000)*1.01

  # sigma_E is the annualized per cent standard deviation of returns
  # --> I follow Bharath and Shumway (2008, RFS, p.1351, section 3))
  sigma_eq=sqrt(252)*sd(diff(log(data.select/1000000)))

  exp_date=1     # Expiration date debt in years
  face_value=0.9 # Debt face value in millions CHF(from Lykke Report Business2016)
  rf_rate=0.0109 # Average 1-year mortgage rate in Switzerland in December 2016
                 # https://en.comparis.ch/hypotheken/zinssatz/zinsentwicklung
# Starting values for the two parameters V_t and sigma
  asset_value=equity+face_value
  sigma_asset=sigma_eq
  startval=c(asset_value,sigma_asset)


## 1) FIRST METHOD TO SOLVE THE NONLINEAR SYSTEM OF THE MERTON'S MODEL
# Function with the two-equation nonlinear system
merton.model = function(sol){
  #d1 and d2 from Merton’s model
  d1=(log(sol[1]/face_value)+
        (rf_rate+0.5*sol[2]^2)*exp_date)/(sol[2]*sqrt(exp_date))
  d2=(log(sol[1]/face_value)+
        (rf_rate-0.5*sol[2]^2)*exp_date)/(sol[2]*sqrt(exp_date))
  # Nonlinear equations to be solved
  f1=sol[1]*pnorm(d1)-face_value*exp(-rf_rate*exp_date)*pnorm(d2) - equity
  f2=pnorm(d1)*sol[2]*sol[1] - sigma_eq*equity
  c(f1 = f1, f2 = f2)
}
# Solve the nonlinear system using rootSolve::multiroot
sol = rootSolve::multiroot(merton.model, startval)$root

#Estimated Asset Value at time t (V_t):
sol[1]
#[1] 2.345131
# Estimated Asset volatility (sigma)
sol[2]
#[1] 0.5289149
# Estimated (Risk-Neutral) Default probability
d2=(log(sol[1]/face_value)+(rf_rate-0.5*sol[2]^2)*exp_date)/(sol[2]*sqrt(exp_date))
RN_def_prob = pnorm(-d2)
# Estimated Expected Default Frequency (EDF)
mu = mean(diff(log(data.select)))
EDF= pnorm(- (log(sol[1]/face_value)+
                (mu-0.5*sol[2]^2)*exp_date)/(sol[2]*sqrt(exp_date)) )
c(RN_def_prob, EDF)
#[1] 0.05857582 0.06036874


## 2) SECOND METHOD TO SOLVE THE NONLINEAR SYSTEM OF THE MERTON'S MODEL
# Function with the two-equation nonlinear system expressed as a single eq.
merton.model2 = function(sol){
  #d1 and d2 from Merton’s model
  d1=(log(sol[1]/face_value)+
        (rf_rate+0.5*sol[2]^2)*exp_date)/(sol[2]*sqrt(exp_date))
  d2=(log(sol[1]/face_value)+
        (rf_rate-0.5*sol[2]^2)*exp_date)/(sol[2]*sqrt(exp_date))
  # Nonlinear equations to be solved
  f1=sol[1]*pnorm(d1)-face_value*exp(-rf_rate*exp_date)*pnorm(d2) - equity
  f2=pnorm(d1)*sol[2]*sol[1] - sigma_eq*equity
  return(f1^2+f2^2)
}
# Solve the nonlinear equation using stats::optim
sol2=stats::optim(startval,merton.model2)$par
#Estimated Asset Value at time t (V_t):
sol2[1]
#[1] 2.345159
# Estimated Asset volatility (sigma)
sol2[2]
#[1] 0.52889


# Use available data for Lykke from coinmarketcap.com
data.lyk <-  read.csv("D:/Dean/LEZIONI/Basic R for finance ICEF/Lectures_notes_4_lectures/4_fourth_lecture/Credit_risk2/lkk-usd-max.csv")
data.lyk$date<- as.Date(data.lyk$snapped_at)
# extract close prices
data.select <-data.lyk$price

# ZPP with the Naive Constant Variance model
zpp.NCV<-bitcoinFinance::zpp_norm(data.select)
zpp.NCV
#[1] 0.7143673
# ZPP with Normal GARCH(1,1) model
zpp.norm.garch11<-bitcoinFinance::zpp_norm_garch11(data.select)
zpp.norm.garch11
#[1] 0.1351545
# ZPP with Student's T GARCH(1,1) model
set.seed(1)
zpp.t.garch<-bitcoinFinance::zpp_general(data.select, distribution.model = "std")
zpp.t.garch
#[1] 0.04975

system.time(bitcoinFinance::zpp_norm(data.select))
##    user  system elapsed
##       0       0       0
system.time(bitcoinFinance::zpp_norm_garch11(data.select))
##    user  system elapsed
##    0.21    0.00    0.20 
system.time(bitcoinFinance::zpp_general(data.select,
              distribution.model="std", scenarios=5000))
##    user  system elapsed
##    1.77    0.01    1.82


zpp.roll.500.days<-zoo::rollapply( data.select, width=500,
        FUN=bitcoinFinance::zpp_general, distribution.model="std",
        scenarios =5000, by.column=F, align='right')
zpp.roll.500.days<-xts::xts( zpp.roll.500.days,
        order.by = data.lyk$date[500:NROW(data.lyk$date)] )
plot(zpp.roll.500.days)




# Download the daily time series data of the CoinDesk Bitcoin Price Index (BPI)
xbp_data<-bitcoinFinance::coindesk_download(start="2011-01-01",end="2018-06-20")
zpp.roll.1000.days<-zoo::rollapply( as.numeric(xbp_data), width=1000,
                   FUN=bitcoinFinance::zpp_general, distribution.model="std",
                   scenarios =2500, by.column=F, align='right')
zpp.roll.1000.days<-xts::xts( zpp.roll.1000.days,
                           order.by = zoo::index(xbp_data)[1000:NROW(xbp_data)] )
plot(zpp.roll.1000.days)


zpp.roll.1000.days.norm<-zoo::rollapply( as.numeric(xbp_data), width=1000,
                   FUN=bitcoinFinance::zpp_norm, by.column=F, align='right')
zpp.roll.1000.days.norm<-xts::xts(zpp.roll.1000.days.norm,
                           order.by = zoo::index(xbp_data)[1000:NROW(xbp_data)] )
plot(zpp.roll.1000.days.norm)



# Use available data for noirshares from coinmarketcap.com
data.ns <-  read.csv("D:/Dean/LEZIONI/Basic R for finance ICEF/Lectures_notes_4_lectures/4_fourth_lecture/Credit_risk2/noirshares.csv", sep=";")
# extract close prices
data.select <-data.ns$Close
data.select <- xts::xts(data.select, order.by = lubridate::mdy(data.ns$Date))
# Compute ZPP
zpp.roll.100.days.norm<-zoo::rollapply( as.numeric(data.select), width=100,
                   FUN=bitcoinFinance::zpp_norm, by.column=F, align='right')
# Impose the upper bound of 1: due to parameters biases, the probability is
# not granted to be bounded between 0 and 1.
zpp.roll.100.days.norm<-pmin(zpp.roll.100.days.norm,1)
zpp.roll.100.days.norm<-xts::xts(zpp.roll.100.days.norm,
              order.by = zoo::index(data.select[100:NROW(data.ns)]) )
plot(zpp.roll.100.days.norm)



# Use available data for noirshares from coinmarketcap.com
data.BBQ <-  read.csv("D:/Dean/LEZIONI/Basic R for finance ICEF/Lectures_notes_4_lectures/4_fourth_lecture/Credit_risk2/BBQCoin.csv", sep=";")
# extract close prices
data.select <-data.BBQ$Close
data.select <- xts::xts(data.select, order.by = lubridate::mdy(data.BBQ$Date))
# Compute ZPP
zpp.roll.100.days.norm<-zoo::rollapply( as.numeric(data.select), width=100,
                   FUN=bitcoinFinance::zpp_norm, by.column=F, align='right')
# Impose the upper bound of 1: due to parameters biases, the probability is
# not granted to be bounded between 0 and 1.
zpp.roll.100.days.norm<-pmin(zpp.roll.100.days.norm,1)
zpp.roll.100.days.norm<-xts::xts(zpp.roll.100.days.norm,
              order.by = zoo::index(data.select[100:NROW(data.BBQ)]) )
plot(zpp.roll.100.days.norm)


