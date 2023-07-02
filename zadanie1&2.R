
if (!require("pacman")) install.packages("pacman")
pacman::p_load(RCurl,RcmdrMisc,openxlsx,readxl,data.table,quantmod,xts,DistributionUtils,nortest,Jmisc)

#####Zadanie1######
#####Symulacja GARCH(2,2) #####

theta <- c(0.1, 0.1, 0.06, 0.06, 0.03) #omega, alfa1, beta1, alfa2, beta2

n<-1000
set.seed(123)
e <- rnorm(n)  
xi <- double(n)
sigma2 <- double(n)
xi[1] <- rnorm(1, sd = sqrt(theta[1]/(1.0-theta[2]-theta[3]-theta[4]-theta[5])))
sigma2[1] = xi[1]^2 / rnorm(1)^2

for(i in 3:n)
{
  sigma2[i] = theta[1] + theta[2]*xi[i-1]^2 + theta[3]*sigma2[i-1] + theta[4]*xi[i-2]^2 + theta[5]*sigma2[i-2]
  xi[i] <- e[i] * sqrt(sigma2[i])
}

ts.plot(xi)
hist(xi, breaks = "FD", col = "#31766B", main = "Histogram xi")
qqplot(qnorm(ppoints(n)), col="#4B8465",xlab="", xi, main = "Wykres kwantylowy xi")


#####testowanie stacjonarności i występowania efektu ARCH#####

library(tseries)

#stacjonarność szeregu
adf.test(xi)

# testujemy efekty ARCH - Ljungi-Boxa dla kwadratów
Box.test(xi^2, type="Ljung-Box")
Box.test(xi^2, lag=2, type="Ljung-Box")
Box.test(xi^2, lag=5, type="Ljung-Box")
Box.test(xi^2, lag=10, type="Ljung-Box")

#####Estymacja MNW i funkcja wiarygodności#####
theta
#estymacja MNW
# Log-likelihood GARCH(2,2)
  
n=length(xi)
logL_GARCH <- function(theta) {
  xi_0<-xi[1]
  sigma2_0<-xi[1]^2
  
  sigma2[1] <- theta[1]+theta[2]*xi_0^2+theta[3]*sigma2_0
  sum <- 0.5 * log(sigma2[1]) + xi[1]^2 / (2 * sigma2[1])
    
  sigma2[2] <- theta[1] + theta[2] * (xi[1])^2 + theta[3] * sigma2[1]+theta[4]*(xi_0)^2+theta[5]*sigma2_0
  sum <- sum + 0.5 * log(sigma2[2]) + xi[2]^2 / (2 * sigma2[2])

  for (i in 3:n) {
    sigma2[i] <- theta[1] + theta[2] * xi[i-1]^2 + theta[3] * sigma2[i-1] + theta[4] * xi[i-2]^2 + theta[5] * sigma2[i-2]
    sum <- sum + 0.5 * log(sigma2[i]) + xi[i]^2 / (2 * sigma2[i])
  }
  loglike<--n/2*log(2*pi)-sum
  return(c(loglike))
}


# Liczymy funkcję log-wiarogodności
logL_GARCH(c(0.1, 0.1, 0.06, 0.06, 0.03))

x1<-seq(0,2,by=0.01)
logL_GARCH1 <- x1
for(i in 1:length(x1)){
  logL_GARCH1[i]<-logL_GARCH(c(x1[i], 0.1, 0.06, 0.06, 0.03))
}
x2<-seq(0,2,by=0.01)
logL_GARCH2 <- x2
for(i in 1:length(x2)){
  logL_GARCH2[i]<-logL_GARCH(c(0.1, x2[1], 0.06, 0.06, 0.03))
}
x3<-seq(0,2,by=0.01)
logL_GARCH3 <- x3
for(i in 1:length(x3)){
  logL_GARCH3[i]<-logL_GARCH(c(0.1, 0.1, x3[1], 0.06, 0.03))
}
x4<-seq(0,2,by=0.01)
logL_GARCH4 <- x4
for(i in 1:length(x4)){
  logL_GARCH4[i]<-logL_GARCH(c(0.1, 0.1, 0.06, x4[i], 0.03))
}
x5<-seq(0,2,by=0.01)
logL_GARCH5 <- x5
for(i in 1:length(x5)){
  logL_GARCH5[i]<-logL_GARCH(c(0.1, 0.1, 0.06, 0.06, x5[i]))
}
par(mfrow=c(1,3))
par(mar = c(2, 2, 2, 1))
plot(x1,logL_GARCH1, type="l", main="x1")
plot(x4,logL_GARCH4, type="l", main="x4")
plot(x5, logL_GARCH5, type="l", main="x5")

plot(x3,logL_GARCH3, type="l")
plot(x2,logL_GARCH2, type="l") #funkcja stała

#Numeryczna optymalizacja

hat_theta1 <- optim(par = c(0.25, 0.15, 0.25, 0.2, 0.1), fn =logL_GARCH, control = list(fnscale = -1))
hat_theta1$par

hat_theta2 <- optim(par = c(0.001, 0.001, 0.001, 0.001, 0.001), fn =logL_GARCH, control = list(fnscale = -1))
hat_theta2$par

hat_theta3 <- optim(par = c(0.001, 0.003, 0.002, 0.002, 0.0015), fn =logL_GARCH, control = list(fnscale = -1))
hat_theta3$par

theta

# zmieniamy metodę optymalizacji na "L-BFGS", aby dostać hessian
hat_theta <- optim(par = c(0.001, 0.003, 0.002, 0.002, 0.0015), method = "L-BFGS", fn = logL_GARCH, control = list(fnscale = -1),
                   lower = c(0,0,0,0,0), hessian=TRUE )

hat_theta$par

#####Ugarch#####

garch22_spec <- ugarchspec(
  variance.model = list(model = "sGARCH", garchOrder = c(2, 2)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "norm")

#Estymacja GARCH(2,2) dla danych symulowanych
modelfit <- ugarchfit(spec = garch22_spec, data = xi)
par(mar = c(4, 4, 2, 1))
plot(modelfit_2, which = "all")

#####wariancja warunkowa, zmienność warunkowa, prognozy######

layout(1)
# Wykresy wariancji i zmienności warunkowej GARCH(2,2)
plot(tail(sigma2,998), main = "Wariancja warunkowa GARCH(2,2), theta = (0.1, 0.1, 0.06, 0.06, 0.03)",
     type="line", ylab="",xlab="")
plot(sqrt(sigma2), main = "Zmienność warunkowa GARCH(2,2),  theta = (0.1, 0.1, 0.06, 0.06, 0.03)",
     type="line", ylab="",xlab="")

# Prognozy z modelu garch22 dla danych xi
prog_2 <- ugarchforecast(modelfit, data = xi, n.ahead = 10)
prog_2
plot(prog_2)


#####Zadanie2######

if (!require("pacman")) install.packages("pacman")
pacman::p_load(RCurl,RcmdrMisc,openxlsx,readxl,data.table,quantmod,xts,DistributionUtils,nortest,Jmisc)
#setwd - folder w którym znajduje się import danych
source('import_danych.R', encoding = 'UTF-8')
library(dplyr)

#####stopy zwrotu, wykresy######

log_stooq_wig20 <- xts(stooq_wig20$Close, order.by = as.Date(stooq_wig20$Date)) %>% periodReturn(., period = 'daily', type = 'log') %>%
  tail(., 1000) %>% ts(.)

log_stooq_snk <- xts(stooq_snk$Close, order.by = as.Date(stooq_snk$Date)) %>% periodReturn(., period = 'daily', type = 'log') %>%
  tail(., 1000) %>% ts(.)

log_stooq_mil <- xts(stooq_mil$Close, order.by = as.Date(stooq_mil$Date)) %>% periodReturn(., period = 'daily', type = 'log') %>%
  tail(., 1000) %>% ts(.)

#wykresy

ts.plot(log_stooq_wig20, col="#3A4A3F")
ts.plot(log_stooq_snk, col="#3A4A3F")
ts.plot(log_stooq_mil, col="#3A4A3F")
hist(log_stooq_wig20, breaks = "FD", col = "#2e8b57", main = "Histogram log_stooq_wig20")
hist(log_stooq_snk, breaks = "FD", col = "palegreen3", main = "Histogram loq_stooq_snk")
hist(log_stooq_mil, breaks = "FD", col = "#6D9B7C", main = "Histogram loq_stooq_mil")
qqplot(qnorm(ppoints(length(stooq_wig20))), xlab = "",log_stooq_wig20, col = "#2e8b57", main = "Wykres kwantylowy log_stooq_wig20")
qqplot(qnorm(ppoints(length(stooq_snk))),xlab = "", log_stooq_snk, col = "palegreen3", main = "Wykres kwantylowy loq_stooq_snk")
qqplot(qnorm(ppoints(length(stooq_mil))), xlab = "",log_stooq_mil, col = "#6D9B7C", main = "Wykres kwantylowy loq_stooq_mil")

#####testowanie stacjonarności i występowania efektu ARCH#####

  library(tseries)

  #stacjonarność szeregu
  adf.test(log_stooq_wig20)
  adf.test(log_stooq_snk)  
  adf.test(log_stooq_mil)
  
  #testujemy efekty ARCH - test Ljunga-Boxa dla kwadratów reszt
  Box.test(log_stooq_wig20^2, type="Ljung-Box")
  Box.test(log_stooq_wig20^2, lag=2, type="Ljung-Box")
  Box.test(log_stooq_wig20^2, lag=5, type="Ljung-Box")
  Box.test(log_stooq_wig20^2, lag=10, type="Ljung-Box")
  
  Box.test(log_stooq_snk^2, type="Ljung-Box")
  Box.test(log_stooq_snk^2, lag=2, type="Ljung-Box")
  Box.test(log_stooq_snk^2, lag=5, type="Ljung-Box")
  Box.test(log_stooq_snk^2, lag=10, type="Ljung-Box")
  
  Box.test(log_stooq_mil^2, type="Ljung-Box")
  Box.test(log_stooq_mil^2, lag=2, type="Ljung-Box")
  Box.test(log_stooq_mil^2, lag=5, type="Ljung-Box")
  Box.test(log_stooq_mil^2, lag=10, type="Ljung-Box")

#####AR-GARCH######
library(rugarch)

# Estymacja modelu AR-GARCH 
spec1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(1, 0)))
spec2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 2)), mean.model = list(armaOrder = c(1, 0)))
spec3<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), mean.model = list(armaOrder = c(2, 0)))

fit1_wig20 <- ugarchfit(spec1, data = log_stooq_wig20)
fit2_wig20 <- ugarchfit(spec2, data = log_stooq_wig20)
fit3_wig20 <- ugarchfit(spec3, data = log_stooq_wig20)
fit1_snk <- ugarchfit(spec1, data = log_stooq_snk)
fit2_snk <- ugarchfit(spec2, data = log_stooq_snk)
fit3_snk <- ugarchfit(spec3, data = log_stooq_snk)
fit1_mil <- ugarchfit(spec1, data = log_stooq_mil)
fit2_mil <- ugarchfit(spec2, data = log_stooq_mil)
fit3_mil <- ugarchfit(spec3, data = log_stooq_mil)

#####warunkowe wariancje, zmienności i prognozy######

modelfit_wig <- ugarchfit(spec = spec2, data = log_stooq_wig20)
modelfit_snk <- ugarchfit(spec = spec1, data = log_stooq_snk)
modelfit_mil <- ugarchfit(spec = spec1, data = log_stooq_mil)

# Wykresy wariancji warunkowej i zmienności warunkowej
sigma2_wig<- (modelfit_wig@fit$sigma)^2
plot(sigma2_wig, main = "Wariancja warunkowa", type="line")
plot(sqrt(sigma2_wig), main = "Zmienność warunkowa", type="line")

sigma2_snk<- (modelfit_snk@fit$sigma)^2
plot(sigma2_snk, main = "Wariancja warunkowa", type="line")
plot(sqrt(sigma2_snk), main = "Zmienność warunkowa", type="line")

sigma2_mil<- (modelfit_mil@fit$sigma)^2
plot(sigma2_mil, main = "Wariancja warunkowa", type="line")
plot(sqrt(sigma2_mil), main = "Zmienność warunkowa", type="line")

# Prognozy z modelu
prog_wig <- ugarchforecast(modelfit_wig, data = sigma2_wig, n.ahead = 10)
prog_snk <- ugarchforecast(modelfit_snk, data = log_stooq_snk, n.ahead = 10)
prog_mil <- ugarchforecast(modelfit_mil, data = log_stooq_mil, n.ahead = 10)
plot(prog_wig) #Wybór : 1
plot(prog_snk)
plot(prog_mil)
