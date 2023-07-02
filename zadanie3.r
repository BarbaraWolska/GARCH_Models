##### Zadanie3 ####

rm(list=ls())
# biblioteki 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl,xts,zoo,rugarch,rmgarch,quantmod,MTS)
library(dplyr)
library(magrittr)

getSymbols('BUR.PA', src='yahoo',from=Sys.Date()-5*365-1)
getSymbols('^FCHI', src='yahoo',from=Sys.Date()-5*365-1)

y1 <- na.omit(`BUR.PA`) %>% xts(.$`BUR.PA`, order.by = as.Date(index(.))) %>% 
  .$`BUR.PA.Close` %>% set_colnames(.,'Burelle') %>% as.data.frame(.) %>% 
  mutate(log_return=c(NA, 100 * diff(log(Burelle)))) %>% na.omit(.) %>% tail(., 1000) %>% 
  rename(log_Burelle=log_return)

y2 <- na.omit(`FCHI`) %>% xts(.$`FCHI`, order.by = as.Date(index(.))) %>% 
  .$`FCHI.Close` %>% set_colnames(.,'FCHI') %>% as.data.frame(.) %>% 
  mutate(log_return=c(NA, 100 * diff(log(FCHI)))) %>% na.omit(.) %>% tail(., 1000) %>% 
  rename(log_FCHI=log_return)
y<-cbind(y1,y2)
daty <- as.Date(rownames(y1))

#Wykresy
ts.plot(y$log_Burelle,col="#3A4A3F", ylab="log_Burelle", xlab="")
ts.plot(y$log_FCHI,col="#3A4A3F", ylab="log_FCHI", xlab="")

##### 2. Estymacja odch. std, oraz korelacji z próby rolowanej ####

w <- 60

Rsd1  <- rollapply(y[,2], width=w, sd, align="right")
Rsd2  <- rollapply(y[,4], width=w, sd, align="right")
Rcor  <- rollapply(y, width=w, function(x) cor(x[,2],x[,4]), by.column=FALSE, align="right")
Rcov  <- rollapply(y, width=w, function(x) cov(x[,2],x[,4]), by.column=FALSE, align="right")

windows()
par(mfrow=c(2,2), cex=0.5, bty="l", pin=c(2,2), lwd=1.5)
plot(Rsd1, main="Estymator T=60 odch. std. dla log_Burelle", xlab="", ylab="", col = "#3A4A3F", type="line")
plot(Rsd2, main="Estymator T=60 odch. std. dla log_FCHI", xlab="", ylab="", col = "#3A4A3F", type="line")
plot(Rcor, main="Estymator T=60 korelacji dla log_Burelle vs log_FCHI", xlab="", ylab="", col = "#3A4A3F", type="line")
abline(h=0)
plot(Rcov, main="Estymator T=60 kowariancji dla log_Burelle vs log_FCHI", xlab="", ylab="", col = "#3A4A3F", type="line")
abline(h=0)

#########3. DCC-GARCH ###########
y<-cbind(y[,2],y[,4]) %>% set_colnames(.,c("log_Burelle","log_FCHI"))
spec1<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), mean.model = list(armaOrder = c(1, 0)))
spec2<-ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(2, 2)), mean.model = list(armaOrder = c(2, 0)))

dcc_spec1 <- dccspec(multispec(c( spec1, spec1)), dccOrder = c(1,1), model = "DCC",distribution = "mvnorm")
dcc_spec2 <- dccspec(multispec(c( spec1, spec1)), dccOrder = c(1,0), model = "DCC",distribution = "mvnorm")
dcc_spec3 <- dccspec(multispec(c( spec2, spec2)), dccOrder = c(1,1), model = "DCC",distribution = "mvnorm")
dcc_spec4 <- dccspec(multispec(c( spec2, spec2)), dccOrder = c(1,0), model = "DCC",distribution = "mvnorm")

fit1<-dccfit(dcc_spec1, y)
fit1
# macierz Sigmat, zmienności, kowariancje i korelacje

sd1  <- rcov(fit1)  
sd1  <- as.xts(sd1[1,1,]^0.5, daty)
sd2  <- rcov(fit1)
sd2  <- as.xts(sd2[2,2,]^0.5, daty)
cov  <- rcov(fit1) 
cov  <- as.xts(cov[1,2,], daty)
cor  <- rcor(fit1) 
cor  <- as.xts(cor[1,2,], order.by=daty)

T      <- length(daty)

windows()
par(mfrow=c(2,2), cex=0.5, bty="l", pin=c(2,2), lwd=1)
plot(sd1, main="DCC odch. std. dla log_Burelle", xlab="", ylab="")
plot(sd2, main="DCC odch. std. dla log_FCHI", xlab="", ylab="")
plot(cor, main="DCC korelacja dla log_Burelle vs log_FCHI", xlab="", ylab="")
abline(h=0)
plot(cov, main="DCC kowariancja dla log_Burelle vs log_FCHI", xlab="", ylab="")
abline(h=0)

#####4.PROGNOZY#####

library(forecast)
prog1 <- forecast(sd1, h = 50)
prog2 <- forecast(sd2, h = 50)
prog3 <- forecast(cov, h = 50)
prog4 <- forecast(cor, h = 50)

windows()
par(mfrow=c(2,2), cex=0.5, bty="l", pin=c(2,2), lwd=1)
plot(prog1, main="Prognoza warunkowego odchylenia std. Burelle", xlab="", ylab="")
plot(prog2, main="Prognoza warunkowego odchylenia std. FCHI", xlab="", ylab="")
plot(prog3, main="Prognoza warunkowej kowariancji", xlab="", ylab="")
abline(h=0)
plot(prog4, main="Prognoza warunkowej korelacji", xlab="", ylab="")
abline(h=0)
