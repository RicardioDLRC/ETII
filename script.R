
getwd()

datos = read.csv(
  "C:/Users/Ricardo/Documents/Ecología/ENES/6to sem/Ec_Teórica/Gerardo/ETII/HeV-survival.csv")

d.4 <- subset(datos, Temp == 4)
d.22 <- subset(datos, Temp == 22)
d.56 <- subset(datos, Temp == 56)

m4 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
          data = d.4, 
          start = list(rho = 1, kappa = 0.9),
          lower = c(0.00001, 0.1),
          upper = c(1, 1),
          algorithm =  "port")

summary(m4)

m22 <- nls(ln.S ~ -( rho * Time.h )^ (kappa), 
           data = d.22, 
           start = list(rho = 1, kappa = 0.9),
           lower = c(0.0001, 0.1),
           upper = c(1, 1),
           algorithm =  "port")

summary(m22)

m56 <- nls(ln.S ~ -( rho * Time )^ (kappa), 
           data = d.56, 
           start = list(rho = 1, kappa = 0.9),
           lower = c(0.1, 0.1),
           upper = c(50, 1.5),
           algorithm =  "port")

summary(m56)

#Identificando efecto de temperatura sobre rho y kappa

par.estim <- data.frame(rbind(coef(m4), coef(m22), coef(m56)))

par.estim$Temp <-c(4, 22, 56)

library(ggplot2)

ggplot(par.estim) + geom_point(aes(x = Temp, y = log(rho))) +
  geom_smooth(aes(x = Temp, y = log(rho)), method = "lm")

ggplot(par.estim) + geom_point(aes(x = Temp, y = log(kappa))) +
  geom_smooth(aes(x = Temp, y = log(kappa)), method = "lm")

m.rho <- lm(log(rho) ~ Temp, data = par.estim)

m.rho
summary(m.rho)

m.kappa <- lm(log(kappa) ~ Temp + Temp^2, data = par.estim)

summary(m.kappa)

rho.coef <- coef(m.rho)
kappa.coef <- coef(m.kappa)

coef.pk <- data.frame(par = c("ap", "Bp", "ak", "Bk"), 
                      valor = c(rho.coef, kappa.coef))

coef.pk

write.csv(coef.pk, file = "coef_pk.csv")

pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par
mod <- nls(ln.S ~ -(exp(ap  + Bp * Temp) * Time.h)^exp(ak + Bk * Temp), 
           data = datos,
           start = pars,
           lower = c(-12, 0, -1, -0.05),
           upper = c(-3, 0.5, 0.1, 0.1),
           algorithm = "port")
summary(mod)

########################################################################

source("weibull.R")

datos.nuevos <- expand.grid(Tiempo = seq(0, 12, len = 50),
                            Temp = seq(4, 56, len = 50))
knitr::kable(head(datos.nuevos))
weib <- function(Tiempo = NA, Temp = NA, pars = NA){
  ap <- pars$ap
  Bp <- pars$Bp
  ak <- pars$ak
  Bk <- pars$Bk
  
  p <- exp(ap + Bp * Temp)
  k <- exp(ak + Bk * Temp)
  
  S <- exp(- (p * Tiempo) ^ k)
  return(S)
}
pars.1 <- as.list(coef(mod))

Sup <- weib(Tiempo = datos.nuevos$Tiempo,
            Temp = datos.nuevos$Temp,
            pars = pars.1)

datos.nuevos$Sup <- Sup

library(lattice)

wireframe(Sup ~ Tiempo + Temp, 
          data = datos.nuevos,
          drape = T,
          screen = list(z = -135, x = -70, y = 3))

library(terra)

bio1 = rast("C:/Users/Ricardo/Documents/Ecología/ENES/6to sem/Ec_Teórica/Gerardo/ETII/Bio1.tif")
plot(bio1)

coef.pk <- read.csv("C:\\Users\\Ricardo\\Documents\\Ecología\\ENES\\6to sem\\Ec_Teórica\\Gerardo\\ETII\\coef_pk.csv") #Importando coeficientes
source("C:\\Users\\Ricardo\\Documents\\Ecología\\ENES\\6to sem\\Ec_Teórica\\Gerardo\\ETII\\weibull.R") #Importando la función

bio1.df <- as.data.frame(bio1, xy = T)
pars <- as.list(coef.pk$valor)
names(pars) <- coef.pk$par

Sup.bio1 <- weib(Tiempo = 2,
                 Temp = bio1.df$Bio1,
                 pars = pars)

bio1.df$Sup <- Sup.bio1
Sup.r <- rast(bio1.df[, -3])
plot(Sup.r)

weib.ode <- function(t, y, params){
  with(params, {
    S <- y
    
    Temp <- Tmax - Tdif * cos(pi * t/24)^2
    
    rho <- exp(ap + Bp * Temp)
    
    kappa <- exp(ak + Bk * Temp)
    
    dS <- - rho * kappa * (-log(S))^(1 - 1/kappa) * S
    
    list(c(dS))
  })
}

setwd("C:\\Users\\Ricardo\\Documents\\Ecología\\ENES\\6to sem\\Ec_Teórica\\Gerardo\\ETII")
getwd()

Tmax <- rast("Tmax-01.tif")
Tmin <- rast("Tmin-01.tif")

Tdif <- Tmax - Tmin

Tmax.df <- as.data.frame(Tmax, xy = T)
Tdif.df <- as.data.frame(Tdif, xy = F)

Temps.df <- data.frame(Tmax.df, Tdif = Tdif.df$`Tmax-01`)
names(Temps.df) <- c("x", "y", "Tmax", "Tdif")

library(deSolve)

y <- 0.9999
t <- seq(0, 12, length(100))

library(tidyr)

library(doParallel)

registerDoParallel(cores = 6)

sims <- foreach(i = 1:nrow(Temps.df), .combine = c) %dopar% {
  library(deSolve)
  params <- pars
  params$Tmax <- Temps.df$Tmax[i]
  params$Tdif <- Temps.df$Tdif[i]
  out <- lsoda(y = y, times = t,
               func = weib.ode,
               parms = params)
  return(out[nrow(out), 2])
}

Sup.fluc <- data.frame(Temps.df[, c("x", "y")], Sup = sims)
Sup.fluc.r <- rast(Sup.fluc)
plot(Sup.fluc.r)

