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