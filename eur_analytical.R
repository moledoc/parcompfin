d1 <- function(){return((log(s/e)+(r+sigma**2/2)*t)/(sigma*sqrt(t)))}
d2 <- function(){return(d1()-sigma*sqrt(t))}
v <- function(s,e,r,sigma,t){s*pnorm(d1())-e*exp(-r*t)*pnorm(d2())}
