s <- 100
e <- 110
sigma <- 0.75
t <- 1
r <- 0.02

d1 <- function(){return((log(s/e)+(r+sigma**2/2)*t)/(sigma*sqrt(t)))}
d2 <- function(){return(d1()-sigma*sqrt(t))}
v <- function(s,e,r,sigma,t){s*pnorm(d1())-e*exp(-r*t)*pnorm(d2())}
print(v(s,e,r,sigma,t))
