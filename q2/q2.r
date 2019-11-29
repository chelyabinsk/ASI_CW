strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)
library(mvtnorm)

t <- strength.data$strength
x <- strength.data$length

nll <- function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  -sum(dweibull(y,shape = 1/b,scale = L^(-b*xi)*sigma,log = T))
}

expr <- expression(log((b/(sigma*L^(-b*xi)))*(y/(sigma*L^(-b*xi)))^(b-1)*exp(-(y/(sigma*L^(-b*xi)))^(b))))
expr_ <- deriv(expr,c("b","sigma","xi"),function.arg=c("b","sigma","L","y","xi"),hessian = T)


nll_ <- function(theta,y,L){
  # Gradient
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  res <- expr_(b,sigma,y,L,xi)
  apply(attr(res,"gradient"),2,sum)
}



max_ll <- 2^31-1
counter = 0
th <- c(0,0,0)
while(counter < 10){
  # Initialise random starting point
  th0 <- rmvnorm(1,c(-1,1,1),diag(c(1,1,1)))
  m1 <- try(optim(th0,nll,gr = nll_,y=t,L=x,hessian = T))#,method="BFGS")) # 229.0827
  if(typeof(m1) != "list") next
  
  
  print(cat(counter,max_ll,m1$value))
  if(m1$value < max_ll){
    counter = 1
    th <- m1$par
    max_ll <- m1$value
    print(cat(m1$value,th0,""))
  }else if(round(m1$value,2) == round(max_ll,2)){
    counter = counter + 1
  }
}

b <- exp(m1$par[1])
sigma <- exp(m1$par[2])
xi <- exp(m1$par[3])/(1+exp(m1$par[3]))

xx <- seq(0,7,length.out = 100)
val <- 50
m <- mean(strength.data$strength[strength.data$length==val])
m <- quantile(strength.data$strength[strength.data$length==val],c(0.25,0.5,0.75))
L <- val
y_ <- 1-exp(-L^(xi)*(xx/sigma)^(1/b))
plot(xx,y_,type="o")
lines(xx,rep(0.5,100))
lines(xx,rep(0.25,100))
lines(xx,rep(0.75,100))
lines(rep(m[2],100),seq(0,0.5,length.out = 100))
lines(rep(m[1],100),seq(0,0.25,length.out = 100))
lines(rep(m[3],100),seq(0,0.75,length.out = 100))


