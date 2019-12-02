# https://www.rdocumentation.org/packages/hydroApps/versions/0.1-1/topics/pBurrXII

set.seet(23121)
strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)
library(mvtnorm)

y <- strength.data$strength
L <- strength.data$length

nll <- function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  -sum(dweibull(y,shape = 1/b,scale = L^(-b*xi)*sigma,log = T))
}

expr <- expression( 
  
  -log( (-L^(exp(th3)/(1+exp(th3)))) /(exp(th1)*exp(th2)) ) -
    (y/(exp(th2)))^((1/(exp(th1)))-1) +
    (L^(exp(th3)/(1+exp(th3))))*(y/(exp(th2)))^(1/(exp(th1)))
  
)
expr_ <- deriv(expr,c("th1","th2","th3"),function.arg=c("th1","th2","th3","y","L"),hessian = T)


nll_ <- function(theta,y,L){
  # Gradient
  res <- expr_(theta[1],theta[2],theta[3],y,L)
  apply(attr(res,"gradient"),2,sum)
}



max_ll <- 2^31-1
counter = 0
th <- c(0,0,0)
while(counter < 10){
  # Initialise random starting point
  th0 <- rmvnorm(1,c(-1,1,1),diag(c(1,1,2)))
  m1 <- try(optim(th0,nll,gr=nll_,y=y,L=L,hessian = T,method="BFGS")) # 229.0827
  if(typeof(m1) != "list") next  # Error handling: rmvnorm can pick a very bad starting point where nll is infinite
  # when that happens R quits the loop. 
  # Instead, make R pick another (hopefully) better starting point
  if(m1$value < max_ll){
    counter = 1
    th <- m1
    max_ll <- m1$value
    print(c(counter,m1$value,m1$par))
  }else if(round(m1$value,2) == round(max_ll,2)){
    counter = counter + 1
    print(c(counter,max_ll,m1$par))
  }
}

b <- exp(th$par[1])
sigma <- exp(th$par[2])
xi <- exp(th$par[3])/(1+exp(th$par[3]))

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


