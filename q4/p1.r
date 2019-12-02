set.seed(2)
s_data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)
library(mvtnorm)
library(hydroApps)
#library(CoSMoS)

y <- s_data$strength
L <- s_data$length

#cdf <- expression(1- ((1+(1/k)*(y/(sigma*L^nu))^(1/b))^(-k)) )
#pdf <- deriv(cdf,c("y"),function.arg=c("y","k","sigma","L","nu","b"))

pdf <- function(y,k,sigma,L,nu,b){
  -(1/b)*(1/(sigma*L^nu))*(y/(sigma*L^nu))^(1/b-1)*(1+(1/k)*(y/(sigma*L^nu))^(1/b))^(-k-1)
}

nll <- function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  k <- exp(theta[3])
  nu <- theta[4]
  lambda <- sigma*L^nu
  -1/k -> k
  c <-  1/b
  #print(c(lambda,k,c))
  #print(b)
  #-sum(log(dBurrXII(y,lambda,k,c)))
  #probs <- dBurrXII(y,lambda,k,c)
  #probs <- replace(probs,probs>1,1)
  #-sum(probs)
  res <- pdf(y,k,sigma,L,nu,b)
  #print(res)
  #-sum(log(attr(res,"gradient")))
  -sum(log(res))
}

#lambda = 1.514316e-03
#k = -6.977443e-19
#c = 6.628443e-01

#lambda = 1.574455e-29
#k = -3.919489e-36
#c = 6.392594e-01

#dBurrXII(y,lambda,k,c)
#(1-k*(y[256]/lambda)^c)^(1/k)

nll(c(1,1,1,1),y[2],L[2])

expr <- expression( 
  -log(- exp(th3) *(1+ (1/exp(th3))*(y/(exp(th2)*L^(th4)))^(exp(-th1)))^(-exp(th3)-1) +
         (1/(exp(th1+th2+th3) *L^(th4)) *(y/(exp(th2)*L^(th4)))^(exp(-th1)-1) )   )
)
expr_ <- deriv(expr,c("th1","th2","th3","th4"),function.arg=c("th1","th2","th3","th4","y","L"),hessian = T)


nll_ <- function(theta,y,L){
  # Gradient
  res <- expr_(theta[1],theta[2],theta[3],theta[4],y,L)
  apply(attr(res,"gradient"),2,sum)
}


max_ll <- 2^31-1
counter = 0
th <- c(0,0,0,0)
while(counter < 10){
  # Initialise random starting point
  th0 <- rmvnorm(1,c(-1,1,1,0),diag(c(1,1,2,1)))
  m1 <- try(optim(th0,nll,y=y,L=L,hessian = T,method="BFGS")) # 229.0827
  if(typeof(m1) != "list") next  # Error handling: rmvnorm can pick a very bad starting point where nll is infinite
  # when that happens R quits the loop. 
  # Instead, make R pick another (hopefully) better starting point
  if(m1$value < max_ll){
    counter = 1
    th <- m1$par
    max_ll <- m1$value
    print(c(counter,m1$value,m1$par))
  }else if(round(m1$value,8) == round(max_ll,8)){
    counter = counter + 1
    print(c(counter,max_ll,m1$par))
  }
}





# Plot the data
b <- exp(th[1])
sigma <- exp(th[2])
k <- exp(th[3])
nu <- th[4]
k_ <- -1/k
c <-  1/b

xx <- seq(0,7,length.out = 100)
col_list <- c("red","orange","green","purple")  # List of colours
loop_count <- 1  # Counter for the loop
for(i in unique(s_data$length))
{
  lambda <- sigma*i^nu
  m <- quantile(s_data$strength[s_data$length==i],seq(1,0,length.out = 100))
  y_ <- (1-k_*(xx/lambda))^(1/k_)
  if(i==1){
    plot(xx,y_,type="l",col=col_list[loop_count],xlab="y (Stress in giga-pascals)",ylab=bquote(S[L](y)),main="Visual validation of the model")
    points(as.vector(m),seq(0,1,length.out = 100),col=col_list[loop_count])
  }
  else{
    lines(xx,y_,type="l",col=col_list[loop_count])
    points(as.vector(m),seq(0,1,length.out = 100),col=col_list[loop_count])
  }
  loop_count <- loop_count + 1
}

legend('topright',c('Model','Data'),lty=c(1,NA),pch=c(NA,'o'),bg='white',ncol=1,col=c("black"))
legend('topleft',c('Length = 1mm','Length = 10mm','Length = 20mm','Length = 50mm'),pch=15,bg='white',ncol=1,col=col_list)
