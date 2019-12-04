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


# Plot the data
b <- exp(th$par[1])
sigma <- exp(th$par[2])
xi <- exp(th$par[3])/(1+exp(th$par[3]))

xx <- seq(0,7,length.out = 100)
col_list <- c("red","orange","green","purple")  # List of colours
loop_count <- 1  # Counter for the loop
for(i in unique(s_data$length))
{
  m <- quantile(s_data$strength[s_data$length==i],seq(1,0,length.out = 100))
  y_ <- exp(-i^(xi)*(xx/sigma)^(1/b))
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
