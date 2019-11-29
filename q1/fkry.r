set.seed(1230)  # Reproducibility

library(mvtnorm)

s_data <-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)

y <- s_data$strength
L <- s_data$length

nll <- function(theta,y,L){
  # Negative log-likelihood
  b <- exp(theta[1])  # Recover b>0
  sigma <- exp(theta[2])  # Recover sigma>0
  -sum(dweibull(y,shape = 1/b,scale = L^(-b)*sigma,log = T))
}

# https://stats.stackexchange.com/questions/22787/em-maximum-likelihood-estimation-for-weibull-distribution
expr <- expression(-log(-L/(exp(th2)*exp(th1)) ) + (1-1/(exp(th1) ))*log(y/exp(th2)) +  L*(y/exp(th2))^(1/exp(th1) ))
expr_ <- deriv(expr,c("th1","th2"),function.arg=c("th1","th2","y","L"),hessian = T)

nll_ <- function(theta,y,L){
  # Gradient
  res <- expr_(theta[1],theta[2],y,L)
  apply(attr(res,"gradient"),2,sum)
}



max_ll <- 2^31-1
counter = 0
th <- c(0,0)
while(counter < 10){
  # Initialise random starting point
  th0 <- rmvnorm(1,c(-1,1),diag(c(1,1)))
  m1 <- try(optim(th0,nll,gr=nll_,y=y,L=L,hessian = T,method="BFGS")) # 229.0827
  if(typeof(m1) != "list") next  # Error handling: rmvnorm can pick a very bad starting point where nll is infinite
                                 # when that happens R quits the loop. 
                                 # Instead, make R pick another (hopefully) better starting point
  if(m1$value < max_ll){
    counter = 1
    th <- m1$par
    max_ll <- m1$value
    print(c(counter,m1$value,m1$par))
  }else if(m1$value == max_ll){
    counter = counter + 1
    print(c(counter,max_ll,m1$par))
  }
}



# Plot the data
b <- exp(m1$par[1])
sigma <- exp(m1$par[2])

xx <- seq(0,7,length.out = 100)
col_list <- c("red","orange","green","purple")  # List of colours
loop_count <- 1  # Counter for the loop
for(i in unique(s_data$length))
{
  m <- quantile(s_data$strength[s_data$length==i],seq(1,0,length.out = 100))
  y <- exp(-i*(xx/sigma)^(1/b))
  if(i==1){
    plot(xx,y,type="l",col=col_list[loop_count],xlab="y (Stress in giga-pascals)",ylab=bquote(S[L](y)),main="Visual validation of the model")
    points(as.vector(m),seq(0,1,length.out = 100),col=col_list[loop_count])
  }
  else{
    lines(xx,y,type="l",col=col_list[loop_count])
    points(as.vector(m),seq(0,1,length.out = 100),col=col_list[loop_count])
  }
  loop_count <- loop_count + 1
}

legend('topright',c('Model','Data'),lty=c(1,NA),pch=c(NA,'o'),bg='white',ncol=1,col=c("black"))
legend('topleft',c('Length = 1mm','Length = 10mm','Length = 20mm','Length = 50mm'),pch=15,bg='white',ncol=1,col=col_list)


#lines(xx,rep(0.5,100))
#lines(xx,rep(0.25,100))
#lines(xx,rep(0.75,100))
#lines(rep(m[2],100),seq(0,0.5,length.out = 100))
#lines(rep(m[1],100),seq(0,0.75,length.out = 100))
#lines(rep(m[3],100),seq(0,0.25,length.out = 100))


#library(survival)

#fit.surv<-survreg(Surv(strength) ~length,data=s_data, dist="weibull")
#summary(fit.surv)
#b_ <- 1/fit.surv$scale
#a_ <- exp(coef(fit.surv))
