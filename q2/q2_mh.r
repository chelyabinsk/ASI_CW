set.seed(42)
# Second part of the second question
strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)
library(mvtnorm)

y <- strength.data$strength
L <- strength.data$length

theta1 <- seq(-3,0,length=100)
theta2 <- seq(-0,5,length=100)

# Plot log-likelihood countours

log.likelihood<-function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  S <- dweibull(y,shape = 1/b,scale = L^(-b*xi)*sigma,log = T)
  #print(sum(S[abs(S)!=Inf]))
  sum(S[abs(S)!=Inf])
}


rl<-matrix(0,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    rl[i,j]<-log.likelihood(c(theta1[i],theta2[j],2.151824),y,L)
  }
}

rang<-quantile(rl,probs = c(0.8,1))
levs.l<-seq(rang[1],rang[2],length.out = 10)
contour(theta1,theta2,rl,levels=levs.l, xlab=expression(theta[1]),
        ylab=expression(theta[2]),col="blue",main="Contours of the log likelihood")


# MH

log.post <- function(theta,y,prior.mean) {
  ## theta = c(log(b),log(sigma),log(xi/(1-xi)));
  b <- exp(theta[1])
  sigma <- exp(theta[2]) #we use log sigma
  xi <- exp(theta[3])/(1-exp(theta[3]))
  # Now the priors
  log_pi0_b<-log(1)
  log_pi0_log.sigma<-log(1)
  #log_pi0_xi<-dexp(theta[3],prior.mean,log=TRUE)
  log_pi0_xi<-dexp(-log(xi),prior.mean,log=TRUE)
  # we add the priors since we assume they are independent
  log_pi0<-log_pi0_b+log_pi0_log.sigma+log_pi0_xi
  # Now the log likelihood
  log.lik<-sum(dweibull(y,shape = 1/b,scale = L^(-b*xi)*sigma,log = T))
  # now the log posterior = log likelihood +log prior
  log.lik+log_pi0
}


theta0 <- c(1,1,1) # initial value

n.rep <- 10^8 # total Number of iterations
burn.in<-1:10^7 # iterations to be discarded
MH<-function(theta0,sigma.prop,n.rep,y,prior.mean){
  # nu.prop is a positive integer which defines the
  # support of the proposal. nu prop is like a standard deviation
  # in the sense that, the higher its value the more spread
  # there is in the proposal
  p<-3
  accept <- rep(0,n.rep)
  theta <- matrix(0,n.rep,p)
  theta[1,] <- theta0
  lp0 <- log.post(theta0,y,prior.mean)
  for (i in 2:n.rep){
    current_theta<-theta[i-1,]
    print(c(current_theta,round(i/n.rep*100,2)))
    proposed_theta <- current_theta+c(rnorm(p,sd=sigma.prop))
    lp1 <- log.post(proposed_theta,y,prior.mean)
    if(is.na(lp1)){i=i-1;next}
    acc <- exp(min(0,lp1-lp0))
    if (runif(1) <= acc) { # accept
      theta[i,]<-proposed_theta
      accept[i] <- 1
      lp0 <- lp1 ## Keep ll0 in sync with th
    } else { ## reject
      theta[i,] <- current_theta
      lp1 <- lp0 ## Keep ll1 in sync with th
    }
  }
  list(theta=theta,ar=mean(accept))
}


mh<-MH(theta0,sigma.prop=c(0.9,0.9,0.9),n.rep,y=nhtemp,prior.mean=0.7)
mh$ar






show.plot<- (n.rep-1000):n.rep
par(mfrow=c(3,1),mar=c(3,4,1,1))
plot(mh$theta[show.plot,1],type="l",ylab=expression(mu))
plot(mh$theta[show.plot,2],type="l",ylab=expression(log(sigma)))
plot(mh$theta[show.plot,3],type="l",ylab=expression(nu))

#pairs(mh$theta,pch=".")

par(mfrow=c(2,2),mar=c(4,4,1,1))
acf(mh$theta[-burn.in,1],xlab=expression(mu))
acf(mh$theta[-burn.in,2],xlab=expression(log(sigma)))
acf(mh$theta[-burn.in,3],xlab=expression(nu))