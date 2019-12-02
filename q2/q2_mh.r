set.seed(42)
# Second part of the second question
strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)
library(mvtnorm)

y <- strength.data$strength
L <- strength.data$length

theta1 <- seq(-3,1,length=100)
theta2 <- seq(0,3,length=100)

# Plot log-likelihood countours
par(mfrow=c(1,1))

log.likelihood<-function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  S <- dweibull(y,shape = 1/b,scale = L^(-b*xi)*sigma,log = T)
  #print(sum(S[abs(S)!=Inf]))
  sum(S[abs(S)!=Inf])
}

# 2.184956 -3.927016 -9.902739

rl<-matrix(0,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    rl[i,j]<-log.likelihood(c(theta1[i],theta2[j],2.55397),y,L)
  }
}

rang<-quantile(rl,probs = c(0.9,1))
levs.l<-seq(rang[1],rang[2],length.out = 10)
contour(theta1,theta2,rl,levels=levs.l, xlab=expression(theta[1]),
        ylab=expression(theta[2]),col="blue",main="Contours of the log likelihood")


# MH

log.prior <- function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  
  log_pi0_log.b<-log(1)
  log_pi0_log.sigma<-log(1)
  log_pi0_xi<-dexp(-log(xi),1,log=TRUE) # dexp(-log(xi),2,log = T)
  log_pi0_log.b + log_pi0_log.sigma + log_pi0_xi
}


theta1 <- seq(0,10,length=100)  # x
theta2 <- seq(-0.5,2,length=100)   # y

rl<-matrix(0,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    rl[i,j]<-log.prior(c(theta1[i],2.55397,theta2[j]),y,L)
  }
}

rang<-quantile(rl,probs = c(0.9,1))
levs.l<-seq(rang[1],rang[2],length.out = 10)
contour(theta1,theta2,rl,levels=levs.l, xlab=expression(theta[1]),
        ylab=expression(theta[2]),col="blue",main="Contours of the log likelihood")



log.post <- function(theta,y,L,mean_guess) {
  # Now the priors
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  
  log_pi0_log.b<-log(1)
  log_pi0_log.sigma<-log(1)
  log_pi0_xi<-dexp(-log(xi),mean_guess,log=TRUE) # dexp(-log(xi),2,log = T)
  
  # we add the priors since we assume they are independent
  log_pi0<-log_pi0_log.b + log_pi0_log.sigma + log_pi0_xi
  # Now the log likelihood
  
  log.lik <- -sum(dweibull(y,shape = 1/b,scale = L^(-b*xi)*sigma,log = T))

  # now the log posterior = log likelihood +log prior
  log.lik+log_pi0
}

MH<-function(theta0,sigma.prop,n.rep,y,L,mean_guess){
  # nu.prop is a positive integer which defines the
  # support of the proposal. nu prop is like a standard deviation
  # in the sense that, the higher its value the more spread
  # there is in the proposal
  p<-3
  accept <- rep(0,n.rep)
  theta <- matrix(0,n.rep,p)
  theta[1,] <- theta0
  lp0 <- log.post(theta0,y,L,mean_guess)
  for (i in 2:n.rep){
    current_theta<-theta[i-1,]
    proposed_theta <- current_theta+rnorm(p,sd=sigma.prop)
    if (proposed_theta[3] <=0 ) {
      # definitely reject the whole vector if df becomes negative
      theta[i,] <- current_theta
    } else { # otherwise proceed with Metropolis-Hastings step
      lp1 <- log.post(proposed_theta,y,L,mean_guess)
      acc <- exp(min(0,lp1-lp0))
      if(is.na(acc)) next
      if (runif(1) <= acc) { # accept
        theta[i,]<-proposed_theta
        accept[i] <- 1
        lp0 <- lp1 ## Keep ll0 in sync with th
      } else { ## reject
        theta[i,] <- current_theta
        lp1 <- lp0 ## Keep ll1 in sync with th
      }
    }
  }
  list(theta=theta,ar=mean(accept))
}

theta0 <- c(0,0,0) # initial value
sigma.prop <- 10 # initial proposal std.dev.
n.rep <- 100000 # total Number of iterations
mean_guess <- 2
burn.in<-1:10000 # iterations to be discarded

#MH<-function(theta0,sigma.prop,n.rep,y,L,mean_guess)
o <- MH(theta0,sigma.prop,n.rep,y,L,mean_guess)
o
