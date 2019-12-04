#https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

set.seed(42)
# Second part of the second question
strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)

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
  
  sum(S[abs(S)!=Inf])
}


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

log.prior <- function(theta,y,L,test_mean){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  
  log_pi0_log.b<-log(1)
  log_pi0_log.sigma<-log(1)
  log_pi0_xi<-dexp(-log(xi),test_mean,log=TRUE) # dexp(-log(xi),2,log = T)
  log_pi0_log.b + log_pi0_log.sigma + log_pi0_xi
}

posterior <- function(theta,y,L,test_mean){
  log.likelihood(theta,y,L) + log.prior(theta,y,L,test_mean)
}

##MH

proposonalFunction <- function(theta,sigma.prop){
  rnorm(3,mean=theta,sigma.prop)
}

run_metropolis_MCMC <- function(startvalue, iterations,y,L,sigma.prop,test_mean){
  chain = array(dim = c(iterations+1,3))
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    proposal = proposonalFunction(chain[i,],sigma.prop)
    #print(exp(posterior(proposal,y,L,test_mean)))
    #break
    probab = exp(posterior(proposal,y,L,test_mean) - posterior(chain[i,],y,L,test_mean))
    if(is.na(probab)) next
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

###############################################
startvalue = c(-1,0,1)
sigma.prop = c(0.01,0.08,0.25)
test_mean = 3  # mu0>0 (fixed) second bullet point
N <- 90000
burnIn = N*0.9
################################################


chain = run_metropolis_MCMC(startvalue, N,y,L,sigma.prop,test_mean)


acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

## He said something about plotting log(xi) instead of just xi
## when checking for convergence

plot.mcmc<-function(mcmc.out)
{
  op=par(mfrow=c(2,2))
  plot((ts(mcmc.out)),col=2)
  hist(mcmc.out,30,col=3)
  qqnorm(mcmc.out,col=4)
  abline(0,1,col=2)
  acf(mcmc.out,col=2,lag.max=100)
  par(op)
}

plot.mcmc(chain[-(1:burnIn),1])
plot.mcmc(chain[-(1:burnIn),2])
plot.mcmc(chain[-(1:burnIn),3])

acceptance
