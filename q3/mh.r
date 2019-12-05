#https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

set.seed(42)
# Second part of the second question
strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)

y <- strength.data$strength
L <- strength.data$length

theta1 <- seq(1.53217965-0.2,1.53217965+0.32,length=100)
theta2 <- seq(-0.01393004-0.3,-0.01393004+0.2,length=100)

# Plot log-likelihood countours
par(mfrow=c(1,1))

log.likelihood<-function(theta,y,L){
  eta <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  tau <- theta[4]
  b <- eta*L^tau
  sum(dweibull(y,shape = b,scale = sigma*L^(-xi/b),log = T),na.rm = T)
}

sol <- c(1.70484974,1.53217965,2.03697199,-0.01393004)
rl<-matrix(0,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    # 1.70484974   1.53217965   2.03697199  -0.01393004
    rl[i,j]<-log.likelihood(c(1.70484974,theta1[i],2.03697199,theta2[j]),y,L)
  }
}

rang<-quantile(rl,probs = c(0.3,0.7))
levs.l<-seq(rang[1],rang[2],length.out = 10)
contour(theta1,theta2,rl,levels=levs.l, xlab=expression(theta[1]),
        ylab=expression(theta[2]),col="blue",main="Contours of the log likelihood")

log.prior <- function(theta,y,L,test_mean){
  eta <- exp(theta[1])
  sigma <- exp(theta[2])
  xi <- exp(theta[3])/(1+exp(theta[3]))
  tau <- theta[4]
  b <- eta*L^tau
  
  log_pi0_log.eta<-log(1)
  log_pi0_log.sigma<-log(1)
  log_pi0_log.xi <- dexp(-log(xi),test_mean[1],log=T)
  log_pio_log.tau <- -abs(tau)/test_mean[2]
  
  #print(log_pi0_log.eta + log_pi0_log.sigma + log_pi0_log.xi + log_pio_log.tau)
  
  log_pi0_log.eta + log_pi0_log.sigma + log_pi0_log.xi + log_pio_log.tau
}

posterior <- function(theta,y,L,test_mean){
  #print(log.likelihood(theta,y,L))
  log.likelihood(theta,y,L) + log.prior(theta,y,L,test_mean)
}

#posterior(c(-0.7356172,2.852293,1.730333,-1.795238),y,L,test_mean)

##MH

proposonalFunction <- function(theta,sigma.prop){
  mvtnorm::rmvnorm(1,mean=theta,sigma=sigma.prop)
}


run_metropolis_MCMC <- function(startvalue, iterations,y,L,sigma.prop,test_mean){
  chain = array(dim = c(iterations+1,4))
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    proposal = proposonalFunction(chain[i,],sigma.prop)
    
    #print(exp(posterior(proposal,y,L,test_mean)))
    #break
    probab = exp(posterior(proposal,y,L,test_mean) - posterior(chain[i,],y,L,test_mean))
    #print(proposal)
    if(is.na(probab)){next}
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}

#########################################################################
# 1.70484974   1.53217965   2.03697199  -0.01393004
startvalue = c(1.70484974,1.53217965,2.03697199,-0.01393004)
sigma.prop = matrix(c(1,0,0,0,
                      0,1,0,0,
                      0,0,1,0,
                      0,0,0,1), ncol=4)
test_mean = c(3,200)
N <- 50000
burnIn = N*0.75
#########################################################################


chain = run_metropolis_MCMC(startvalue, N,y,L,sigma.prop,test_mean)


acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))



plot.mcmc<-function(mcmc.out)
{
  op=par(mfrow=c(2,2))
  plot(ts(mcmc.out),col=2)
  hist(mcmc.out,30,col=3)
  qqnorm(mcmc.out,col=4)
  abline(0,1,col=2)
  acf(mcmc.out,col=2,lag.max=100)
  par(op)
}

plot.mcmc(chain[-(1:burnIn),1])
plot.mcmc(chain[-(1:burnIn),2])
plot.mcmc(chain[-(1:burnIn),3])
plot.mcmc(chain[-(1:burnIn),4])

acceptance
