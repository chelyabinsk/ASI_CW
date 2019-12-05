#https://theoreticalecology.wordpress.com/2010/09/17/metropolis-hastings-mcmc-in-r/

set.seed(42)
# Second part of the 4th question
strength.data<-
  read.table(url("http://people.bath.ac.uk/kai21/ASI/CW2019/strength.txt"),header = TRUE)

y <- strength.data$strength
L <- strength.data$length

sol <- c(-1.5116566,1.5887773,-2.2711463,-0.1815788)

theta1 <- seq(sol[2]-0.2,sol[2]+0.1,length=100)
theta2 <- seq(sol[4]-0.1,sol[4]+0.1,length=100)

# Plot log-likelihood countours
par(mfrow=c(1,1))

pdf <- function(y,k,sigma,L,nu,b){
  (1/b)*(1/(sigma*L^nu))*(y/(sigma*L^nu))^(1/b-1)*(1+(1/k)*(y/(sigma*L^nu))^(1/b))^(-k-1)
}

log.likelihood<-function(theta,y,L){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  k <- exp(theta[3])
  nu <- theta[4]
  lambda <- sigma*L^nu
  -1/k -> k
  res <- pdf(y,k,sigma,L,nu,b)
  sum(log(res),na.rm = T)  # Not sure about this na.rm thing ... :|
}
rl<-matrix(0,nrow=100,ncol=100)
for (i in 1:100){
  for (j in 1:100){
    rl[i,j]<-log.likelihood(c(sol[1],theta1[i],sol[3],theta2[j]),y,L)
  }
}

rang<-quantile(rl,probs = c(0.5,1),na.rm = T)
levs.l<-seq(rang[1],rang[2],length.out = 10)
contour(theta1,theta2,rl,levels=levs.l, xlab=expression(theta[1]),
        ylab=expression(theta[2]),col="blue",main="Contours of the log likelihood")

log.prior <- function(theta,y,L,test_mean){
  b <- exp(theta[1])
  sigma <- exp(theta[2])
  k <- exp(theta[3])
  nu <- theta[4]
  
  log_pi0_log.b<-log(1)
  log_pi0_log.sigma<-log(1)
  log_pi0_log.k <- dexp(1/k,test_mean[1],log=T)
  log_pio_log.nu <- dexp(nu+b,test_mean[2],log=T)
  
  log_pi0_log.b + log_pi0_log.sigma + log_pi0_log.k + log_pio_log.nu
}

posterior <- function(theta,y,L,test_mean){
  log.likelihood(theta,y,L) + log.prior(theta,y,L,test_mean)
}


##MH

proposonalFunction <- function(theta,sigma.prop){
  mvtnorm::rmvnorm(1,mean=theta,sigma=sigma.prop)
}


run_metropolis_MCMC <- function(startvalue, iterations,y,L,sigma.prop,test_mean){
  chain = array(dim = c(iterations+1,4))
  chain[1,] = startvalue
  for (i in 1:iterations){
    
    proposal = proposonalFunction(chain[i,],sigma.prop)
    
    probab = exp(posterior(proposal,y,L,test_mean) - posterior(chain[i,],y,L,test_mean))

    if(is.na(probab)) {next}
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(chain)
}


#################################################
# -1.5116566   1.5887773  -2.2711463  -0.1815788
startvalue = c(-1,1,-2,0)
sigma.prop = matrix(c(1,0,0,0,
                      0,1,0,0,
                      0,0,1,0,
                      0,0,0,1), ncol=4)
N <- 80000
burnIn = N*0.75

test_mean = c(3,3)
#################################################


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
plot.mcmc(chain[-(1:burnIn),4]+chain[-(1:burnIn),1])  # Dunno

acceptance


