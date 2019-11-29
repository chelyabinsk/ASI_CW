# https://ysquared2.wordpress.com/2014/03/27/metropolis-hasting-algorithm-an-example-with-cc-and-r/
# http://dirk.eddelbuettel.com/code/rcpp/RcppGSL.pdf


library(Rcpp)
library(gsl)
library(RcppGSL)
library(RcppArmadillo)

library(Rcpp)
sourceCpp("H:/Applied Statistical Inference/Coursework/q2/MHCpp.cpp")

library(MASS)
target <- function(x){
  exp(-x[1]^2-x[1]^2*x[2]^2-x[2]^2)
}

MHR<- function(n.sim, sd){
  mat <- matrix(0,ncol=2,nrow=n.sim)
  obs <- c(0,0)
  sig <- matrix(0, 2,2)
  diag(sig) <- c(1,1)*sd
  mat[1,] <- obs
  prior <- function(x,mu){
    (2*pi)^(-1)/sd*exp(-0.5/sd^2*sum((x-mu)^2))
  }
  for (i in 2:n.sim){
    can <- mvrnorm(n = 1, mu=obs, Sigma=sig)
    accept <- min(1,target(x=can)/target(x=obs))
    foo <- runif(1)
    if (foo <=accept){
      obs <- can
      mat[i,] <- can
    } else{
      mat[i,] <- obs
    }
  }
  mat <- data.frame(mat)
  names(mat) <- c("x","y")
  mat
}

MHcpp(1e8)
#system.time(sim.mhr <- MHR(n.sim=1e6,sd=1.5))
