# https://ysquared2.wordpress.com/2014/03/27/metropolis-hasting-algorithm-an-example-with-cc-and-r/
# http://dirk.eddelbuettel.com/code/rcpp/RcppGSL.pdf


library(Rcpp)
library(gsl)
library(RcppGSL)
library(RcppArmadillo)

library(Rcpp)
sourceCpp("C:\\Users\\Seva\\Documents\\R\\ASI\\ASI_CW-master\\q2\\mh.cpp")

MHcpp(1e2)

