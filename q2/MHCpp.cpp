#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <time.h>
#include <iostream>

using namespace Rcpp;
using namespace std;
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

//n_sim: number of simulations
// [[Rcpp::export]]
arma::mat MHcpp(const unsigned int n_sim){
  int des_size = 1000;
  arma::mat out(des_size, 2);
  arma::mat sigma(2, 2);
  sigma.zeros(2, 2);
  //set up the covariance matrix
  sigma(0,0) = pow(1.5, 2);
  sigma(1,1) = pow(1.5, 2);
  arma::vec can(2), update(2);
  arma::mat update_dummy;
  //start with {0,0}
  can(0) = 0;
  can(1) = 0;
  double uni, accept, percent;
  
  for (int i = 0; i < n_sim; i++) {
    
    update_dummy = mvrnormArma(1, can, sigma);
    
    update(0) = update_dummy(0);
    update(1) = update_dummy(1);
    //generate a random number from U(0,1) using built-in R function runif()
    uni = ::Rf_runif (0, 1);
    accept = exp(-pow(update(0), 2)-pow(update(0) * update(1), 2)-pow(update(1), 2))/exp(-pow(can(0), 2)-pow(can(0) * can(1), 2)-pow(can(1), 2));
    
    if(i % 1000000 == 0){
      percent = (1.0*i)/(1.0*n_sim)*100.0;
      cout << "i : " << i << " : " << percent << endl;  
    }
    
    if (uni < accept) {
      can(0) = update(0);
      can(1) = update(1);
      if(i > n_sim - des_size){
        out(i, 0) = update(0);
        out(i, 1) = update(1);
      }
      
    } else {
      if(i > n_sim - des_size){
        out(i, 0) = can(0);
        out(i, 1) = can(1);
      }
    }
  }
  return (out);
}

//random multivariate normal sample generator using RcppArmadillo
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}