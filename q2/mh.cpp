#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <time.h>
#include <iostream>

using namespace Rcpp;
using namespace std;
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);

// [[Rcpp::export]]
double dweibull_cpp(double x, double shape,double scale){
  double l = scale;
  double k = shape;
  return (k/l)*pow(x/l,k-1)*exp(-pow(x/l,k));
}

double dexp_cpp(double x, double rate){
  return rate*exp(-rate*x);
}

// [[Rcpp::export]]
double log_likelihood(arma::mat theta, arma::vec y, arma::vec L){
  double b = exp(theta(0,0));
  double sigma = exp(theta(0,1));
  double xi = exp(theta(0,2))/(1+exp(theta(0,2)));
  double S =	0;
  double tmp = 0;
  
  for(unsigned int i=0; i<y.size(); i++){
    tmp = log(dweibull_cpp(y(i),1/b,pow(L(i),-b*xi)*sigma));
     if(!ISNAN(tmp) && tmp != arma::math::inf()){
       S = S + tmp;
     }
  }
  return S;
}

// [[Rcpp::export]]
double log_prior(arma::mat theta, double test_mean){
  double xi = exp(theta(0,2))/(1+exp(theta(0,2)));
  
  double log_pi0_log_b = log(1);
  double log_pi0_log_sigma = log(1);
  double log_pi0_xi = log(dexp_cpp(-log(xi),test_mean));
  return log_pi0_log_b + log_pi0_log_sigma + log_pi0_xi;
}

// [[Rcpp::export]]
double posterior_cpp(arma::mat theta, arma::vec y, arma::vec L, double test_mean){
  return log_prior(theta,test_mean) + log_likelihood(theta,y,L);
}

//random multivariate normal sample generator using RcppArmadillo
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
} 

// [[Rcpp::export]]
arma::mat proposalFunction(arma::vec theta, arma::mat sigma_prop){
  return mvrnormArma(1,theta, sigma_prop);
}

// [[Rcpp::export]]
arma::mat run_MCMC_cpp(arma::vec startValue, 
                       unsigned int iterations,
                       unsigned int burning,
                       arma::vec y,
                       arma::vec L,
                       arma::mat sigma_prop,
                       double test_mean){
  // Define output chain 
  arma::mat out(iterations-burning,3);
  // Define a chain for burning
  arma::mat out2(1,3);
  out2(0,0) = startValue(0);
  out2(0,1) = startValue(1);
  out2(0,2) = startValue(2);
  
  arma::mat update_dummy;
  double probab,uni,percent;
  for(unsigned int i=0; i< iterations; i++){
    update_dummy = mvrnormArma(1, startValue, sigma_prop);
    probab = exp(posterior_cpp(update_dummy,y,L,test_mean) - posterior_cpp(out2,y,L,test_mean));
    
    if(ISNAN(probab)){continue;}
    
    if(i % 1000 == 0){
      percent = (1.0*i)/(1.0*iterations)*100.0;
      cout << "i : " << i << " : " << percent << endl;  
    }
    
    uni = ::Rf_runif (0, 1);
    
    if(uni < probab){
      out2(0,0) = update_dummy(0,0);
      out2(0,1) = update_dummy(0,1);
      out2(0,2) = update_dummy(0,2);
    }
    
    if(i >= burning){
      out(i-burning,0) = out2(0,0);
      out(i-burning,1) = out2(0,1);
      out(i-burning,2) = out2(0,2);
    }
    
    startValue(0) = out2(0,0);
    startValue(1) = out2(0,1);
    startValue(2) = out2(0,2);
  }
  return out;
}

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

