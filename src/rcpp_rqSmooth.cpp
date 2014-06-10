//#include <Rcpp.h>
#include "rcpp_smoothness_penalty.h"
using namespace std;

//' Extends the parameter vector to the natural dimension
//' 
//' @param beta the parameter vector
//' @param blockLength the vector of blockLengths corresponding to each entry in beta
//' @return the extended parameter vector
//[[Rcpp::export]]
NumericVector completeBetaCpp (NumericVector beta, NumericVector blockLength) {
  if (beta.length() != blockLength.length()) stop("beta and blockLength must have same size");
  int n = beta.length();
  NumericVector sumIndex = cumsum(blockLength);
  vector<int> indices (1, 0);// = Rcpp::as<std::vector<int>>(tmp);
  NumericVector tmp(sumIndex(n-1));
  for(int i = 0; i < n; i++) {
    indices.push_back(sumIndex[i]);
  }
  for(int i = 0; i < n; i++) {
    for(int j = indices[i]; j < indices[i+1]; j++) {
      tmp[j] = beta[i];
    }
    //tmp[Range(indices[i], indices[i+1])] = beta[i];
  }
  NumericVector newBeta = Rcpp::wrap(tmp);
  return newBeta;
}

//' Computes the smoothed objective criterion after Muggeo.etal2012
//' 
//' @param y the data
//' @param beta the parameter vector
//' @param blockLength the blocklength vector
//' @param tau the quantile 
//[[Rcpp::export]]
double rqSmooth(NumericVector y, NumericVector beta, NumericVector blockLength, float tau, float lambda = 0.0,
    std::string method = "1order", float c = 0.1) {

    NumericVector cBeta  = completeBetaCpp(beta, blockLength);
    NumericVector r      = y - cBeta;
    float         convex = c * tau * (1-tau) / 2;
    double        sumRes = 0.0;

    std::for_each(r.begin(), r.end(), [&tau, &c, &convex, &sumRes](double x) {
            if (x <= -c * tau) sumRes += (tau - 1) * x;
            else if (x >= c * ( 1-tau )) sumRes += tau * x;
            else if (x > -c * tau && x <= 0) sumRes += (1 - tau) * std::pow(x, 2) / (2 * c * tau) + convex ;
            else sumRes += tau * std::pow(x, 2) / (2 * c *(1 - tau)) + convex ;
            });

    return sumRes + penaltyCpp(beta, lambda, method);
}

//' Computes the smoothed objective criterion for multiple samples
//' 
//' @param y the data
//' @param beta the parameter vector
//' @param blockLength the blocklength vector
//' @param tau the quantile 
//[[Rcpp::export]]
double rqSmoothRegional(NumericMatrix y, NumericVector beta, NumericVector blockLength, float tau, float lambda = 0.0,
    std::string method = "1order", float c = 0.1) {
  
  int    ncols = y.ncol();
  double sumRes = 0.0;
  for (int i = 0; i < ncols; i++) {
    sumRes += rqSmooth(y(_,i), beta, blockLength, tau, 0, method, c);
  }
  return sumRes + penaltyCpp(beta, lambda, method);
}

NumericVector d1_rqSmooth(NumericVector y, NumericVector beta, NumericVector blockLength, float tau,  float lambda,
  std::string method, float c) {
  
  int n = beta.length();
  NumericVector sumIndex = cumsum(blockLength);
  vector<int> indices (1, 0);// = Rcpp::as<std::vector<int>>(tmp);
  //NumericVector tmp(sumIndex(n-1));
  for(int i = 0; i < n; i++) {
    indices.push_back(sumIndex[i]);
  }

  NumericVector cBeta = completeBetaCpp(beta, blockLength);
  NumericVector r = clone(y) - cBeta;
  NumericVector x(n);
  for (int i = 0; i < n; i++) {
    for (int j = indices[i]; j < indices[i+1]; j++) {
      if (r[j] <= -c * tau) x[i] += 1 - tau;
      else if (r[j] >= c * (1-tau)) x[i] += - tau;
      else if (r[j] > -c * tau && r[j] <= 0) x[i] += (1 - tau) * (beta[i] - y[j]) / (c * tau);
      else if (r[i] >  0 && r[i] < c * (1 - tau)) x[i] += tau * (beta[i] - y[j]) / (c * (1 - tau));
    }
  }
  
  return x + d1Penalty(beta, lambda, method);
}

NumericVector d1_rqSmoothRegional(NumericMatrix y, NumericVector beta,
   NumericVector blockLength, float tau,  float lambda, std::string method, float c) {
  
  int           ncols = y.ncol();
  int           nrows = beta.length();
  NumericMatrix temp(nrows, ncols);
  NumericVector derivative(nrows);

  for (int i = 0; i < ncols; i++) {
    temp(_,i) = d1_rqSmooth(y(_,i), beta, blockLength, tau, 0, method = method, c = c);
  }
  
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      derivative[i] += temp(i, j);
    }
  }//*/
  
  return derivative + d1Penalty(beta, lambda, method);
}