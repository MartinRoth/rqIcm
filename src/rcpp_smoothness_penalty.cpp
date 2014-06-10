#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double penalty2order(NumericVector beta) {
  return sum(pow(diff(diff(beta)),2));
}

double penalty1order(NumericVector beta) {
  int n = beta.length()-1;
  return beta[n] - beta[0];
}

//' Computes the penalty term for a given parameter vector
//' 
//' @param beta parameter vector
//' @param lambda multiplicator
//' @param method (eiter "1order" = TV or "2order" = sum Second Derivative)
//' @return penalty
//[[Rcpp::export]]
double penaltyCpp(NumericVector beta, float lambda, std::string method = "1order") {
  if (lambda == 0.0) return 0.0;
  //else return lambda * penalty1order(beta);
  else if (method.compare("1order")   == 0) return lambda * penalty1order(beta);
  else if (method.compare("2order")   == 0) return lambda * penalty2order(beta);
  else Rcpp::stop("Penalty method not defined");
} 

NumericVector d1Penalty1order(NumericVector beta) {
  int n = beta.length();
  vector<double> p;
  p.push_back(-1);
  for (int i = 0; i < n-2; i++) p.push_back(0);
  p.push_back(1);
  NumericVector penalty = Rcpp::wrap(p);
  return penalty;
}

NumericVector d1EilersPenalty2order(NumericVector beta) {
  int n = beta.length();
  vector<double> a;
  vector<double> b;
  vector<double> p(n);
  for (int i = 0; i < n-2; i++) {
    a.push_back( 2 * beta[i+2] - 4 * beta[i+1] + 2 * beta[i]);
    b.push_back(-4 * beta[i+2] + 8 * beta[i+1] - 4 * beta[i]);
  }
  p[0] = a[0];
  p[1] = b[0] + a[1];
  for (int i = 2; i < n-2; i++) {
    p[i] = a[i] + b[i-1] + a[i-2];
  }
  p[n-2] = a[n-4] + b[n-3];
  p[n-1] = a[n-3];
  NumericVector penalty = Rcpp::wrap(p);
  return penalty;
}

NumericVector d1Penalty(NumericVector beta, float lambda, std::string method) {
  if (lambda == 0.0) {
    NumericVector x(beta.length(), 0.0);
    return x;
  }
  //else return lambda * d1Penalty1order(beta);
  else if (method.compare("1order")   == 0) return lambda * d1Penalty1order(beta);
  else if (method.compare("2order")   == 0) return lambda * d1EilersPenalty2order(beta);
}
