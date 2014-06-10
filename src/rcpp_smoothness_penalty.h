#ifndef PENALTY_H
#define PENALTY_H

#include <Rcpp.h>
using namespace Rcpp;

double penalty2order(NumericVector beta);
double penalty1order(NumericVector beta);
double penaltyCpp(NumericVector beta, float lambda, std::string method = "1order");
NumericVector d1Penalty1order(NumericVector beta);
NumericVector d1EilersPenalty2order(NumericVector beta);
NumericVector d1Penalty(NumericVector beta, float lambda, std::string method);

#endif


