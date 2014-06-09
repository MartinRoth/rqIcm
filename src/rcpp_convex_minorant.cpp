#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include "convexHull.h"


using namespace std;
using namespace Rcpp;

//Computes the slopes of the convex minorant given the points of the lower convex hull
//The arguments could be changed to double vectors
NumericVector compute_slopes(NumericVector x, NumericVector y) {
  double        startValue = y[0];
  NumericVector xDiff      = diff(x);
  NumericVector yDiff      = diff(y);
  int           nxd        = xDiff.length();

  vector<double> tmp;
  tmp.push_back(startValue);
  for (int i = 0; i < nxd; i++) {
    for (int j = 0; j < xDiff[i]; j++) {
      tmp.push_back( yDiff[i] / xDiff[i]);
    }
  }
  vector<double> tmp2;
  partial_sum(tmp.begin(), tmp.end(), back_inserter(tmp2));
  return Rcpp::wrap(tmp2);
}

//' Returns the convex minorant of a polygon.
//'
//' @param str input character vector
//' @return characters in each element of the vector
// [[Rcpp::export]]
List rcpp_convex_minorant(NumericVector x, NumericVector y) {
//NumericVector convexMinorant(NumericVector y, NumericVector x) {
  
  //Rcout << "Convex minorant started" << std::endl;
  int ny = y.length();
  //Rcout << "ny equals" << ny << std::endl;

  NumericVector XX = x; //cumsum(x);
  NumericVector XY = y; //cumsum(y);
  
  vector<Point> P(ny);
  //P[0].x = 0.0, P[0].y =0.0;
  
  for (int i = 0; i < ny; i++) {
    P[i].x = XX[i];
    P[i].y = XY[i];
  }
  
  //Rcout << "Points are initialized" << std::endl;
  vector<Point> convHull = convex_hull(P);
  //Rcout << "Convex minorant computed" << std::endl;
  
  int            nP = convHull.size();
  //Rcout << "nP equals " << nP << std::endl;
  vector<int>    convHullX(nP); 
  vector<double> convHullY(nP); 
  //convHullX[0] = 0;                 
  //convHullY[0] = 0;                 
  for (int i = 0; i < nP; i++) {
    convHullX[i] = convHull.at(i).x;
    convHullY[i] = convHull.at(i).y;
  }
  //Rcout << "Points separated" << std::endl;
  
  NumericVector XXX  = Rcpp::wrap(convHullX);
  NumericVector XYY  = Rcpp::wrap(convHullY);
  NumericVector newBeta = compute_slopes(XXX, XYY);
  //return newBeta;
  return List::create(Named("slopes") = newBeta, Named("x") = XXX, Named("y") = XYY);
}