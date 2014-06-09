#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include "convexHull.h"


using namespace std;
using namespace Rcpp;

//' The length of a string (in characters).
//'
//' @param str input character vector
//' @return characters in each element of the vector
// [[Rcpp::export]]
NumericVector rcpp_convex_minorant(NumericVector x, NumericVector y) {
//NumericVector convexMinorant(NumericVector y, NumericVector x) {
  
  //Rcout << "Convex minorant started" << std::endl;
  int ny = y.length();
  //Rcout << "ny equals" << ny << std::endl;

  NumericVector XX = x; //cumsum(x);
  NumericVector XY = y; //cumsum(y);
  
  vector<Point> P(ny + 1);
  P[0].x = 0.0, P[0].y =0.0;
  
  for (int i = 1; i < ny + 1; i++) {
    P[i].x = XX[i-1];
    P[i].y = XY[i-1];
  }
  
  //Rcout << "Points are initialized" << std::endl;
  vector<Point> convHull = convex_hull(P);
  //Rcout << "Convex minorant computed" << std::endl;
  
  int            nP = convHull.size();
  //Rcout << "nP equals " << nP << std::endl;
  vector<int>    convHullX(nP + 1); 
  vector<double> convHullY(nP + 1); 
  convHullX[0] = 0;                 
  convHullY[0] = 0;                 
  for (int i = 0; i < nP; i++) {
    convHullX[i+1] = convHull.at(i).x;
    convHullY[i+1] = convHull.at(i).y;
  }
  //Rcout << "Points separated" << std::endl;
  
  NumericVector XXX  = Rcpp::wrap(convHullX);
  NumericVector XYY  = Rcpp::wrap(convHullY);
  NumericVector xDiff = diff(XXX);
  NumericVector yDiff = diff(XYY);
  
  //Rcout << "size xDiff =" << xDiff.length() << std::endl;
  //Rcout << "np =" << nP << std::endl;
  
  //Rcout << "Some computations" << std::endl;
  vector<double> tmp;
  for (int i = 0; i < nP; i++) {
    for (int j = 0; j < xDiff[i]; j++) {
      //Rcout << "yDiff is" << yDiff[i] << std::endl;
      //Rcout << "xDiff is" << xDiff[i] << std::endl;
      tmp.push_back(yDiff[i] / xDiff[i]);
    }
  }
  //if(tmp[ny-1] < tmp[ny-2]) tmp[ny-1] = tmp[ny-2];
  
  //Rcout << "Convex minorant algoritm finished" << std::endl;
  
  NumericVector newBeta = Rcpp::wrap(tmp);
  return newBeta;
  //return List::create(Named("slopes") = newBeta, Named("x") = XXX, Named("y") = XYY);
}