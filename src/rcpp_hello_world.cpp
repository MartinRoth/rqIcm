
#include <Rcpp.h>
using namespace Rcpp;

//' I can do my own documentation
// [[Rcpp::export]]
List rcpp_hello_world() {
   
    CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
    NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z            = List::create( x, y ) ;
    
    return z ;
}
