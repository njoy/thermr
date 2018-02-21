#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "legndr.h"


TEST_CASE( "generate legendre polynomials x by recursion" ){
  std::vector<double> p {2.5, 8.8, 2.3 };
  double x = 3.4;
  legndr( x, p );
  REQUIRE( 1.00  == Approx( p[0] ).epsilon(1e-6) );
  REQUIRE( 3.40  == Approx( p[1] ).epsilon(1e-6) );
  REQUIRE( 16.84 == Approx( p[2] ).epsilon(1e-6) );
  REQUIRE( true );  

}
 
