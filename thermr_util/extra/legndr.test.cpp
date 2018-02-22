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


  p = { 0.0, 0.1, 1.2, 1.3, 2.4, 2.5, 3.6, 3.7, 4.8, 4.9 };
  x = 3.4;
  legndr( x, p );

  std::vector<double> correctOut { 1.0, 3.4, 16.84, 93.16000812, 541.672063, 
    3240.505113, 19747.75571, 121915.39950, 759931.407419, 4772079.264976 };
  for ( size_t i = 0; i < correctOut.size(); ++i ){ 
    REQUIRE( correctOut[i]  == Approx( p[i] ).epsilon(1e-6) );
  }
  
  x = 8.5;
  legndr( x, p );
  correctOut = { 1.0, 8.5, 107.875, 1522.5625, 22567.2109375, 344060.27734375, 
    5342799.97949, 84045005.15283, 1334792319.64120, 21356125571.881271 };
  for ( size_t i = 0; i < correctOut.size(); ++i ){ 
    REQUIRE( correctOut[i]  == Approx( p[i] ).epsilon(1e-6) );
  }


}
 
