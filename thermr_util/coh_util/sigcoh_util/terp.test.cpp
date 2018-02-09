#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "terp.h"


TEST_CASE( "do260" ){
  std::vector<double> x, y, argVals, trueVals;
  double arg;
  GIVEN( "x, y vectors of equal length" ){
  x = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 },
  y = { 2.1997, 2.7448, 3.2912, 3.851, 4.421, 4.9969, 6.1624, 7.3387, 
          9.6287, 11.992 };

    WHEN( "arg inputs that are not in x" ){
    AND_WHEN( "inputs are in range" ){
      
      argVals = { 312, 385, 423, 945, 1345, 1995 };
      trueVals = { 2.283561, 2.666179, 2.865350, 5.601333, 7.697872, 11.10474 };
      int l = 1, il = 2;
    THEN( "outputs" ){
      for ( size_t i = 0; i < argVals.size(); ++i ){
        REQUIRE( trueVals[i] == Approx( do260( x, y, argVals[i], l, il ) ).epsilon(1e-6) );
      }
    } // THEN
    } // AND WHEN

    AND_WHEN( "inputs are not in range" ){
      argVals = { 0, 5, 50, 185, 204, 2001, 2450, 3000, 10000 };
      trueVals = { 0.6482615, 0.6744682, 0.9103288, 1.617910, 1.717496, 11.13619, 13.48955, 16.3723, 53.06172 };
      int l = 1, il = 2;
    THEN( "outputs" ){
      for ( size_t i = 0; i < argVals.size(); ++i ){
        REQUIRE( trueVals[i] == Approx( do260( x, y, argVals[i], l, il ) ).epsilon(1e-6) );
      }
    } // THEN
    } // AND WHEN
  } //  WHEN 

  WHEN( "inputs that are in x" ){
      argVals = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 };
      trueVals = { 2.1997, 2.7448, 3.268934, 3.793069, 4.317203, 4.841338, 5.889607, 6.937876, 9.034415, 11.13095 };
      arg = 296;
      int l = 1, il = 2;
      auto sum = do260( x, y, arg, l, il );
    THEN( "outputs" ){
      for ( size_t i = 0; i < argVals.size(); ++i ){
        REQUIRE( trueVals[i] == Approx( do260( x, y, argVals[i], l, il ) ).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  } // GIVEN
} // TEST CASE

TEST_CASE( "terp" ){
  GIVEN( "" ){
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
