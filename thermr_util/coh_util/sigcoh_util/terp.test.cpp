#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "terp.h"


TEST_CASE( "do260" ){
  std::vector<double> x, y, argVals, trueVals;
  double arg;
  GIVEN( "inputs" ){
    THEN( "outputs" ){
      x = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 },
      y = { 2.1997, 2.7448, 3.2912, 3.851, 4.421, 4.9969, 6.1624, 7.3387, 
          9.6287, 11.992 };
      
      argVals = { 204, 296, 312, 1000 };
      trueVals = { 1.7174961538461542, 2.1997, 2.2835615384615386, 5.8896076923076954 };
      arg = 296;
      int l = 1, il = 2;
      auto sum = do260( x, y, arg, l, il );
      for ( size_t i = 0; i < argVals.size(); ++i ){
        REQUIRE( trueVals[i] == Approx( do260( x, y, argVals[i], l, il ) ).epsilon(1e-6) );
      }
      //REQUIRE( 2.1997 == Approx( sum ).epsilon(1e-6) );
    } // THEN
  } // GIVEN
} // TEST CASE

TEST_CASE( "terp" ){
  GIVEN( "" ){
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
