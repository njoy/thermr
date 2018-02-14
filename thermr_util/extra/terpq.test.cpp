#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "terpq.h"


TEST_CASE( "terpq" ){
  double out;
  GIVEN( "desired x value is below provided range (x < x1)" ){
    WHEN( "x < x1 and y1 <= y2" ){
      THEN( "terp1 is called, interpolation code 3 is used" ){
        //                x1    x2    x3         x       y1     y2     y3    y
        double in[8] = { 2.25, 2.30, 2.35, 5.0428643E-2, 0.10, 1.40, 2.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( -224.55032636 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "x < x1 and y1 > y2" ){
      THEN( "y = y1 is returned" ){
        //                x1    x2    x3         x       y1     y2     y3    y
        double in[8] = { 2.25, 2.30, 2.35, 5.0428643E-2, 4.10, 3.40, 2.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( 4.1 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
  GIVEN( "desired x value is above provided range (x > x3)" ){

    WHEN( "x > x3 and y3 <= y2" ){
      THEN( "terp1 is called, interpolation code 2 is used" ){
        //                x1    x2    x3    x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 15.0, 4.10, 3.40, 2.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( 1.1553953878936571 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "x > x3 and y3 > y2" ){
      THEN( "y is set to y3 and returned" ){
        //                x1    x2    x3    x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 15.0, 1.10, 3.40, 8.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( 8.10 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
  GIVEN( "desired x value is within range (x1 < x < x3)" ){
    WHEN( "big steps and x >= x2" ){
      THEN( "terp1 is used, interpolating between x2 and x3" ){
        //                x1    x2    x3   x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 5.0, 1.10, 3.40, 8.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( 4.7525182 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "big steps and x < x2" ){
      THEN( "terp1 is used, interpolating between x1 and x2" ){
        //                x1    x2    x3   x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 2.0, 1.10, 3.40, 8.10, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( 2.2205129 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
    WHEN( "small steps" ){
      THEN( "terp1 is used, interpolating between x1 and x2" ){
        //                x1    x2    x3   x    y1    y2     y3    y
        double in[8] = { 1.05, 3.00, 9.95, 2.0, 1.10, 1.42, 1.73, 0.0 };
        out = terpq( in[0], in[1], in[2], in[3], in[4], in[5], in[6], in[7] );
        REQUIRE( 1.268652857 == Approx( out ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
 
