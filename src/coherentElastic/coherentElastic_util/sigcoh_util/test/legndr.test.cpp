#include "catch.hpp"
#include "coherentElastic/coherentElastic_util/sigcoh_util/legndr.h"


TEST_CASE( "legndr" ){
  std::vector<double> p { 0.1, 0.2, 0.3, 0.4, 0.5 }, correctP(5);
  double x; int np;
  GIVEN( "np is small" ){
    WHEN( "np is zero" ){
      x = 0.5; np = 0;
      correctP = { 1.0, 0.5, 0.3, 0.4, 0.5 };
      legndr( x, p, np );
      THEN( "only the first two np values are changed" ){
        for ( size_t i = 0; i < p.size(); ++i ){
          REQUIRE( correctP[i] == Approx( p[i] ).epsilon(1e-6) );
        }
      } // THEN
    } // WHEN

    WHEN( "np is three" ){
      x = 0.5; np = 3;
      correctP = { 1.0, 0.5, -0.125, -0.4375, 0.5 };
      legndr( x, p, np );
      THEN( "only the first four np values are changed" ){
        for ( size_t i = 0; i < p.size(); ++i ){
          REQUIRE( p[i] == Approx( correctP[i] ).epsilon(1e-6) );
        }
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "np is large" ){
    np = 4;
    WHEN( "x is in a valid range" ){
      x = 0.5; np = 4;
      correctP = { 1.0, 0.5, -0.125, -0.4375, -0.289062 };
      legndr( x, p, np );
      for ( size_t i = 0; i < p.size(); ++i ){
        REQUIRE( p[i] == Approx( correctP[i] ).epsilon(1e-6) );
      }
    } // WHEN
    /*
    WHEN( "x is not in a valid range" ){
      AND_WHEN( "x < -1" ){
	x = -1.1;
        REQUIRE_THROWS( legndr( x, p, np ) );
      } // AND WHEN
      AND_WHEN( "x > 1" ){
	x = 1.1;
        REQUIRE_THROWS( legndr( x, p, np ) );
      } // AND WHEN
    } // WHEN
    */
  } // GIVEN
} // TEST CASE
