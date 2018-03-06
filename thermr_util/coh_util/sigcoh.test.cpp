#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sigcoh.h"


TEST_CASE( "sigcoh" ){
  GIVEN( "" ){
    double e = 0, enext = 0, temp = 296, emax = 1.2;
    int nl = 1, lat = 1, natom = 1;
    std::vector<double> s (6, 0.0);
    sigcoh( e, enext, s, nl, lat, temp, emax, natom );
    WHEN( "" ){
      REQUIRE( true );
      
    } // WHEN
  } // GIVEN
} // TEST CASE
