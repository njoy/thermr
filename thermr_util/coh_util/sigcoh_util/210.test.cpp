#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "210.h"


TEST_CASE( "computeCrossSections" ){
  GIVEN( "" ){
    std::vector<double> s(10, 0), fl { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1 };
    double e = 0.015;
    int nl = 10, k = 10, nw = 10;
    double emax = 0.625, scon = 0.2, recon = 0.5;
      
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
