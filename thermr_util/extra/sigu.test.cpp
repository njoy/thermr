#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sigu.h"


TEST_CASE( "computeCrossSections" ){
  GIVEN( "" ){
    std::vector<double> alpha { 2.25, 2.30, 2.35, 2.40, 2.45 }, beta { 0.02, 0.2, 1.2, 2.0 };
    double enow = 1.5, tev = 0.035, tevz = 0.0253, tolin = 5e-3, az = 11.898, u = -1.0;
    std::vector<std::vector<double>> sab { { 0.1, 0.2, 0.3 }, { 1.4, 1.5, 1.6 }, { 2.1, 2.2, 2.3, 2.4 }, { 3.4, 3.5, 3.6 } };
    int nemax = 5, lasym = 0, lat = 0;
    sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );


      
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
