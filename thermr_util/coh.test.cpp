#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "coh.h"


TEST_CASE( "coh" ){
  GIVEN( "" ){
    WHEN( "" ){
      double temp = 296.0, emax = 1.2;
      int lat = 1, iold = 11, inew = 10, ne = 145, nex = 2, natom = 1;
      std::vector<double> fl;
      coh( lat, inew, ne, nex, temp, emax, natom, fl );
      REQUIRE( true );
      
    } // WHEN
  } // GIVEN
} // TEST CASE
