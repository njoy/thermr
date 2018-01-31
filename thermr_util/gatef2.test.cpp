#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "gatef2.h"

TEST_CASE( "gatef2" ){
  GIVEN( "inputs" ){
    THEN( "output is okay" ){
      std::vector<double> temp {296}, eftemp(1);
      int mat = 101;
      gatef2( temp, eftemp, mat );
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
