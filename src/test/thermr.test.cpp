#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"


TEST_CASE( "thermr" ){
  GIVEN( "" ){
    WHEN( "" ){
      REQUIRE( 0 == thermr() );
      
    } // WHEN
  } // GIVEN
} // TEST CASE
