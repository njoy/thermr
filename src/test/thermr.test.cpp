#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"


TEST_CASE( "thermr" ){
  GIVEN( "" ){
    WHEN( "" ){
        std::string fileName = "/Users/ameliajo/thermr/src/test/tape24";
        int mat = 101;
        int iform = 0;
        thermr(mat,fileName,iform);

        REQUIRE( 0 == 0 );


    } // WHEN
  } // GIVEN
} // TEST CASE
