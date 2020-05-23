#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"


TEST_CASE( "thermr" ){
  GIVEN( "" ){
    WHEN( "" ){
        //std::string fileName = "/Users/ameliajo/thermr/src/test/tape24";
        std::string fileName = "/Users/ameliajo/thermr/src/test/tape24_smaller_alpha_beta_grid";
        int mat = 101;
        int iform = 0;
        int iinc = 2;
        int nbin = 8;
        int natom = 1;
        std::vector<double> temperatures {296.0};
        double tol = 0.05, emax = 0.625;
        thermr(mat,fileName,iform,iinc,nbin,temperatures,tol,emax,natom);

        REQUIRE( 0 == 0 );


    } // WHEN
  } // GIVEN
} // TEST CASE
