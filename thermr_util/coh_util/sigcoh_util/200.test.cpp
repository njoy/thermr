#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "200.h"


TEST_CASE( "read bragg parameters" ){
  GIVEN( "" ){
    std::vector<double> fl { 1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9, 10.0,
    11.11, 12.12, 13.13, 14.14, 15.15, 16.16, 0, 0, 0, 0, 0, 0 };
    int l = 1;
    double econ = 1.93038465E+19;
    readBraggParameters( fl, l, econ );
    std::vector<double> correctFl { 1.9110807E20, 10, 2.1446572E20, 2.12, 
      2.5345950E20, 2.02, 2.924532E20, 2.02, 0, -16.16, 0, 0, 13.13, 14.14, 
      15.15, 16.16 };
    for ( size_t i = 0; i < correctFl.size(); ++i ){ 
      REQUIRE( correctFl[i] == Approx( fl[i] ).epsilon(1e-6) ); 
    }
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
  GIVEN( "" ){
    std::vector<double> fl { 11.1, 12.2, 13.3, 14.4, 15.5, 16.6, 17.7, 18.8, 
      19.9, 120.0, 121.11, 122.12, 123.13, 124.14, 125.15, 126.16, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0 };
    int l = 1;
    double econ = 1.93038465E+19;
    readBraggParameters( fl, l, econ );
    std::vector<double> correctFl { 3.841465E20, 120, 2.337888E21, 2.12, 
      2.3768825E21, 2.02, 2.4158764E21, 2.02, 0, -126.16, 0, 0, 0, 0, 0, 0 };
    for ( size_t i = 0; i < correctFl.size(); ++i ){ 
      REQUIRE( correctFl[i] == Approx( fl[i] ).epsilon(1e-6) ); 
    }

  } // GIVEN
} // TEST CASE
