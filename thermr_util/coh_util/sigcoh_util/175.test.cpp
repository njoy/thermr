#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "175.h"


TEST_CASE( "ending things with 175 / 180 / 185" ){
  GIVEN( "input values" ){
  std::vector<double> fl, correctFl, wrk (10000,0.0);
  std::vector<double> wrkVals { 1.1, 1.2, 2.3, 2.4, 3.5, 3.6, 4.7, 4.8, 5.9, 16.1, 16.2, 16.3, 16.4, 26.5, 26.6, 26.7, 26.8, 36.9, 37.0 };
  for ( size_t i = 0; i < wrkVals.size(); ++i ){ wrk[i] = wrkVals[i]; }
  std::vector<double> correctVals { 1.1, 1.2, 2.3, 2.4, 3.5, 3.6, 4.7, 4.8, 5.9, 16.1, 1.2e19, 16.1, 0, 0, 0, 0, 0, 0 };

  int k = 5;
  double recon = 5.18e-20, ulim = 1.2e19;
  fl = end_175_180_185( wrk, k, recon, ulim );
  //std::cout << fl[0] << std::endl;
  for ( size_t i = 0; i < correctVals.size(); ++i ){ REQUIRE( correctVals[i] == Approx( fl[i] ).epsilon(1e-6) ); }
    THEN( "correct" ){
      REQUIRE( true );
    } 
  } // GIVEN
GIVEN( "input values" ){
  std::vector<double> fl, correctFl, wrk (10000,0.0);
  std::vector<double> wrkVals { 1.1, 1.2, 2.3, 2.4, 3.5, 3.6, 4.7, 4.8, 5.9, 16.1, 16.2, 16.3, 16.4, 26.5, 26.6, 26.7, 26.8, 36.9, 37.0 };
  for ( size_t i = 0; i < wrkVals.size(); ++i ){ wrk[i] = wrkVals[i]; }
  std::vector<double> correctVals { 1.1, 1.2, 2.3, 2.4, 3.5, 3.6, 1.2e19, 3.6,
     3.5999999046325684,        2.0000000000000000,
   2.0000000000000000,        9.9999999999057137E-004,   2.0000000000000000E-003,
   3.5000000000000001E-003,   5.0000000000000001E-003,   7.4999999999999997E-003,
   1.0000000000000000E-002,   1.2999999999999999E-002,   1.6500000000000001E-002};

  int k = 3;
  double recon = 5.18e-20, ulim = 1.2e19;
  fl = end_175_180_185( wrk, k, recon, ulim );
  //std::cout << fl[0] << std::endl;
  for ( size_t i = 0; i < correctVals.size(); ++i ){ 
    std::cout << fl[i] << "     " << correctVals[i] << std::endl;
    REQUIRE( correctVals[i] == Approx( fl[i] ).epsilon(1e-6) );
  }
    THEN( "correct" ){
      REQUIRE( true );
    } 
  } // GIVEN

} // TEST CASE
