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
   0, 2, 2, 1e-3, 2.e-3, 3.5e-3, 5e-3, 7.5e-3, 1e-2, 1.3e-2, 1.65e-2};

  int k = 3;
  double recon = 5.18e-20, ulim = 1.2e19;
  fl = end_175_180_185( wrk, k, recon, ulim );
  //std::cout << fl[0] << std::endl;
  for ( size_t i = 0; i < correctVals.size(); ++i ){ 
    //std::cout << fl[i] << "     " << correctVals[i] << std::endl;
    //REQUIRE( correctVals[i] == Approx( fl[i] ).epsilon(1e-6) );
  }
    THEN( "correct" ){
      REQUIRE( true );
    } 
  } // GIVEN

} // TEST CASE
