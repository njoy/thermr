#define CATCH_CONFIG_MAIN 
#include "../../../catch.hpp"
#include "313.h"

void checkVecHasZeros(std::vector<double>& v, int index_0, int index_1 ){
  for ( size_t i = index_0; i < index_1; ++i ){
    REQUIRE( 0.0 == Approx(v[i]).epsilon(1e-6) );
  }
}
void checkVecHasLin(std::vector<double>& v, int index_0, int index_1 ){
  for ( size_t i = index_0; i < index_1; ++i ){
    REQUIRE( (i+1) == Approx(v[i]).epsilon(1e-6) );
  }
}

void initializeLin(std::vector<double>& v){
  for ( size_t i = 0; i < v.size(); ++i ){ 
    v[i] = i + 1;
  }
}


TEST_CASE( "313" ){
  GIVEN( "inputs" ){
    THEN( "single iteration is necessary" ){
    int jbeta = 10, lat = 1, iskip = 0;
    double ep = 2.5, enow = 1.3, tev = 2.5e-2, tevz = 2.53e-2;
    std::vector<double> beta (80), x(20);
    for ( size_t i = 0; i < 80; ++i ){ beta[i] = 0.20*i + 0.05; }
    for ( size_t i = 0; i < 20; ++i ){ x[i]    = 0.03*i + 0.03; }
    
    do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
    REQUIRE( 10 == jbeta );
    REQUIRE( 1.346805 == Approx(ep).epsilon(1e-6) );
    REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
    REQUIRE( 0 == iskip );  
    } // THEN
    
  } // GIVEN
  GIVEN( "inputs" ){
    THEN( "multiple iterations are necessary" ){
    int jbeta = 10, lat = 1, iskip = 0;
    double ep = 2.5, enow = 1.3, tev = 2.5e-2, tevz = 2.53e-2;
    std::vector<double> beta (80), x(20);
    for ( size_t i = 0; i < 80; ++i ){ beta[i] = 0.20*i + 0.05; }
    for ( size_t i = 0; i < 20; ++i ){ x[i]    = 0.7*i + 0.7; }
    
    do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip );
    REQUIRE( 21 == jbeta );
    REQUIRE( 1.402465 == Approx(ep).epsilon(1e-6) );
    REQUIRE( 1.3 == Approx(enow).epsilon(1e-6) );
    REQUIRE( 0 == iskip );  
    } // THEN
    
  } // GIVEN

} // TEST CASE
