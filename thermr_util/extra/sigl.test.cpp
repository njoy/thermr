#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sigl.h"


TEST_CASE( "sigl" ){
  GIVEN( "" ){
    double e = 1.5, ep = 2.3, tev = 3.5e-2;
    std::vector<double> alpha { 2.25, 2.30, 2.35, 2.40 }, 
      beta { 0.02, 0.2, 1.2 }, s ( 65, 1.2 );
    std::vector<std::vector<double>> sab { {0.1, 0.2, 0.3}, {1.4, 1.5, 1.6}, 
      {2.1, 2.2, 2.3}, {3.4, 3.5, 3.6} };
    double tolin = 5e-3, az = 11.898, az2 = 0.0, teff = 300, teff2 = 350, tevz = 2.53e-2, cliq = 0.0, sb = 0.0, sb2 = 0.0;

    int nlin = 2, nlmax = 65, lat = 0, lasym = 0, iinc = 2;
    sigl( e, ep, tev, alpha, beta, s, sab, tolin, nlin, az, az2, teff, teff2, lat,tevz, lasym, iinc, cliq, sb, sb2 );
    //std::cout << "   "<< std::endl;
    //for ( auto entry : s ) { std::cout << entry << std::endl; }
    REQUIRE( 0.0 == Approx( s[0] ).epsilon(1e-6) ); 
    for ( size_t i = 1; i < s.size(); ++i ){
      REQUIRE( 1.2 == Approx( s[i] ).epsilon(1e-6) );
    }
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
