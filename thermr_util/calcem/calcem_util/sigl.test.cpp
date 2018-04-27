#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "sigl.h"


TEST_CASE( "sigl" ){
  GIVEN( "inputs" ){
    int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 2, nlmax = 65, 
        nlin = 10;

    std::cout << std::setprecision(10) ;
    double e = 1.0e-6, ep = 1.2e-4, tev = 1.5e-4, bbm = 0.0, az = 11.9,
      tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2;

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
      s(65,0.0);

    std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
    for ( int i = 0; i < nalpha; ++i ){
      for ( int j = 0; j < nbeta; ++j ){
        sab[i][j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 


    sigl( nlin, nalpha, nbeta, nlmax, e, ep, tev, alpha, beta, sab, s, tolin, az, tevz, iinc, lat, bbm, lasym, az2, teff2, cliq, sb, sb2, teff );

    for ( int i = 0; i < 10; ++i ){std::cout << s[i] << std::endl; }




    REQUIRE( true );
  } // GIVEN
} // TEST CASE

