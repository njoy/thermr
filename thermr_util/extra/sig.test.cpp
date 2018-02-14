#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sig.h"


TEST_CASE( "compute differential scattering ccross sections" ){
  GIVEN( "" ){
    double e = 1.5, ep = 0, u = -1, tev =  3.5E-2, tevz = 2.53E-2, 
      az = 11.898, az2 = 0;
    std::vector<double> alpha { 2.25, 2.30, 2.35, 2.40 }, beta { 0.02, 0.2, 1.2 };
    int lat = 0, iinc = 2, lasym = 0;
    double cliq = 0, sb = 0, sb2 = 0, teff = 0, teff2 = 0;
    std::vector<std::vector<double>> sab { {0.1, 1.4, 2.1, 3.4}, {0.2, 1.5, 2.2, 3.5}, {0.3, 1.6, 2.3, 3.6 } };


    //WHEN( "a > alpha[nalpha-1]" ){
    //  double sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
    //  REQUIRE( 0 == Approx( sigVal ).epsilon(1e-6) ); 
    //} // WHEN

    WHEN( "temperature is high" ){
      tev = 2.5;
      double sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
    } // WHEN


    REQUIRE( true );


  } // GIVEN
} // TEST CASE
