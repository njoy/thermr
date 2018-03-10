#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sigl.h"


TEST_CASE( "sigl" ){
  int nlin, nlmax, lat, lasym, iinc, nbin;
  /*
  GIVEN( "110 and 120 oscillations, no 130 (sum = 0)" ){
    double e, ep, tev, tolin, az, az2, teff, teff2, tevz, cliq, sb, sb2;
    std::vector<double> alpha(4), beta(3), s(65, 1.2);
    std::vector<std::vector<double>> sab { {0.1, 0.2, 0.3}, {1.4, 1.5, 1.6}, 
        {2.1, 2.2, 2.3}, {3.4, 3.5, 3.6} };


    nlin = 2, nlmax = 65, lat = 0, lasym = 0, iinc = 2, nbin = 8;
    e = 1.5, ep = 2.3, tev = 3.5e-2, tolin = 5e-3, az = 11.898, 
    az2 = 0.0, teff = 300, teff2 = 350, tevz = 2.53e-2, cliq = 0.0, 
    sb = 0.0, sb2 = 0.0;
    alpha = { 2.25, 2.30, 2.35, 2.40 };
    beta = { 0.02, 0.2, 1.2 };

    sigl( e, ep, tev, alpha, beta, s, sab, tolin, nlin, az, az2, teff, 
        teff2, lat,tevz, lasym, iinc, cliq, sb, sb2, nbin );

    REQUIRE( 0.0 == Approx( s[0] ).epsilon(1e-6) ); 

    THEN( "output s vector is correct" ){
      for ( size_t i = 1; i < s.size(); ++i ){
        REQUIRE( 1.2 == Approx( s[i] ).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  std::cout << "\n" << std::endl;
  GIVEN( "110 and 120 oscillations, and 130 (sum != 0)" ){
    double e, ep, tev, tolin, az, az2, teff, teff2, tevz, cliq, sb, sb2;
    std::vector<double> alpha(4), beta(3), s(65, 1.2);
    std::vector<std::vector<double>> sab { {0.1, 0.2, 0.3}, {1.4, 1.5, 1.6}, 
      {2.1, 2.2, 2.3}, {3.4, 3.5, 3.6} };


    nlin = 2, nlmax = 65, lat = 0, lasym = 0, iinc = 2, nbin = 8;
    e = 1.5, ep = 2.3, tev = 3.5e-2, tolin = 5e-3, az = 11.898, 
    az2 = 0.0, teff = 300, teff2 = 350, tevz = 2.53e-2, cliq = 0.0, 
    sb = 200, sb2 = 0.0;
    alpha = { 2.25, 2.30, 2.35, 2.40 };
    beta = { 0.02, 0.2, 1.2 };

    sigl( e, ep, tev, alpha, beta, s, sab, tolin, nlin, az, az2, teff, 
        teff2, lat,tevz, lasym, iinc, cliq, sb, sb2, nbin );

    REQUIRE( 1.328733E-9 == Approx( s[0] ).epsilon(1e-6) ); 
    REQUIRE( 0.358778963 == Approx( s[1] ).epsilon(1e-6) ); 

    THEN( "output s vector is correct" ){
      for ( size_t i = 2; i < s.size(); ++i ){
      REQUIRE( 1.2 == Approx( s[i] ).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  */
  /*
  GIVEN( "oscillations down to 190, where 190 goes back to 160 for a bit, hitting 170 as well" ){
    double e, ep, tev, tolin, az, az2, teff, teff2, tevz, cliq, sb, sb2;
    std::vector<double> alpha(4), beta(3), s(65, 1.2);
    std::vector<std::vector<double>> sab { {0.1, 0.2, 0.3}, {1.4, 1.5, 1.6}, 
        {2.1, 2.2, 2.3}, {3.4, 3.5, 3.6} };


    nlin = 4, nlmax = 65, lat = 0, lasym = 0, iinc = 2, nbin = 8;
    e = 1.5, ep = 2.3, tev = 3.5e-2, tolin = 5e-3, az = 11.898, 
    az2 = 0.0, teff = 300, teff2 = 350, tevz = 2.53e-2, cliq = 0.0, 
    sb = 200, sb2 = 0.0;
    alpha = { 2.25, 2.30, 2.35, 2.40 };
    beta = { 0.02, 0.2, 1.2 };



    sigl( e, ep, tev, alpha, beta, s, sab, tolin, nlin, az, az2, teff, 
        teff2, lat,tevz, lasym, iinc, cliq, sb, sb2, nbin );

    REQUIRE( 1.3287332e-9 == Approx( s[0] ).epsilon(1e-6) ); 
    REQUIRE( 0.3587789636 == Approx( s[1] ).epsilon(1e-6) ); 
    REQUIRE( 0.1070770554 == Approx( s[2] ).epsilon(1e-6) ); 
    REQUIRE( 0.1463496778 == Approx( s[3] ).epsilon(1e-6) ); 

    THEN( "output s vector is correct" ){
      for ( size_t i = 4; i < s.size(); ++i ){
        REQUIRE( 1.2 == Approx( s[i] ).epsilon(1e-6) );
      }
    } // THEN
  } // GIVEN
  */

  GIVEN( "170 --> 190" ){
    double e, ep, tev, tolin, az, az2, teff, teff2, tevz, cliq, sb, sb2;
    std::vector<double> alpha(4), beta(3), s(65, 1.2);
    std::vector<std::vector<double>> sab { {0.1, 0.2, 0.3}, {1.4, 1.5, 1.6}, 
        {2.1, 2.2, 2.3}, {3.4, 3.5, 3.6} };


    nlin = 34, nlmax = 65, lat = 0, lasym = 0, iinc = 2, nbin = 8;
    e = 1.5, ep = 2.3, tev = 3.5e-2, tolin = 5e-3, az = 11.898, 
    az2 = 0.0, teff = 300, teff2 = 350, tevz = 2.53e-2, cliq = 0.0, 
    sb = 200, sb2 = 0.0;
    alpha = { 2.25, 2.30, 2.35, 2.40 };
    beta = { 0.02, 0.2, 1.2 };



    sigl( e, ep, tev, alpha, beta, s, sab, tolin, nlin, az, az2, teff, 
        teff2, lat,tevz, lasym, iinc, cliq, sb, sb2, nbin );

    std::vector<double> correctS { 1.3287332e-9, 0.358778963, 0.182154203, 
      7.9932519e-2, 1.5802056e-2, -7.1904468e-3, -1.6204572e-2, -3.8658941e-4, 
      -2.4271956e-3, 1.1677628e-2, -4.5363482e-3, 5.6039616e-3, -1.0806403e-2, 
      4.1381041e-3, -5.2553740e-3, 2.0454716e-3, -7.3896342e-4, -7.1686680e-3, 
      3.7516430e-3, -1.4007032e-2, 8.0315285e-3, -1.4431481e-2, 3.7936759e-3, 
      -8.4267224e-3, -8.7385347e-3, 2.6585217e-3, -1.8869074e-2, 7.6194707e-3, 
      -1.7909766e-2, -1.8522299e-3, -6.1935417e-3, -1.4653166e-2, 1.0003996e-3, 
      -1.3574744e-2 };

    for ( size_t i = 0; i < 34; ++i ){ 
      REQUIRE( correctS[i] == Approx( s[i] ).epsilon(1e-6) ); 
    } 
    THEN( "output s vector is correct" ){
      for ( size_t i = 35; i < s.size(); ++i ){
        REQUIRE( 1.2 == Approx( s[i] ).epsilon(1e-6) );
      }
    } // THEN


    //REQUIRE( 1.3287332e-9 == Approx( s[0] ).epsilon(1e-6) ); 
    //REQUIRE( 0.3587789636 == Approx( s[1] ).epsilon(1e-6) ); 
    //REQUIRE( 0.1070770554 == Approx( s[2] ).epsilon(1e-6) ); 
    //REQUIRE( 0.1463496778 == Approx( s[3] ).epsilon(1e-6) ); 

  } // GIVEN



} // TEST CASE
