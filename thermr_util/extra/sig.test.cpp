#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sig.h"


TEST_CASE( "compute differential scattering ccross sections" ){
  GIVEN( "" ){
    double e = 1.5, ep = 0, u = -1, tev =  3.5E-2, tevz = 2.53E-2, 
      az = 11.898, az2 = 0, sigVal;
    std::vector<double> alpha { 2.25, 2.30, 2.35, 2.40 }, beta { 0.02, 0.2, 1.2 };
    int lat = 0, iinc = 2, lasym = 0;
    double cliq = 0, sb = 0, sb2 = 0, teff = 0, teff2 = 0;
    std::vector<std::vector<double>> sab { {0.1, 0.2, 0.3}, {1.4, 1.5, 1.6}, {2.1, 2.2, 2.3}, {3.4, 3.5, 3.6} };

/*

    WHEN( "temperature is high" ){
      AND_WHEN( "energy e is relatively low" ){
        // 150A 155 160
        THEN( "sigVal < sigmin, so it's set to zero" ){
          tev = 2.5;
          sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
          REQUIRE( 0.0 == Approx( sigVal ).epsilon(1e-6) );
        } // THEN
      } // AND WHEN

      AND_WHEN( "energy is higher" ){
        // 150A 155 160
        THEN( "sigVal is not set to zero" ){
          e = 15;
          ep = 14.999;
          tev = 2.5;
          sb = 5;
          sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
          REQUIRE( 1.45409235E-3 == Approx( sigVal ).epsilon(1e-6) );

          ep = 12.05;
          sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
          REQUIRE( 6.30183328E-6 == Approx( sigVal ).epsilon(1e-6) );

          tev = 3.0; 
          sb = 15;
          sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
          REQUIRE( 3.0756984E-10 == Approx( sigVal ).epsilon(1e-6) );
          
        } // THEN
      } // AND WHEN
    } // WHEN

    // 170
    WHEN( "straight to 170" ){
      e = 15;
      ep = 10;
      tev = 1.0;
      sb = 15;
      teff = 1.5;
      sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
      REQUIRE( 0.6722681752 == Approx( sigVal ).epsilon(1e-6) );
      //std::cout << sigVal << std::endl;
      //
      sb = 10;
      teff = 80.5;
      sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
      REQUIRE( 6.289982889E-2 == Approx( sigVal ).epsilon(1e-6) );

      teff = 280.5;
      sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
      REQUIRE( 3.370886372E-2 == Approx( sigVal ).epsilon(1e-6) );
    } // WHEN

    */

    WHEN( "iinc != 2 " ){
      AND_WHEN( "iinc == 1" ){
        iinc = 1;
        e = 12; 
        ep = 9;
        tev = 1.0;
        sb = 7;
        teff = 320.5;
        sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
        REQUIRE( 0.4478383035 == Approx( sigVal ).epsilon(1e-6) );
      //  std::cout << sigVal << std::endl;
        e = 17; 
        ep = 16;
        tev = 0.2;
        sb = 13;
        teff = 320.5;
        sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
        REQUIRE( 1.60255718E-2 == Approx( sigVal ).epsilon(1e-6) );


      } // AND WHEN
      AND_WHEN( "iinc != 1" ){
        iinc = 3;
        sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
        REQUIRE( 0.0 == Approx( sigVal ).epsilon(1e-6) );
      //  std::cout << sigVal << std::endl;

      } // AND WHEN

    } // WHEN

    // 170C
    WHEN( "lasym = 1" ){
      e = 17;
      ep = 16;
      tev = 5.2;
      sb = 13;
      teff = 320.5;
      lasym = 1;
      sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
      REQUIRE( 4.20709698E-2 == Approx( sigVal ).epsilon(1e-6) );
    } // WHEN

    // 170C
    WHEN( "lasym = 1 and 170C" ){
      e = 17;
      ep = 22;
      tev = 5.2;
      sb = 13;
      teff = 20.5;
      lasym = 1;
      beta = { 0.01, 0.1, 0.6 };
      sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
      REQUIRE( 6.86308387285E-2== Approx( sigVal ).epsilon(1e-6) );

    } // WHEN
    // 170D
    WHEN( "lasym = 0 and 170C" ){
      e = 17;
      ep = 22;
      tev = 5.2;
      sb = 13;
      teff = 20.5;
      lasym = 1;
      beta = { 0.01, 0.1, 0.6 };
      sigVal = sig( e, ep, u, tev, tevz, alpha, beta, sab, az, az2, lat, iinc, lasym, cliq, sb, sb2, teff, teff2 );
      REQUIRE( 6.863083872E-2 == Approx( sigVal ).epsilon(1e-6) );

    } // WHEN

  } // GIVEN
} // TEST CASE
