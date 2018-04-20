#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "sig.h"


TEST_CASE( "sig" ){
  int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc;

  double e = 1.0e-5, ep = 0.0, u = -1.0, tev = 2.5e-2, bbm = 0.0, az = 11.9,
    tevz = 2.53e-2, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2;

  double sigVal1, sigVal2;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
    beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
  for ( int i = 0; i < nalpha; ++i ){
    for ( int j = 0; j < nbeta; ++j ){
      sab[i][j] = 0.01*((j+1) + 0.1*(i+1));

    } 
  } 


  GIVEN( "invalid input requested (iinc != 1,2)" ){
    iinc = 0;
    THEN( "throw exception" ){
      REQUIRE_THROWS( sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
        bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc ) );
    } // THEN
  } // GIVEN

  GIVEN( "free gas option selected (iinc = 1)" ){
    iinc = 1;

    WHEN( "final neutron energy E' is zero (ep = 0)" ){
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 0.0 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN

    ep = 1.2e-6;

    WHEN( "temperature is extremely small" ){
      tev = 1.5e-8;
      sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
        bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
      THEN( "alpha, beta --> big, so exponential in Eq. 229 dies" ){
        REQUIRE( 0.0 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "nonextreme values are provided" ){
      sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
        bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
      THEN( "nonzero sig output is provided" ){
        REQUIRE( 1384.063418 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN ("alpha, beta values outside of provided range" ){
    iinc = 2; tev = 1.5e-1; tevz = 2.2e-4;
    WHEN( "this is because of a > alpha(nalpha)" ){
      e = 1.0e-2, ep = 1.2e-2;
      sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
        bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );

      e = 1.0e-3, ep = 1.2e-3; alpha = { 0.1, 0.2, 0.3, 0.5, 0.8 };
      sigVal2 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
        bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );

      THEN( "SCT Approximation is used" ){
        REQUIRE( 55.838635 == Approx( sigVal1 ).epsilon(1e-6) );
        REQUIRE( 179.2164755 == Approx( sigVal2 ).epsilon(1e-6) );
      } // THEN
    } // WHEN

    WHEN( "this is because of being out of range with beta" ){
      e = 1.0e-5;

      AND_WHEN( "lasym != 1" ){ 
        e = 1.0e-4, ep = 5.2e-3, tev = 1.5e-1;
        lasym = 0;
        sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 795.18855771477308 == Approx( sigVal1 ).epsilon(1e-6) );
      } // AND WHEN

      AND_WHEN( "lasym == 1" ){
        AND_WHEN( "bbm is smaller than min value in beta vector" ){ 
          ep = 1.2e-6, tev = 1.5e-8;
          alpha = { 0.1, 0.2, 0.3, 0.5, 0.8 };
          lasym = 1;
          sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
            bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
          REQUIRE( 883.34692414 == Approx( sigVal1 ).epsilon(1e-6) );
        } // AND WHEN
        AND_WHEN( "bbm is greater than max value in beta vector" ){ 
            ep = 1.2e-1;
            alpha = { 10.1, 20.2, 30.3, 40.5, 50.8 };

          AND_WHEN( "-arg < sabflg" ){ 
             tev = 1.5e-8;
            sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
              bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( 0.0 == Approx( sigVal1 ).epsilon(1e-6) );
          } // AND WHEN

          AND_WHEN( "-arg > sabflg" ){ 
            sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
              bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );

            e = 1.0e-4;
            sigVal2 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
              bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
            THEN( "SCT Approximation is used" ){
              REQUIRE( 12.9236152 == Approx( sigVal1 ).epsilon(1e-6) );
              REQUIRE( 5.011738726 == Approx( sigVal2 ).epsilon(1e-6) );
            } // THEN
          } // AND WHEN
        } // AND WHEN
      } // AND WHEN
    } // WHEN
  } // GIVEN

  GIVEN( "150" ){
    iinc = 2;
    e = 1.0e-6, ep = 1.2e-4, u = -1.0, tev = 1.5e-1, tevz = 2.2e-4;
    WHEN( "final neutron energy E' is zero (ep = 0)" ){
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 201.87960468 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
    tev = 1.5e-4;
    WHEN( "150 to 155 to 160" ){
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 135829.64964 == Approx( sigVal1 ).epsilon(1e-6) );

        cliq = 1.0;
        sigVal2 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 135829.64964 == Approx( sigVal2 ).epsilon(1e-6) );

      } // THEN
    } // WHEN
    WHEN( "straight to 160" ){
    cliq = 1.0;
      ep = 3.1e-5;
      THEN( "output cross section is 0 (see Eq. 225) " ){
        sigVal1 = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 248176.610043 == Approx( sigVal1 ).epsilon(1e-6) );
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE

