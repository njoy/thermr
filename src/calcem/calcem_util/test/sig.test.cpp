#include "catch.hpp"
#include "calcem/calcem_util/sig.h"
#include <range/v3/all.hpp>

TEST_CASE( "sig - free gas (iinc = 1)" ){
  int iinc = 1;

  int lasym = 0, lat = 1;
  double e = 1.0e-5, ep = 0.0, u = -1.0, tev = 2.5e-2, az = 11.9,
    tevz = 2.53e-2, cliq = 0.0, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2, sigVal;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab (alpha.size()*beta.size());

  GIVEN( "invalid input requested (iinc != 1,2)" ){
    iinc = 0;
    THEN( "throw exception" ){
      REQUIRE_THROWS( sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, 
                           lat, cliq, sb, sb2, teff, iinc ) );
    } // THEN
  } // GIVEN


  GIVEN( "valid input values" ){
    WHEN( "many different E' are considered" ){
      std::vector<double> ePrimeVals {0.0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,5e-1},
                          correct_XS (ePrimeVals.size());

      AND_WHEN( "scattering cosine of -1.0" ){
        u = -1.0;

        correct_XS = { 0.0, 1292.387848, 2690.599752, 4057.900264, 4351.096726, 
                            1386.996889, 5.696163E-3, 0.0};
        THEN( "all output cross sections are correct" ) {
          for (size_t i = 0; i < ePrimeVals.size(); ++i ){
            sigVal = sig( e, ePrimeVals[i], u, tev, alpha, beta, sab, az, tevz, 
                           lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( correct_XS[i] == Approx( sigVal ).epsilon(1e-6) );
          } 
        } // THEN
      } // AND WHEN

      AND_WHEN( "scattering cosine of -0.5" ){
        u = -0.5;
        correct_XS = { 0.0, 1429.23802569, 3106.86309110, 4482.6813516, 
                            4503.96605993, 1361.12947151, 5.1043503E-3, 0.0};
        THEN( "all output cross sections are correct" ) {
          for (size_t i = 0; i < ePrimeVals.size(); ++i ){
            sigVal = sig( e, ePrimeVals[i], u, tev, alpha, beta, sab, az, tevz, 
                           lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( correct_XS[i] == Approx( sigVal ).epsilon(1e-6) );
          } 
        } // THEN
      } // AND WHEN

      AND_WHEN( "scattering cosine of +0.1" ){
        u = 0.1;

        correct_XS =  { 0.0, 1670.038013, 4010.983453, 5226.3490715, 
                             4709.2476163, 1327.738591, 4.4614069E-3, 0.0};
        THEN( "all output cross sections are correct" ) {
          for (size_t i = 0; i < ePrimeVals.size(); ++i ){
            sigVal = sig( e, ePrimeVals[i], u, tev, alpha, beta, sab, az, tevz, 
                           lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( correct_XS[i] == Approx( sigVal ).epsilon(1e-6) );
          } 
        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN
} // TEST CASE





TEST_CASE( "sig - bound scattering (iinc != 1)" ){
  int lasym = 0, lat, iinc = 2;

  double e = 1.0e-5, ep = 0.0, u = -0.2, tev = 2.5e-2, az = 11.9,
    tevz = 2.53e-2, cliq = 0.0, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2;

  double sigVal;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab (alpha.size()*beta.size());

  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 


  GIVEN ("alpha, beta values not scaled by 0.0253 (lat = 0)"){
    std::vector<double> uVals(4), eVals(4), epVals(4), xsVals(4), temps(4);
    lat = 0;

    WHEN( "either alpha or beta is outside of given range" ){
      tev = 2.53e-2; e = 1e-5, ep = 1e-1; u = -0.2;
      AND_WHEN( "various scattering cosines" ){
        uVals  = {-1.0,-0.5,0.1,0.9};
        xsVals = {1.233432018,1.1830012526,1.12387306,1.047426894 };
        THEN( "SCT Approximation is used" ){
          for (size_t i = 0; i < uVals.size(); ++i){
            sigVal = sig( e, ep, uVals[i], tev, alpha, beta, sab, az, tevz,
                          lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "requested alpha, beta values are within given bounds" ){
      e = 1.0e-2, ep = 1.2e-2; u = -0.2;
      AND_WHEN( "various temperatures" ){
        temps  = {1e-2,2.53e-2,1e-1,5e-1};
        xsVals = {279.23697546,115.536361596,29.457086622,5.7889762195};
        THEN( "interpolation is used" ){
          for (size_t t = 0; t < temps.size(); ++t){
            sigVal = sig( e, ep, u, temps[t], alpha, beta, sab, az, 
                          tevz, lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[t] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN 





  GIVEN ("alpha, beta values scaled by 0.0253 (lat = 1)"){
    std::vector<double> uVals(4), eVals(4), epVals(4), xsVals(4), temps(4);
    lat = 1;

    WHEN( "either alpha or beta is outside of given range" ){
      tev = 2.53e-2; e = 1e-5, ep = 1e-1; u = -0.2;
      AND_WHEN( "various scattering cosines" ){
        uVals  = {-1.0,-0.5,0.1,0.9};
        xsVals = {1.233432018,1.1830012526,1.12387306,1.047426894 };
        THEN( "SCT Approximation is used" ){
          for (size_t i = 0; i < uVals.size(); ++i){
            sigVal = sig( e, ep, uVals[i], tev, alpha, beta, sab, az, tevz,
                          lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN


    WHEN( "requested alpha, beta values are within given bounds" ){
      tev = 1.5e-1; e = 1.0e-2, ep = 1.2e-2; u = -0.2;

      AND_WHEN( "various scattering cosines" ){
        uVals  = {-1.0,-0.5,0.1,0.9};
        xsVals = {20.152903,20.144560,20.1297739287,20.067111072};
        THEN( "interpolation is used" ){
          for (size_t i = 0; i < uVals.size(); ++i){
            sigVal = sig( e, ep, uVals[i], tev, alpha, beta, sab, az, tevz, 
                          lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

      AND_WHEN( "various initial energies" ){
        eVals  = {1e-4,1e-3,1e-2,1e-1};
        xsVals = {194.41246264,61.754536926,20.13809684854,9.2971130024};
        THEN( "interpolation is used" ){
          for (size_t i = 0; i < eVals.size(); ++i){
            sigVal = sig( eVals[i], ep, u, tev, alpha, beta, sab, az, tevz, 
                          lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

      AND_WHEN( "various final energies" ){
        epVals = {1e-4,1e-3,1e-2,1e-1};
        xsVals = {1.9139438455,6.0459491756,0.0,46.71527444};
        THEN( "interpolation is used" ){
          for (size_t i = 0; i < epVals.size(); ++i){
            sigVal = sig( e, epVals[i], u, tev, alpha, beta, sab, az, 
                          tevz, lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

      AND_WHEN( "various temperatures" ){
        temps  = {1e-2,2.53e-2,1e-1,5e-1};
        xsVals = {275.153823,115.536361,30.106623,6.0696885};
        THEN( "interpolation is used" ){
          for (size_t t = 0; t < temps.size(); ++t){
            sigVal = sig( e, ep, u, temps[t], alpha, beta, sab, az, 
                          tevz, lasym, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( xsVals[t] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

    } // WHEN


    WHEN( "weird approx #1 (lasym = 0, cliq != 0, smallish alpha,beta" ){ 
      THEN( "weird approx invoked" ){
        lasym = 0; cliq  = 0.5; 
        
        e = 1.0e-4, ep = 1e-3, tev = 1.5e-1; u = -0.2;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      cliq, sb, sb2, teff, iinc );
        REQUIRE( 826.6527224 == Approx( sigVal ).epsilon(1e-6) );

        e = 1.1e-3, ep = 1e-3, tev = 3.e-2; u = -0.8;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      cliq, sb, sb2, teff, iinc );
        REQUIRE( 832.725737808 == Approx( sigVal ).epsilon(1e-6) );

      } // THEN
    } // WHEN


    WHEN( "weird approx #2 (largeish alpha,beta" ){ 

       for (double& alphaVal : alpha){ alphaVal *= 100; }
       for (double& betaVal  : beta ){ betaVal  *= 100; }


       sab[1] = -226.0;

        
        e = 1.0e-1, ep = 6e-1, tev = 2e-1; u = -0.2;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      cliq, sb, sb2, teff, iinc );

        REQUIRE( 2.78824367E-5 == Approx( sigVal ).epsilon(1e-6) );

      } // WHEN







  } // GIVEN
} // TEST CASE



