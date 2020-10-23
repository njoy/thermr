#include "catch.hpp"
#include "inelastic/sig.h"
#include <range/v3/all.hpp>



TEST_CASE( "get alpha, beta indices" ){
  std::vector<double> alphas { 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0 },
                      betas  { 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0 };
  double alpha = 6e-2, beta = 0.8;
  REQUIRE( getIndices(alphas,alpha) == 4 );
  REQUIRE( getIndices(betas, beta ) == 6 );
  
  alphas = { 1.1, 2.2, 3.3, 4.5, 5.8 };
  betas  = { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };

  alpha = 1e-2, beta = 0.5;
  REQUIRE( getIndices(alphas,alpha) == 1 );
  REQUIRE( getIndices(betas, beta ) == 2 );

  alpha = 2.6, beta = 1.35;
  REQUIRE( getIndices(alphas,alpha) == 2 );
  REQUIRE( getIndices(betas, beta ) == 3 );

  alpha = 2.6, beta =-1.35;
  REQUIRE( getIndices(alphas,alpha) == 2 );
  REQUIRE( getIndices(betas, beta ) == 1 );

  alpha = 8.6, beta =-5.35;
  REQUIRE( getIndices(alphas,alpha) == 4 );
  REQUIRE( getIndices(betas, beta ) == 1 );

} // TEST CASE





TEST_CASE( "sig - free gas (iinc = 1)" ){
  int iinc = 1, lasym = 0, lat = 1;
  double e = 1.0e-5, ep = 0.0, u = -1.0, tev = 2.5e-2, az = 11.9,
    tevz = 2.53e-2, sb = 5.53, sb2 = 0.0, teff = 6.14e-2, sigVal;
  std::vector<double> boundXsVec {sb,sb2};

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab (alpha.size()*beta.size());

  GIVEN( "invalid input requested (iinc != 1,2)" ){
    iinc = 0;
    THEN( "throw exception" ){
      REQUIRE_THROWS( sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, 
                           lat, boundXsVec, teff, iinc ) );
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
                           lasym, lat, boundXsVec, teff, iinc );
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
                           lasym, lat, boundXsVec, teff, iinc );
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
                           lasym, lat, boundXsVec, teff, iinc );
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
    tevz = 2.53e-2, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2;

  double sigVal;
  std::vector<double> boundXsVec {sb,sb2};

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab (alpha.size()*beta.size());

  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 




  GIVEN ("alpha, beta values not scaled by 0.0253 (lat = 0)"){
      std::cout.precision(15);
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
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN
    WHEN( "requested alpha, beta values are within given bounds" ){
      e = 1.0e-2, ep = 1.2e-2; u = -0.2;
      AND_WHEN( "various temperatures" ){
        temps  = {1e-2,2.53e-2,1e-1,5e-1};
        xsVals = {279.23697546,420.408985357,214.631127663,96.373620249};

        temps  = {1e-2,2.53e-2,1e-1,5e-1};
        xsVals = {279.236975,115.5363615,29.4570866,5.7889762};

        THEN( "interpolation is used" ){
          for (size_t t = 0; t < temps.size(); ++t){
            sigVal = sig( e, ep, u, temps[t], alpha, beta, sab, az, 
                          tevz, lasym, lat, boundXsVec, teff, iinc );
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
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN

    WHEN( "requested alpha, beta values are within given bounds" ){
      tev = 1.5e-1; e = 1.0e-2, ep = 1.2e-2; u = -0.2;

      AND_WHEN( "various scattering cosines" ){
        uVals  = {-1.0,-0.5,0.1,0.9};
        xsVals = {20.1529032,20.14456083,20.129773,20.0671110};
        THEN( "interpolation is used" ){
          for (size_t i = 0; i < uVals.size(); ++i){
            sigVal = sig( e, ep, uVals[i], tev, alpha, beta, sab, az, tevz, 
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN


      AND_WHEN( "various initial energies" ){
        eVals  = {1e-5,1e-4,1e-3,1e-1};
        xsVals = {614.50162602,194.41246264,61.754536926,9.2971130024};
        THEN( "interpolation is used" ){
          for (size_t i = 0; i < eVals.size(); ++i){
            sigVal = sig( eVals[i], ep, u, tev, alpha, beta, sab, az, tevz, 
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

      AND_WHEN( "various final energies" ){
        epVals = {1e-5,1e-4,1e-3,1e-1};
        xsVals = {0.6052961363074,1.9139438455,6.0459491756,46.71527444};
        THEN( "interpolation is used" ){
          for (size_t i = 0; i < epVals.size(); ++i){
            sigVal = sig( e, epVals[i], u, tev, alpha, beta, sab, az, tevz,
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( xsVals[i] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

      AND_WHEN( "various temperatures" ){
        e = 1e-3;
        temps  = {1e-2,2.53e-2,1e-1,5e-1};
        xsVals = {554.3988417903,305.598988548,90.9490304062,19.00802666573};
        THEN( "interpolation is used" ){
          for (size_t t = 0; t < temps.size(); ++t){
            sigVal = sig( e, ep, u, temps[t], alpha, beta, sab, az, tevz, 
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( xsVals[t] == Approx( sigVal ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } // WHEN



    WHEN( "there might be a S(a,b) singularity at beta=0 (liquids)" ){ 
      THEN( "weird approx #1 invoked (do beta^2/alpha interpolation)" ){
        lasym = 0; 
        
        e = 1.0e-4, ep = 1e-3, tev = 1.5e-1; u = -0.2;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      boundXsVec, teff, iinc );
        REQUIRE( 57.424147 == Approx( sigVal ).epsilon(1e-6) );

        e = 1.1e-3, ep = 1e-3, tev = 3.e-2; u = -0.8;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      boundXsVec, teff, iinc );
        REQUIRE( 84.3981766 == Approx( sigVal ).epsilon(1e-6) );


        tev = 1.5e-1; e = 1.0e-2, ep = 1.2e-2; u = -0.2;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      boundXsVec, teff, iinc );
        REQUIRE( 20.1380968 == Approx( sigVal ).epsilon(1e-6) );


        tev = 1.5e-1; e = 1.0e-2, ep = 1.0e-2; u = -0.2;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      boundXsVec, teff, iinc );

        REQUIRE( 0.0 == Approx( sigVal ).epsilon(1e-6) );

        tev = 1e-2; e = 1e-2; ep = 1.2e-2; u = -0.2;
        sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz, lasym, lat, 
                      boundXsVec, teff, iinc );
        REQUIRE( 275.15382 == Approx( sigVal ).epsilon(1e-6) );




      } // THEN
    } // WHEN


    WHEN( "weird approx #2 (largeish alpha,beta" ){ 
      for (double& alphaVal : alpha){ alphaVal *= 100; }
      for (double& betaVal  : beta ){ betaVal  *= 100; }
      sab[1] = -226.0;
      e = 1.0e-1, ep = 6e-1, tev = 2e-1; u = -0.2;
      sigVal = sig(e,ep,u,tev,alpha,beta,sab,az,tevz,lasym,lat,boundXsVec,teff,iinc);
      REQUIRE( 2.78824367E-5 == Approx( sigVal ).epsilon(1e-6) );
    } // WHEN




    WHEN( "beta = 0" ){ 
      tev = 2.5e-2; e = 1e-1, ep = 1e-1; 
      teff = 0.12; iinc = 2; lasym = 0; az = 0.99917;
      lat = 0;

      std::vector<double> eVec {1e-1, 1e-2, 1e-3, 1e-5},
                         muVec {-1.0,-0.2,0.0,1.0},
                 correctVals = {
        // mu = -1.0      -0.2          0.0          1.0 
        0.27948044746, 0.503687172, 0.599753494, 0.0, // e,ep = 1e-1
        0.00000000000, 0.0        , 0.0        , 0.0, // e,ep = 1e-2
        0.0          , 0.0        , 0.0        , 0.0, // e,ep = 1e-3
        0.0          , 0.0        , 0.0        , 0.0  // e,ep = 1e-5
      };


      sb = 1.0;
      for (size_t i = 0; i < eVec.size(); ++i){
        for (size_t j = 0; j < muVec.size(); ++j){
          boundXsVec = {sb,sb2};
          REQUIRE( correctVals[i*muVec.size()+j] == 
                   Approx(sig(eVec[i],eVec[i],muVec[j],tev,alpha,beta,sab,az,
                          tevz,lasym,lat,boundXsVec,teff,iinc ) ).epsilon(1e-6) );
          boundXsVec = {4*sb,sb2};
          REQUIRE( 4.0*correctVals[i*muVec.size()+j] == 
                   Approx(sig(eVec[i],eVec[i],muVec[j],tev,alpha,beta,sab,az,
                      tevz,lasym,lat,boundXsVec,teff,iinc ) ).epsilon(1e-6) );
        }
      }
    } // WHEN
  } // GIVEN
} // TEST CASE












TEST_CASE( "sig - additional because i'm spooked" ){
  int lasym = 0, lat = 0, iinc = 2;

  double e = 1.0e-2, ep = 1.0e-2, u = 1.0, tev = 2.5507297687999999E-002, az = 0.99917,
    tevz = 2.53e-2, sb = 163.72792237360667, sb2 = 0.0,
    teff = 0.12044192657731301;

  double sigVal;
  std::vector<double> boundXsVec{sb,sb2};

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab (alpha.size()*beta.size());

  sab = { -0.18259619, -0.30201347, -3.93654779, -3.98809174, -4.33545607, 
          -4.39515402, -5.88934921, -0.76225291, -0.81658341, -3.14161459, 
          -3.30566186, -3.90554652, -3.96233362, -5.23696660, -1.19182884, 
          -1.23155471, -2.79610565, -2.95633099, -3.74989225, -3.80837585, 
          -4.93373911, -1.58342860, -1.61310713, -2.71233943, -2.84291608, 
          -3.69699590, -3.75199349, -4.77433858, -1.96121202, -1.98720663, 
          -2.78454600, -2.88531460, -3.71288120, -3.77142141, -4.71158392 };



  GIVEN ("alpha, beta values not scaled by 0.0253 (lat = 0)"){

    WHEN( "either alpha or beta is outside of given range" ){
      AND_WHEN( "various scattering cosines" ){
        THEN( "SCT Approximation is used" ){
            u = 1.0;
            sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz,
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( 2804297.7748918594 == Approx( sigVal ).epsilon(1e-6) );
            u = 0.95;
            sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz,
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( 14157.155384780059 == Approx( sigVal ).epsilon(1e-6) );
            u = -0.95;
            sigVal = sig( e, ep, u, tev, alpha, beta, sab, az, tevz,
                          lasym, lat, boundXsVec, teff, iinc );
            REQUIRE( 2093.6429488489475 == Approx( sigVal ).epsilon(1e-6) );



        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN 
} // TEST CASE







