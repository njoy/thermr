#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "sigl.h"


TEST_CASE( "sigl" ){
  GIVEN( "inputs" ){
    int lasym = 0, lat = 1, iinc = 2, nlmax = 65, 
        nlin = 10;

    double e = 1.0e-6, ep = 1.2e-4, tev = 1.5e-4, az = 11.9,
      tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2;

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
      s(65,0.0);

    std::vector<std::vector<double>> sab(alpha.size(), std::vector<double>(beta.size(),0));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 


    sigl( nlin, nlmax, e, ep, tev, alpha, beta, sab, s, tolin, az, tevz, iinc, 
        lat, lasym, az2, teff2, cliq, sb, sb2, teff );

    std::vector<double> correctS { 271591.653204, -8.776270547E-5 , 
      -6.162495039E-3, 7.1582888474E-5, -2.018780205E-2,  1.0016609703E-4, 
      -3.944407885E-2, 3.0130372017E-5, -5.459079490E-2, -9.1381657730E-5 };
   for ( int i = 0; i < 10; ++i ){ 
     REQUIRE( correctS[i] == Approx(s[i]).epsilon(1e-6) ); 
   }
   for ( int i = 10; i < s.size(); ++i ){
     REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
   }


    REQUIRE( true );

  } // GIVEN



  GIVEN( "inputs 2" ){
    int lasym = 0, lat = 1, iinc = 2, nlmax = 65, 
        nlin = 10;

    double e = 1.0e-6, ep = 1.2e-4, tev = 1.5e-4, az = 11.9,
      tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2;

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
      s(65,0.0);

    std::vector<std::vector<double>> sab(alpha.size(), std::vector<double>(beta.size(),0));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 


    sigl( nlin, nlmax, e, ep, tev, alpha, beta, sab, s, tolin, az, tevz, iinc, 
        lat, lasym, az2, teff2, cliq, sb, sb2, teff );

    std::vector<double> correctS { 271591.653204, -8.776270547E-5 , 
      -6.162495039E-3, 7.1582888474E-5, -2.018780205E-2,  1.0016609703E-4, 
      -3.944407885E-2, 3.0130372017E-5, -5.459079490E-2, -9.1381657730E-5 };
   for ( int i = 0; i < 10; ++i ){ 
     REQUIRE( correctS[i] == Approx(s[i]).epsilon(1e-6) ); 
   }
   for ( int i = 10; i < s.size(); ++i ){
     REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
   }


    REQUIRE( true );

  } // GIVEN



} // TEST CASE

