#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep_util/mainLoop.h"



TEST_CASE( " mu-E' ordering " ){
  int imax = 20, lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297688e-2, tol = 0.05;
  std::vector<double> 
  alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
  betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
  sab {-0.18259619, -0.30201347, -3.93654779, -3.98809174, -4.33545607, 
  -4.39515402, -5.88934921, -0.76225291, -0.81658341, -3.14161459, -3.30566186, 
  -3.90554652, -3.96233362, -5.23696660, -1.19182884, -1.23155471, -2.79610565, 
  -2.95633099, -3.74989225, -3.80837585, -4.93373911, -1.58342860, -1.61310713, 
  -2.71233943, -2.84291608, -3.69699590, -3.75199349, -4.77433858, -1.96121202, 
  -1.98720663, -2.78454600, -2.88531460, -3.71288120, -3.77142141, -4.71158392 };

  double az = 0.99917, sigma_b = 163.72792237, sigma_b2 = 0.0, teff = 0.120441926;
  std::vector<double> eVec, correctEnergies;
  double enow = 1e-5;
  do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);


} // TEST CASE


/*
TEST_CASE( "main loop" ){
  int nep, npage, j, k, ib, nw, ncds, nb;
  double sum;
  std::vector<double> scr_1_60 (60);
  GIVEN( "a nonzero k input value" ){
    THEN( "moreio is not invoked" ){
      {
        std::vector<double> yu(10000,0.0), scr(500000,0.0);
        nep = 64, npage = 306, j = 64, k = 8, ib = 75, nw = 8, ncds = 6, nb = 0;
        sum = 54321.0;
        yu[0] = 35587.0; yu[1] = 65.0;
        for ( size_t i = 2; i < 50; ++i ){ yu[i] = i; }
        mainLoop( nep, npage, j, k, ib, scr, yu, sum, nw, ncds, nb );
  
        REQUIRE( j    == 135  );
        REQUIRE( k    == 0    );
        REQUIRE( ib   == 64   );
        // This isn't working because of tab1io and moreio
        //REQUIRE( nw == 6    );
        REQUIRE( ncds == 1470 );
        REQUIRE( nb   == 0    );

        scr_1_60 = { 0, 0, 0, 0, 0, 0, 0, 0, 2, 1.1045452e-4, 4, 1.8409086e-4, 6, 
        2.5772721e-4, 8, 3.3136356e-4, 10, 4.0499990e-4, 12, 4.7863625e-4, 14, 
        5.5227260e-4, 16, 6.2590894e-4, 18, 6.9954529e-4, 20, 7.7318164e-4, 22, 
        8.4681798e-4, 24, 9.2045433e-4, 26, 9.9409068e-4, 28, 1.0677270e-3, 30, 
        1.1413633e-3, 32, 1.2149997e-3, 34, 1.2886360e-3, 36, 1.3622724e-3, 38, 
        1.4359087e-3, 40, 1.5095451e-3, 42, 1.5831814e-3, 44, 1.6568178e-3, 46, 
        1.7304541e-3, 48, 1.8040904e-3, 0, 0, 0, 0 };
        for ( size_t i = 0; i < scr_1_60.size(); ++i ){ 
          REQUIRE( scr_1_60[i] == Approx(scr[i]).epsilon(1e-6) ); 
        }
      }
      {
        nep = 50, npage = 500, j = 30, k = 12, ib = 40, nw = 8, ncds = 10, nb = 0;
        sum = 54321.0;
        std::vector<double> yu(10000,0.0), scr(500000,0.0);
        yu[0] = 32123.0; yu[1] = 100.0;
        for ( size_t i = 2; i < 50; ++i ){ yu[i] = i; }
        mainLoop( nep, npage, j, k, ib, scr, yu, sum, nw, ncds, nb );

        REQUIRE( j    == 111 );
        REQUIRE( k    == 0   );
        REQUIRE( ib   == 50  );
        // This isn't working because of tab1io and moreio
        //REQUIRE( nw == 112 );
        REQUIRE( ncds == 955 );
        REQUIRE( nb   == 0   );

        scr_1_60 = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1.1045452e-4, 4, 
        1.8409086e-4, 6, 2.5772721e-4, 8, 3.3136356e-4, 10, 4.0499990e-4, 12, 
        4.7863625e-4, 14, 5.5227260e-4, 16, 6.2590894e-4, 18, 6.9954529e-4, 20, 
        7.7318164e-4, 22, 8.4681798e-4, 24, 9.2045433e-4, 26, 9.9409068e-4, 28, 
        1.0677270e-3, 30, 1.1413633e-3, 32, 1.2149997e-3, 34, 1.2886360e-3, 36, 
        1.3622724e-3, 38, 1.4359087e-3, 40, 1.5095451e-3, 42, 1.5831814e-3, 44, 
        1.6568178e-3, 46, 1.7304541e-3, 48, 1.8040904e-3 };
        for ( size_t i = 0; i < scr_1_60.size(); ++i ){ 
          REQUIRE( scr_1_60[i] == Approx(scr[i]).epsilon(1e-6) ); 
        }
      } 
    } // THEN
  } // GIVEN
} // TEST CASE


*/


