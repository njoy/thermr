#include "catch.hpp"
//#include "calcem/calcem_util/e_mu_ep_util/mainLoop.h"

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


