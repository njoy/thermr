#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "upstk.h"

TEST_CASE( "uptsk" ){
  int nl, nx, i = 1;
  std::vector<double> s(6);
  double e;

  GIVEN( "legendre order nl of 1" ){
    nl = 1;
    e = 4.0e-4;

    WHEN( "cycle length is small (nx=2)" ){
      nx = 2;
      std::vector<std::vector<double>> stk ( nx, std::vector<double> (20,0.0) );
      s = {0,0,0,0,0,0};

      AND_WHEN( "stk is initialized to nearly all zeros" ){
        stk[0][0] = 4.5e-4;
        upstk( e, s, stk, nl, nx, i );
        std::vector<double> stkVals { 4.0e-4, 0, 4.5e-4, 0 };

        THEN( "output is correct" ){
          REQUIRE( stk[0][0] == Approx(4.0e-4).epsilon(1e-6) );
          REQUIRE( stk[1][0] == Approx(0.0).epsilon(1e-6) );
          REQUIRE( stk[0][1] == Approx(4.5e-4).epsilon(1e-6) );
          REQUIRE( stk[1][1] == Approx(0.0).epsilon(1e-6) );
          for ( int i = 2; i < 20; ++i ){
            REQUIRE( stk[0][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[1][i] == Approx(0.0).epsilon(1e-6) );
          }
        } // THEN
  
      } // AND WHEN
  
      AND_WHEN( "stk is initialized to some nonzero values" ){

        for ( int i = 0; i < 5; ++i ){
          stk[0][i] = 0.2*( i + 0.5 );
          stk[1][i] = 0.2*( i + 1.0 );
        }


        upstk( e, s, stk, nl, nx, i );
        // These are only the nonzero values in the desired output. In the 
        // for loop down there, we want all values beyond stkVals's range to be 0
        std::vector<double> stkVals {4e-4, 0, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1};

        THEN( "output is correct" ){
          int j = 0;
          for ( int i = 0; i < 20; ++i ){
            if (j < stkVals.size() ){
              REQUIRE( stk[0][i] == Approx(stkVals[j]).epsilon(1e-6) );
              REQUIRE( stk[1][i] == Approx(stkVals[j+1]).epsilon(1e-6) );
              j += 2;
            }
            else {
              REQUIRE( stk[0][i] == Approx(0.0).epsilon(1e-6) );
              REQUIRE( stk[1][i] == Approx(0.0).epsilon(1e-6) );
            }
          }
        } // THEN
      } // AND WHEN
    } // WHEN

    WHEN( "cycle length is moderate (nx=4)" ){

      std::vector<double> s(6,0.0);
      s[0] = 1.5e3;
      int nl = 1, nx = 4, i = 1;
      double e = 1.8e-3;
      std::vector<std::vector<double>> stk ( nx, std::vector<double> (20,0.0) );

      for ( int i = 0; i < 3; ++i ){
        stk[0][i] = i * 0.4 + 0.1;
        stk[1][i] = i * 0.4 + 0.2;
        stk[2][i] = i * 0.4 + 0.3;
        stk[3][i] = i * 0.4 + 0.4;
      }

      std::vector<double> stkVals { 1.8e-3, 1500.0, 0.3, 0.4, 0.1, 0.2, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2 };

        upstk( e, s, stk, nl, nx, i );
        // These are only the nonzero values in the desired output. In the 
        // for loop down there, we want all values beyond stkVals's range to be 0
        int j = 0;
        for ( int i = 0; i < 20; ++i ){
          if (j < stkVals.size() ){
            REQUIRE( stk[0][i] == Approx(stkVals[j]).epsilon(1e-6) );
            REQUIRE( stk[1][i] == Approx(stkVals[j+1]).epsilon(1e-6) );
            REQUIRE( stk[2][i] == Approx(stkVals[j+2]).epsilon(1e-6) );
            REQUIRE( stk[3][i] == Approx(stkVals[j+3]).epsilon(1e-6) );
            j += 4;
          }
          else {
            REQUIRE( stk[0][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[1][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[2][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[3][i] == Approx(0.0).epsilon(1e-6) );
          }
        }
   


    } // WHEN
  } // GIVEN

  GIVEN( "legendre order nl of 3" ){
    std::vector<double> s { 1.5e3, 2.5e3, 3.5e3, 4.5e3, 5.5e3, 6.5e3 };

    WHEN( "cycle length is large (nx=6)" ){
      int nl = 3, nx = 6, i = 1;
      std::vector<std::vector<double>> stk ( nx, std::vector<double> (20,0.0) );
      double e = 1.8e-3;

      for ( int i = 0; i < 3; ++i ){
        stk[0][i] = i * 0.6 + 0.1;
        stk[1][i] = i * 0.6 + 0.2;
        stk[2][i] = i * 0.6 + 0.3;
        stk[3][i] = i * 0.6 + 0.4;
        stk[4][i] = i * 0.6 + 0.5;
        stk[5][i] = i * 0.6 + 0.6;
      }

      std::vector<double> stkVals { 1.8e-3, 1500, 2500, 3500, 0.5, 0.6, 0.1, 0.2, 0.3, 0.4, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8 };
      upstk( e, s, stk, nl, nx, i );
      
      THEN( "output is correct" ){
        int j = 0;
        for ( int i = 0; i < 20; ++i ){
          if (j < stkVals.size() ){
            REQUIRE( stk[0][i] == Approx(stkVals[j]).epsilon(1e-6) );
            REQUIRE( stk[1][i] == Approx(stkVals[j+1]).epsilon(1e-6) );
            REQUIRE( stk[2][i] == Approx(stkVals[j+2]).epsilon(1e-6) );
            REQUIRE( stk[3][i] == Approx(stkVals[j+3]).epsilon(1e-6) );
            REQUIRE( stk[4][i] == Approx(stkVals[j+4]).epsilon(1e-6) );
            REQUIRE( stk[5][i] == Approx(stkVals[j+5]).epsilon(1e-6) );
            j += 6;
          }
          else {
            REQUIRE( stk[0][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[1][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[2][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[3][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[4][i] == Approx(0.0).epsilon(1e-6) );
            REQUIRE( stk[5][i] == Approx(0.0).epsilon(1e-6) );
          }
        }
 

      } // THEN
    } // WHEN


    WHEN( "cycle length is moderate (nx=4)" ){
      int nl = 3, nx = 4, i = 1;
      std::vector<std::vector<double>> stk ( nx, std::vector<double> (20,0.0) );

      for ( int i = 0; i < 3; ++i ){
        stk[0][i] = i * 0.4 + 0.1;
        stk[1][i] = i * 0.4 + 0.2;
        stk[2][i] = i * 0.4 + 0.3;
        stk[3][i] = i * 0.4 + 0.4;
      }

      double e = 1.8e-3;

      upstk( e, s, stk, nl, nx, i );
      std::vector<double> stkVals { 1.8e-3, 1500, 2500, 3500, 0.1, 0.2, 0.3, 0.4, 0.9, 1.0, 1.1, 1.2 };
      int j = 0;
      for ( int i = 0; i < 20; ++i ){
        if (j < stkVals.size() ){
          REQUIRE( stk[0][i] == Approx(stkVals[j]).epsilon(1e-6) );
          REQUIRE( stk[1][i] == Approx(stkVals[j+1]).epsilon(1e-6) );
          REQUIRE( stk[2][i] == Approx(stkVals[j+2]).epsilon(1e-6) );
          REQUIRE( stk[3][i] == Approx(stkVals[j+3]).epsilon(1e-6) );
          j += 4;
        }
        else {
          REQUIRE( stk[0][i] == Approx(0.0).epsilon(1e-6) );
          REQUIRE( stk[1][i] == Approx(0.0).epsilon(1e-6) );
          REQUIRE( stk[2][i] == Approx(0.0).epsilon(1e-6) );
          REQUIRE( stk[3][i] == Approx(0.0).epsilon(1e-6) );
        }
      }
    } // WHEN
  } // GIVEN
} // TEST CASE



