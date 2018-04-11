#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "upstk.h"

TEST_CASE( "uptsk" ){
  GIVEN( "inputs" ){
    WHEN( "inputs" ){
      int nl = 1, nx = 2, i = 1;
      double e = 4.0e-4;
      std::vector<double> s(6,0.0);
      std::vector<std::vector<double>> stk ( 2, 
          std::vector<double> (20,0.0) );
      stk[0][0] = 4.5e-4;

      upstk( e, s, stk, nl, nx, i );
      REQUIRE( stk[0][0] == Approx(4.0e-4).epsilon(1e-6) );
      REQUIRE( stk[1][0] == Approx(0.0).epsilon(1e-6) );
      REQUIRE( stk[0][1] == Approx(4.5e-4).epsilon(1e-6) );
      REQUIRE( stk[1][1] == Approx(0.0).epsilon(1e-6) );
      for ( int i = 2; i < 20; ++i ){
        REQUIRE( stk[0][i] == Approx(0.0).epsilon(1e-6) );
        REQUIRE( stk[1][i] == Approx(0.0).epsilon(1e-6) );
      }
        

      THEN( "inputs" ){


      } // THEN
    } // WHEN
    WHEN( "inputs" ){
      double e = 4.0e-4;
      std::vector<double> s(6,0.0);
      int nl = 1, nx = 2, i = 1;
      std::vector<std::vector<double>> stk ( 2, 
          std::vector<double> (20,0.0) );
      stk[0][0] = 4.5e-4;

      double val = 0.1;
      for ( int i = 0; i < 5; ++i ){
        stk[0][i] = val;
        stk[1][i] = val + 0.1;
        val += 0.2;
      }


      upstk( e, s, stk, nl, nx, i );
      // These are only the nonzero values in the desired output. In the 
      // for loop down there, we want all values beyond stkVals's range to be 0
      std::vector<double> stkVals {4e-4, 0, 0.1, 0.2, 0.5, 0.6, 0.7, 0.8, 0.9, 1};
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
      /*
      std::cout << "0  0   " <<  stk[0][0] << std::endl; 
      std::cout << "1  0   " <<  stk[1][0] << std::endl; 
      std::cout << "0  1   " <<  stk[0][1] << std::endl; 
      std::cout << "1  1   " <<  stk[1][1] << std::endl; 
      std::cout << "0  2   " <<  stk[0][2] << std::endl; 
      std::cout << "1  2   " <<  stk[1][2] << std::endl; 
      */

      THEN( "inputs" ){


      } // THEN
    } // WHEN

  } // GIVEN
} // TEST CASE
