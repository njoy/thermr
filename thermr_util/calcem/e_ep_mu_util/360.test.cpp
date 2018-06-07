#define CATCH_CONFIG_MAIN 
#include "../../../catch.hpp"
#include "360.h"

TEST_CASE( "360" ){
  int j = 0, jmax = 55550, ie = 1, nll = 32673, nl = 9;
  std::vector<double> xsi(95,0.0), x(20,0.0);
  x[0] = 1e-5; x[1] = 5e-5; x[2] = 2.5e-5;
  double ulast = 1.1, u2last = 1.2, u3last = 1.3, xlast = 2.3, ylast = 3.4,
         tolmin = 5.0e-6; 
  std::vector<std::vector<double>> y(65,std::vector<double>(20));
  for ( size_t i = 0; i < y.size(); ++i ){
    for ( size_t j = 0; j < y[0].size(); ++j ){
      y[i][j] = 0.01*(i+j);
    }
  }

  GIVEN( " " ){
    WHEN( "  " ){
      THEN( "  " ){

      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
