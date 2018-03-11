#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "terpa.h"


TEST_CASE( "terpa" ){
  double x, y, xnext;
  int idis, ip, ir;
  std::vector<double> a { 1.1, 2.2, 3.3, 5.5, 8.8, 13.13, 21.21, 34.34 };

  x = 1; y = 2; idis = 3; ip = 2; ir = 1;
  terpa( y, x, xnext, idis, a, ip, ir );
  std::cout << "---------   " <<  xnext << std::endl;
  

  GIVEN( ""  ){

  } // GIVEN

} // TEST CASE
