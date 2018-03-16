#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "contio.h"


TEST_CASE( "contio" ){
  int nin = 26, nout = 0, nscr = 0, nb = 0, nw = 17;
  contio( nin, nout, nscr, nb, nw );
  

  GIVEN( ""  ){

  } // GIVEN

} // TEST CASE
