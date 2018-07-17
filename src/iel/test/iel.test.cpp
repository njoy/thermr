#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "iel/iel.h"


TEST_CASE( "iel" ){
  double x, y, xnext;
  int idis, ip, ir;
  std::vector<double> a(8);

  int mat = 11, itemp = 0, iold = 10, inew = 11, ne = 153, nex = 3, natom = 1, nbin = 8;
  double za = 1001.0, awr= 0.99917, emax = 0.625;
  std::vector<double> tempr = {296.0}, fl;
  //double[500000] scr;
  std::vector<double> scr(500000);
  for ( size_t i = 0; i < 100; ++i ){ scr[i] = i+1; }
  for ( size_t i = 100; i < scr.size(); ++i ){ scr[i] = 0.0; }
  std::vector<double> esi {1.0e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
    1.0e-4, 1.26e-4, 1.6e-4, 2.0e-4, 2.53e-4, 2.97e-4, 3.5e-4, 4.2e-4, 5.06e-4, 
    6.15e-4, 7.5e-4, 8.7e-4, 1.012e-3, 1.23e-3, 1.5e-3, 1.8e-3, 2.03e-3, 
    2.277e-3, 2.6e-3, 3.0e-3, 3.5e-3, 4.048e-3, 4.5e-3, 5.0e-3, 5.6e-3, 
    6.325e-3, 7.2e-3, 8.1e-3, 9.108e-3, 1.0e-2, 1.063e-2, 1.15e-2, 1.2397e-2, 
    1.33e-2, 1.417e-2, 1.5e-2, 1.6192e-2, 1.82e-2, 1.99e-2, 2.0493e-2, 2.15e-2, 
    2.28e-2, 2.53e-2, 2.8e-2, 3.0613e-2, 3.38e-2, 3.65e-2, 3.95e-2, 4.27570e-2, 
    4.65e-2, 5.0e-2, 5.6925e-2, 6.25e-2, 6.9e-2, 7.5e-2, 8.1972e-2, 9.0e-2, 9.6e-2, 
    0.1035, 0.111573, 0.12, 0.128, 0.1355, 0.145728, 0.16, 0.172, 0.184437, 0.2, 
    0.2277, 0.2510392, 0.2705304, 0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 
    0.4170351, 0.45, 0.5032575, 0.56, 0.625, 0.7, 0.0 };

  iel( mat, itemp, iold, inew, ne, nex, tempr, fl, za, scr, awr, emax, natom, 
      nbin, esi );


  GIVEN( ""  ){

  } // GIVEN

} // TEST CASE
