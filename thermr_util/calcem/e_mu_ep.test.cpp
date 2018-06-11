#define CATCH_CONFIG_MAIN 
#include "../../catch.hpp"
#include "e_mu_ep.h"


TEST_CASE( "Branch to handle E-mu-E' ordering" ){
  std::vector<double> alpha(40),beta(80);
  GIVEN( "inputs" ){
    int ncds = 0, iinc = 1, ltt = 0, lasym = 0, nnl = -9, nl = 9, jmax = 55550, 
      nne = 94, iprint=2, mtref = 229, mfh = 7, matdp = 1306, lat = 1, 
      inew = 10, imax = 20, iold = 11, math = 1065;
    double teff = 400.0, teff2 = 290.0, za = 6000.0, awr = 11.8969, emax = 1.2, 
      cliq = 0.0, t = 350.0, tol = 5e-2, az = 11.9, az2 = 0.0, sb = 12.5, 
      sb2 = 6.6, temp = 296.0, break_val = 3000, enow = 0.0, bk = 8.617385e-5, 
      tev = 0.0;
    std::vector<double> scr(500000,0.0), esi(95,0.0), xsi(95,0.0), x(20,0.0);

    for (size_t i = 0; i < alpha.size(); ++i ){ alpha[i] = 0.8*i*(i%2) + 0.01; }
    for (size_t i = 0; i < beta.size();  ++i ){ beta[i]  =     i*(i%5);        }

    std::vector<std::vector<double>> sab ( alpha.size(), std::vector<double> (beta.size()));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.5*i + 0.1*j;
      }
    }

    std::setprecision(15);

    e_mu_ep( ltt, math, matdp, mfh, mtref, tev, teff, teff2, bk, scr, za, awr, 
      ncds, emax, cliq, iinc, lat, esi, xsi, lasym, alpha, beta, sab, t, tol, 
      az, az2, sb, sb2, nnl, nl, jmax, nne, iprint, enow, break_val, temp, x, 
      inew, iold, imax );



  } // GIVEN

} // TEST CASE
