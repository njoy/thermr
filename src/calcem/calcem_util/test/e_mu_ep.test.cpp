#include "catch.hpp"
//#include "calcem/calcem_util/e_mu_ep.h"

/*
TEST_CASE( "Branch to handle E-mu-E' ordering" ){
  GIVEN( "inputs" ){
    THEN( "output is correct" ){

    int matdp = 1301, mtref = 222, ncds = 0, nw = 10, nne = 88, iinc = 2,
        lat = 1, lasym = 0, mumax = 300, imax = 20, nep = 0, npage = 306,
        ib = 75, nb = 0;
    double t = 296.0, teff = 541.8731, teff2 = 0.0, za = 1001.0, awr = 0.99917,
           cliq = 0.0, emax = 0.625, temp = 296.0, breakVal = 3000.0, 
           tevz = 2.53e-2, sb = 163.72792237, sb2 = 0.0, az = 0.99917, 
           tol = 5.0e-2, tolmin = 5.0e-7, yumin = 2e-7;
    std::vector<double> egrid { 1.0e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
      1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 2.53e-4, 2.970e-4, 3.5e-4, 4.2e-4, 5.06e-4, 
      6.15e-4, 7.5e-4, 8.7e-4, 1.012e-3, 1.23e-3, 1.5e-3, 1.8e-3, 2.030e-3, 
      2.277e-3, 2.6e-3, 3.0e-3, 3.5e-3, 4.048e-3, 4.5e-3, 5.0e-3, 5.6e-3, 
      6.325e-3, 7.28e-3, 8.1e-3, 9.108e-3, 1.e-2, 1.063e-2, 1.15e-2, 1.2397e-2, 
      1.33e-2, 1.4170e-2, 1.5e-2, 1.6192e-2, 1.820e-2, 1.990e-2, 2.0493e-2, 
      2.15e-2, 2.280e-2, 2.53e-2, 2.8e-2, 3.0613e-2, 3.38e-2, 3.65e-2, 3.95e-2, 
      4.2757e-2, 4.65e-2, 5.3e-2, 5.6925e-2, 6.25e-2, 6.9e-2, 7.5e-2, 8.1972e-2, 
      9.0e-2, 9.6e-2, 0.1035, 0.1115730, 0.12, 0.128, 0.1355, 0.1457280, 0.16, 
      0.172, 0.184437, 0.20, 0.2277, 0.2510392, 0.2705304, 0.2907501, 0.3011332, 
      0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 0.5032575, 0.56, 0.625, 0.7, 
      0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 1.7, 1.855, 2.02, 2.18, 
      2.36, 2.59, 2.855, 3.120, 3.42, 3.75, 4.07, 4.46, 4.9, 5.35, 5.85, 6.4, 
      7.0, 7.65, 8.4, 9.15, 9.85, 10.0 }, esi(89,0.0), xsi(89,0.0), alpha(65), 
      beta(75), x(20,0.0), yy(20,0.0), yu(10000,0.0), uj(300,0.0), sj(300,0.0), 
      scr(500000,0.0), ubar(118,0.0);
    for ( size_t i = 0; i < alpha.size(); ++i ){ alpha[i] = 0.12*(i+1); }
    for ( size_t i = 0; i < beta.size();  ++i ){ beta[i]  = 0.23*(i+1); }
    std::vector<std::vector<double>> sab(alpha.size(),std::vector<double>(beta.size()));
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i][j] = alpha[i]-beta[j];
      }
    }

    e_mu_ep( matdp, mtref, t, teff, teff2, scr, za, awr, ncds, nw, nne, cliq, 
      iinc, emax, egrid, temp, breakVal, esi, tevz, lat, lasym, yy, yu, sb, sb2, 
      x, alpha, beta, sab, az, uj, sj, tol, tolmin, mumax, imax, ubar, xsi, nep, 
      yumin, npage, ib, nb );

    REQUIRE( 4.6695291e-2 == Approx(teff).epsilon(1e-6) );
    REQUIRE( 0.0 == Approx(teff2).epsilon(1e-6) );
    //std::cout << ncds << std::endl;
    //REQUIRE( 6149748.0 == Approx(ncds).epsilon(1e-6) );


    } // THEN
  } // GIVEN
} // TEST CASE
*/
