#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "beginningLoop.h"

TEST_CASE( "110 120 130" ){
  GIVEN( "inputs" ){
    std::vector<double> x(65,0.0),y(65,0.0),s(65,0.0);
    x[0] = 1.0; x[1] = 0.99; x[2] = -1.0;
    y[0] = 1.35700e5; y[1] = 1.35701e5; y[2] = 1.35809e5;

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
    

    // Initialize S(a,b)
    std::vector<std::vector<double>> sab(alpha.size(), 
        std::vector<double>(beta.size(),0));
    for ( int i = 0; i < alpha.size(); ++i ){
      for ( int j = 0; j < beta.size(); ++j ){
        sab[i][j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 

    int imax = 20, lasym = 0, lat = 1, iinc = 2, nlmax = 65, nl = 10, i = 3, 
        j = 0, nbin = 8;

    double e = 1.0e-6, ep = 1.2e-4, tev = 1.5e-4, bbm = 0.0, az = 11.9,
      tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 1.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2, sum = 0.0, sigmin = 1.0e-32, 
      s1bb = 1.1369180380, tol = 2.5e-2, xtol = 1.0e-5, gral = 0.0, 
      seep = 91200.0, yl = 13500, ymax = 13500, fract = 0.0, xl = -1.0, 
      eps = -1.0e-3;



    do_110_120_130(i, imax, x, y, e, ep, tev, tevz, alpha,beta, sab, bbm, az, az2, lasym, teff, teff2, lat, cliq, sb, sb2, iinc, sum, nl, sigmin, s, nbin, fract, xl, j, ymax, eps, seep, gral,  yl, s1bb, tol, xtol);






  } // GIVEN
} // TEST CASE


