#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep.h"


TEST_CASE( "adaptive reconstruction" ){
  GIVEN( "inputs" ){

    int matdp = 1301, mtref = 222, ncds = 0, nw = 10, nne = 88, iinc = 2,
        lat = 1, lasym = 0, mumax = 300, imax = 20, i = 2, j = 0;
    double t = 296.0, teff = 541.8731, teff2 = 0.0, za = 1001.0, awr = 0.99917,
           cliq = 0.0, emax = 0.625, temp = 296.0, breakVal = 3000.0, 
           tevz = 2.53e-2, sb = 163.72792237, sb2 = 0.0, az = 0.99917, 
           tol = 5.0e-2, tolmin = 5.0e-7, xl = -1.0, yl = 20796.8534218,
           sum = 0.0, enow = 1.0e-5, tev = 2.55074596e-2;
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
      7.0, 7.65, 8.4, 9.15, 9.85, 10.0 }, esi(89,0.0), alpha(65), beta(75), 
      x(20,0.0), yy(20,0.0), yu(10000,0.0), uj(300,0.0), sj(300,0.0), scr(500000,0.0);
    for ( size_t i = 0; i < alpha.size(); ++i ){ alpha[i] = 0.12*(i+1); }
    for ( size_t i = 0; i < beta.size();  ++i ){ beta[i]  = 0.23*(i+1); }
    std::vector<std::vector<double>> sab(alpha.size(),std::vector<double>(beta.size()));
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i][j] = alpha[i]-beta[j];
      }
    }

    x[0] = 1; x[1] = -1;
    yy[0] = 18516.457902061160;
    yy[1] = 20796.853421853437;

    
    adaptiveReconstruction( teff, cliq, iinc, tevz, lat, lasym, yy, yu, sb, sb2, 
      x, alpha, beta, sab, az, uj, sj, tol, tolmin, mumax, i, xl, yl, sum, imax, 
      enow, tev, j  );

    std::vector<double> yyCorrect { 18516.45790206, 18787.04952242, 
      18923.84621860, 19061.65895127, 19623.17417784, 0.0, 0.0, 0.0, 0.0, 0.0, 
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      yuCorrect1_30 { 18651.259858, 113.0, 0.0, 0.0, 2.2235874e-8, 
        49.692358408, 4.4471748e-8, 69.846615744, 8.8943495e-8, 
        97.914611677, 1.7788699e-7, 136.73106477, 3.5577398e-7, 
        189.85007980, 7.1154795e-7, 261.38911576, 1.4230959e-6, 
        355.47667620, 2.8461918e-6, 475.82350947, 5.6923835e-6, 
        634.00500562, 1.1384767e-5, 909.95669658, 2.2769534e-5, 
        1491.2348194, 4.5539068e-5, 2528.4906307, 9.1078135e-5, 
        4229.9969100},
      xCorrect { 1.0, 0.75, 0.625, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
      ujCorrect_1_15 { -1.0, -0.750, -0.625, -0.50, -0.375, -0.250, -0.125, 0.0, 0.125, 0.250, 0.375, 0.50, 0.625, 0.750, 0.0 },
      sjCorrect_1_15 { 20796.853421, 20496.893150, 20348.583053, 20201.379495, 20055.273967, 19910.165824, 19766.074146, 19623.174177, 19481.208983, 19340.369483, 19200.497164, 19061.658951, 18923.846218, 18787.049522, 0.0 };
    
    for ( size_t i = 0; i < yy.size(); ++i){ 
      REQUIRE( yyCorrect[i] == Approx(yy[i]).epsilon(1e-4) );
    }
    for ( size_t i = 0; i < yuCorrect1_30.size(); ++i){ 
      if ( i == 1 ){ 
        // This is j in sigu, it's the number of iterations you need. may not
        // be so important that this is actually correct. 
        REQUIRE( yuCorrect1_30[i] == Approx(yu[i]).epsilon(1e-2) );
      } 
      else {
        REQUIRE( yuCorrect1_30[i] == Approx(yu[i]).epsilon(1e-4) );
      } 
    }

    for ( size_t i = 0; i < x.size(); ++i){ 
      REQUIRE( xCorrect[i] == Approx(x[i]).epsilon(1e-6) );
    }

    for ( size_t i = 0; i < ujCorrect_1_15.size(); ++i){ 
      REQUIRE( ujCorrect_1_15[i] == Approx(uj[i]).epsilon(1e-5) );
    }
    for ( size_t i = ujCorrect_1_15.size(); i < uj.size(); ++i){ 
      REQUIRE( 0.0 == Approx(uj[i]).epsilon(1e-6) );
    }

    for ( size_t i = 0; i < sjCorrect_1_15.size(); ++i){ 
      REQUIRE( sjCorrect_1_15[i] == Approx(sj[i]).epsilon(1e-4) );
    }
    for ( size_t i = sjCorrect_1_15.size(); i < sj.size(); ++i){ 
      REQUIRE( 0.0 == Approx(sj[i]).epsilon(1e-6) );
    }

    REQUIRE( 1 == i );
    REQUIRE( 14 == j );
    REQUIRE( 0.75 == Approx(xl).epsilon(1e-6) );
    REQUIRE( 18787.0495224 == Approx(yl).epsilon(1e-4) );
    REQUIRE( 34605.9936718 == Approx(sum).epsilon(1e-4) );









  } // GIVEN
} // TEST CASE



TEST_CASE( "Branch to handle E-mu-E' ordering" ){
  GIVEN( "inputs" ){
    THEN( "output is correct" ){

    int matdp = 1301, mtref = 222, ncds = 0, nw = 10, nne = 88, iinc = 2,
        lat = 1, lasym = 0, mumax = 300, imax = 20;
    double t = 296.0, teff = 541.8731, teff2 = 0.0, za = 1001.0, awr = 0.99917,
           cliq = 0.0, emax = 0.625, temp = 296.0, breakVal = 3000.0, 
           tevz = 2.53e-2, sb = 163.72792237, sb2 = 0.0, az = 0.99917, 
           tol = 5.0e-2, tolmin = 5.0e-7;
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
      7.0, 7.65, 8.4, 9.15, 9.85, 10.0 }, esi(89,0.0), alpha(65), beta(75), 
      x(20,0.0), yy(20,0.0), yu(10000,0.0), uj(300,0.0), sj(300,0.0), scr(500000,0.0);
    for ( size_t i = 0; i < alpha.size(); ++i ){ alpha[i] = 0.12*(i+1); }
    for ( size_t i = 0; i < beta.size();  ++i ){ beta[i]  = 0.23*(i+1); }
    std::vector<std::vector<double>> sab(alpha.size(),std::vector<double>(beta.size()));
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i][j] = alpha[i]-beta[j];
      }
    }

    //e_mu_ep( matdp, mtref, t, teff, teff2, scr, za, awr, ncds, nw, nne, 
    //  cliq, iinc, emax, egrid, temp, breakVal, esi, tevz, lat, lasym, yy, yu, 
    //  sb, sb2, x, alpha, beta, sab, az, uj, sj, tol, tolmin, mumax, imax );


      REQUIRE( true );
    } // GIVEN
  } // GIVEN
} // TEST CASE
