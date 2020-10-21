#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep_util/adaptiveReconstruction.h"

TEST_CASE( "adaptive reconstruction" ){
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
      7.0, 7.65, 8.4, 9.15, 9.85, 10.0 };
  GIVEN( "inputs" ){

    std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
      betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };

    std::vector<double> sab(alphas.size()*betas.size());
    for ( size_t i = 0; i < alphas.size(); ++i ){
      for ( size_t j = 0; j < betas.size(); ++j ){
        sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 

    int iinc = 2, nbin = 4, lat = 1, lasym = 0;
    double temp = 296.0, emax = 0.625, tol = 0.5, az = 0.99917, sb = 163.72792237360667,
           sb2 = 0.0, teff = 0.12044192657731301;
    std::vector<double> boundXsVec {sb,sb2};

    adaptiveReconstruction( alphas, betas, sab, iinc, egrid, temp, nbin, emax, tol,
      lat, lasym, az, boundXsVec, teff );

    //auto out = adaptiveReconstruction( teff, cliq, iinc, tevz, lat, lasym, yy, 
    //  yu, sb, sb2, x, alpha, beta, sab, az, uj, sj, tol, tolmin, mumax, i, sum, 
    //  imax, enow, tev, j  );


 
    /*
    for ( size_t i = 0; i < yy.size(); ++i){ 
      REQUIRE( yyCorrect[i] == Approx(yy[i]).epsilon(1e-6) );
    }
    for ( size_t i = 0; i < yuCorrect_1_50.size(); ++i){ 
      REQUIRE( yuCorrect_1_50[i] == Approx(yu[i]).epsilon(1e-6) );
    }

    for ( size_t i = 0; i < x.size(); ++i){ 
      REQUIRE( xCorrect[i] == Approx(x[i]).epsilon(1e-6) );
    }

    for ( size_t i = 0; i < ujCorrect_1_50.size(); ++i){ 
      REQUIRE( ujCorrect_1_50[i] == Approx(uj[i]).epsilon(1e-6) );
    }
    for ( size_t i = ujCorrect_1_50.size(); i < uj.size(); ++i){ 
      REQUIRE( 0.0 == Approx(uj[i]).epsilon(1e-6) );
    }

    for ( size_t i = 0; i < sjCorrect_1_50.size(); ++i){ 
      REQUIRE( sjCorrect_1_50[i] == Approx(sj[i]).epsilon(1e-6) );
    }
    for ( size_t i = sjCorrect_1_50.size(); i < sj.size(); ++i){ 
      REQUIRE( 0.0 == Approx(sj[i]).epsilon(1e-6) );
    }
    */




    //double xl = std::get<0>(out);
    //double yl = std::get<1>(out);
    //REQUIRE( 1 == i );
    //REQUIRE( 50 == j );
    //REQUIRE( 0.9999999 == Approx(xl).epsilon(1e-6) );
    //REQUIRE( 31325.541195811933 == Approx(yl).epsilon(1e-6) );
    //REQUIRE( 66810.539071659470 == Approx(sum).epsilon(1e-6) );


  } // GIVEN

} // TEST CASE


/*

*/
