#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/e_ep_mu.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <range/v3/all.hpp>

/*

TEST_CASE( "313" ) {
  int lat = 1, jbeta = -7, iskip = 0;
  double enow = 1e-5, ep = 0.0, tev = 2.5507297688E-2;
  std::vector<double> betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0);


  GIVEN( "Various jbeta values, fixed enow" ){ 
    std::vector<int> jbetaVec(15), final_jbeta(15);
    std::vector<double> final_ep(15);

    jbetaVec    = ranges::view::iota(-7,8);
    final_jbeta = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7};
    final_ep    = {2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 
                   2.54e-3, 2.54e-3, 5.07e-3, 3.29e-2, 3.543e-2, 6.326e-2, 
                   6.579e-2, 9.362e-2};
    enow = 1e-5;
    for ( size_t i = 0; i < jbetaVec.size(); ++i ){
      jbeta = jbetaVec[i];
      ep = do_313( lat, jbeta, enow, betas, x, tev );
      REQUIRE( final_jbeta[i] == jbeta );
      REQUIRE( final_ep[i] == Approx(ep).epsilon(1e-6) );
    }


    final_jbeta = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7};
    final_ep    = {3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 
                   3.53e-3, 3.53e-3, 6.06e-3, 3.389e-2, 3.642e-2, 6.425e-2, 
                   6.678e-2, 9.461e-2};
    enow = 1e-3;
    for ( size_t i = 0; i < jbetaVec.size(); ++i ){
      jbeta = jbetaVec[i];
      ep = do_313( lat, jbeta, enow, betas, x, tev );
      REQUIRE( final_jbeta[i] == jbeta );
      REQUIRE( final_ep[i] == Approx(ep).epsilon(1e-6) );
    } 

    final_jbeta = {-7, -6, -5, -4, -3, -2, -1, 1, 1, 2, 3, 4, 5, 6, 7 };
    final_ep    = {6.39000153e-3, 3.4220001e-2, 3.6750001e-2, 6.4580001e-2, 
                   6.7110001e-2, 9.4940001e-2, 9.7470001e-2, 0.10253,  0.10253, 
                   0.10506, 0.13289, 0.13542, 0.16325, 0.16578, 0.19361 };
    enow = 1e-1;
    for ( size_t i = 0; i < jbetaVec.size(); ++i ){
      jbeta = jbetaVec[i];
      //std::cout << jbeta << std::endl;
      ep = do_313( lat, jbeta, enow, betas, x, tev );
      REQUIRE( final_jbeta[i] == jbeta );
      REQUIRE( final_ep[i] == Approx(ep).epsilon(1e-6) );
    } 

  } // GIVEN



} // TEST CASE
*/





TEST_CASE( "E-E'-mu" ){
  REQUIRE( true );

  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 
  double az = 0.99917, tol = 5e-2;
  int lasym = 0, lat = 1, iinc = 2, jmax = 55550;
  int nne = 88, nnl = -9, nl = 9;
  double T = 296.0, teff = 1397.671, teff2 = 0.0, sigma_b = 163.72792237360667, sigma_b2 = 0.0;

  std::cout.precision(15);

  e_ep_mu( T, teff, teff2, jmax, nne, nnl, nl, tol, sigma_b, sigma_b2, az, lasym, lat, iinc, alphas, betas, sab );


}
