#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/e_ep_mu.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"




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
  double az = 0.99917, tol = 2.5e-2;
  int lasym = 0, lat = 1, iinc = 2; 



  int nne = 88, nnl = -9, nl = 9;
  double T = 296.0, teff = 1397.671, teff2 = 0.0, sigma_b = 163.72792237360667, sigma_b2 = 0.0;
  e_ep_mu( T, teff, teff2, nne, nnl, nl, tol, sigma_b, sigma_b2, az, lasym, lat, iinc, alphas, betas, sab );


}

