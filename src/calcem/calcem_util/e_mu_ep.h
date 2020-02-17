#include "calcem/calcem_util/e_mu_ep_util/mainLoop.h"




template <typename Range, typename Float>
auto e_mu_ep( Range& eVec, const Float& tev, const Float& tol, 
  const int lat, const int iinc, const int lasym, const Range& alphas, 
  const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, 
  const Float& sigma_b2, const Float& teff ){

  for ( size_t energy_i = 0; energy_i < eVec.size(); ++energy_i ){
    mu_ep( eVec[energy_i], tev, tol, lat, iinc, lasym, alphas, betas, sab, az, 
           sigma_b, sigma_b2, teff );
  }
}

