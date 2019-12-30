#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"

#ifndef CALCEM_E_MU_EP_SIGU_DO_113_116
#define CALCEM_E_MU_EP_SIGU_DO_113_116

template <typename Range, typename Float>
inline auto do_113_116( int& jbeta, const int& lat, Range& x, Range& y, 
  const Float& e, const Float& tev, const Float& tevz, const Float& u, 
  const Range& alpha, const Range& beta,  const Range& sab, const Float& az, 
  const int lasym, const Float& teff, const Float& sb, const Float& sb2, 
  const int& iinc){

  Float root1 = (u*sqrt(e)+sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);
  do {
    // 113 continue
    std::cout << " --- 113 --- " << std::endl;
    if (jbeta == 0) jbeta=1;
    if (lat == 1) { x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tevz; }
    else {          x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tev;  }
    x[0] = sigfig(x[0],8,0);
    if ( jbeta <= 0 and x[0] == e ){ x[0] = sigfig(e,8,-1); }
    ++jbeta;

  } while(x[0] <= x[1]);

  jbeta = jbeta-1;
   
  std::cout << " --- 116 --- " << std::endl;
  if (u < 0 and root1*root1 > 1.01*x[1] and root1*root1 < x[0]) {
    x[0] = root1*root1;
  }

  x[0] = sigfig(x[0],8,0);
  y[0] = sig( e, x[0], u, tev, alpha, beta, sab, az, tevz, lasym, lat, sb, sb2, 
              teff, iinc );
} 

#endif


