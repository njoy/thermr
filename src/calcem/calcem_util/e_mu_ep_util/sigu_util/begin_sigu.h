#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"

#ifndef CALCEM_E_MU_EP_SIGU_DO_113_116
#define CALCEM_E_MU_EP_SIGU_DO_113_116

template <typename Range, typename Float>
inline auto do_113_116( int& jbeta, const int& lat, Range& x, 
  Range& y, const Float& e, Float& tev, const Float& tevz,
  const Float& u,
  const Range& alpha, const Range& beta, 
  const Range& sab, Float& az, 
  const int lasym, const Float& teff, 
  const Float& sb, const Float& sb2, const int& iinc){

  Float root1 = (u*sqrt(e)+sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);
  //std::cout << root1 << std::endl;
  //root1 = 2.6720021028992805E-002;
  //std::cout << 113 << std::endl;
  do {
    // 113 continue
    if (jbeta == 0) jbeta=1;
    if (lat == 1) { x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tevz; }
    else {          x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tev;  }
    x[0] = sigfig(x[0],8,0);
    if ( jbeta <= 0 and x[0] == e ){ x[0] = sigfig(e,8,-1); }
    ++jbeta;

  } while(x[0] <= x[1]);

  jbeta = jbeta-1;
   
  //std::cout << 116 << std::endl;
  if (u < 0 and root1*root1 > 1.01*x[1] and root1*root1 < x[0]) {
    x[0]=root1*root1;
    //std::cout << x[0] << std::endl;
  }
  //std::cout << tev << "   " << az << "   " << sb << "    " << sb2 << "   " << teff << std::endl;

  //x[0] = sigfig(x[0],8,0);
  x[0] = sigfig(x[0],8,1);
  //std::cout << x[0] << std::endl;
  //std::cout << std::setprecision(25) << x[0] << "    " << x[1] << "     " << x[2] << std::endl;
  //std::cout << sb2 << "        " << iinc << std::endl;
  //std::cout << std::setprecision(25) << "x[0]    " << x[0] << std::endl;
  //std::cout << x[0] << std::endl;
  //x[0] = 7.1395952000007134E-004;
  y[0] = sig( e, x[0], u, tev, alpha, beta, sab, az, tevz, lasym, /*az2,
      teff2,*/ lat, sb, sb2, teff, iinc );
  //    std::cout << x[0] << "      " << y[0] << std::endl;
} 

#endif


