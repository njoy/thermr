#include "../sig.h"
auto do_113_116( int& jbeta, const int& lat, std::vector<double>& x, 
  std::vector<double>& y, const double& e, const double& tev, const double& tevz,
  const double& root1sq, const double& u,
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, double& bbm, const double& az, 
  const int lasym, const double& az2, const double& teff, const double& teff2, 
  const double& cliq, const double& sb, const double& sb2, const int& iinc){

  do {
    // 113 continue
    std::cout << 113 << std::endl;
    if (jbeta == 0) jbeta=1;

    if (lat == 1) { x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tevz; }
    else {          x[0] = e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tev;  }

    jbeta=jbeta+1;

  } while(x[0] <= x[1]);

  jbeta=jbeta-1;
   
  std::cout << 116 << std::endl;
  if (u < 0 and root1sq > 1.01*x[1] and root1sq < x[0]) {
    x[0]=root1sq;
  }
  y[0] = sig( e, x[0], u, tev, alpha, beta, sab, bbm, az, tevz, lasym, az2,
      teff2, lat, cliq, sb, sb2, teff, iinc );
} 



