#include "../sig.h"
auto do_113_116( int& jbeta, const int& lat, std::vector<double>& x, 
  std::vector<double>& y, const double& e, const double& tev, const double& tevz,
  const double& root1, const double& u, const double& xl, const double& yl, 
  const std::vector<double>& alpha, const std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, double& bbm, const double& az, 
  const int lasym, const double& az2, const double& teff, const double& teff2, 
  const double& cliq, const double& sb, const double& sb2, int& i, const int& iinc){
  while (true){
    // 113 continue
    std::cout << 113 << std::endl;
    if (jbeta == 0) jbeta=1;
      if (lat == 1) {
        x[1-1]=e + jbeta /abs(jbeta) * beta[abs(jbeta)-1]*tevz;
      }
      else { 
        x[1-1]=e + jbeta / abs(jbeta) * beta[abs(jbeta)-1]*tev;
      }
    // if (x[1-1] > x[2-1]) go to 116
    if (x[1-1] > x[2-1]) break;

    jbeta=jbeta+1;
    // go to 113
  } // This is doing 113
   
  // 116 continue
  std::cout << 116 << std::endl;
  if (u < 0 and root1*root1 > 1.01*x[2-1] and root1*root1 < x[1-1]) {
    x[1-1]=root1*root1;
  }
  y[1-1] = sig( e, x[0], xl, tev, alpha, beta, sab, bbm, az, tevz, lasym, az2,
      teff2, lat, cliq, sb, sb2, teff, iinc );
  i = 2;
} 



