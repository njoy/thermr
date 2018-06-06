#include "../calcem_util/sigfig.h"


auto do313( int& jbeta, const int& lat, double& ep, double& enow, 
  const std::vector<double>& beta, const double& tev, const double& tevz,
  const std::vector<double>& x, int& iskip){
  bool do_313 = true;
  while ( do_313 ){
    //std::cout << 313 << std::endl;
    if (jbeta == 0) jbeta=1;
    if (jbeta <= 0) {
      if (lat == 1) {
        ep=enow-beta[-jbeta-1]*tevz;
      }
      else {
        ep=enow-beta[-jbeta-1]*tev;
      } // endif
      if (ep == enow) {
         ep=sigfig(enow,8,-1);
      }
      else {
         ep=sigfig(ep,8,0);
      }

    } 
    else {
      if (lat == 1) {
        ep=enow+beta[jbeta-1]*tevz;
      }
      else {
        ep=enow+beta[jbeta-1]*tev;
      } // endif
      if (ep == enow) {
         ep=sigfig(enow,8,+1);
         iskip=1;
      }
      else { 
         ep=sigfig(ep,8,0);
      }

    } // endif
    if (ep > x[2-1]) { return; } // go to 316
    jbeta=jbeta+1;
    // go to 313
  }

}


