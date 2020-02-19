#ifndef COH_SIGCOH_FORM
#define COH_SIGCOH_FORM
#include <cmath>
#include <vector>

double form( int lat, int l1, int l2, int l3 ){
  /*-------------------------------------------------------------------
   * Compute form factors for the specified lattice.
   *       lat = 1  graphite
   *       lat = 2  be
   *       lat = 3  beo
   *
   * This seems to be approximately the same as leapr's form function in
   * coher/coher_util/formf.h. The only/main difference is that this does
   * not give the option for fcc or bcc lattices
   *-------------------------------------------------------------------
   */
  using std::pow;
  if (lat == 1){      // graphite
    return ( l3 % 2 == 0 ) ? 0.25*(6+10*cos(2*M_PI*(l1-l2)/3))
                           :   std::pow(sin(1*M_PI*(l1-l2)/3),2); }
  else if (lat == 2){ // beryllium
    return 1.0 + cos(2*M_PI*(2*l1+4*l2+3*l3)/6);
  }
  else if (lat == 3){ // beryllium oxide
    double beo1 = 7.54, beo2 = 4.24, beo3 = 11.31;
    return (1+cos(2*M_PI*(2*l1+4*l2+3*l3)/6)) 
         * (beo1+beo2+beo3*cos(M_PI*(3*l3)*0.25));
  }
  std::cout << "invalid lat" << std::endl;
  throw std::exception();
}

#endif
