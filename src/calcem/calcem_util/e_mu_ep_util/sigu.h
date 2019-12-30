#include "calcem/calcem_util/sig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/begin_sigu.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/do150.h"
#include <range/v3/all.hpp>

#ifndef CALCEM_E_MU_EP_SIGU
#define CALCEM_E_MU_EP_SIGU


template <typename Range, typename Float>
auto do_150(int& i, Range& y, Range& x, Float tolmin, const Float& tol, 
  const Float& teff, const Float& e, 
  const Float& u, const Float& tev, const Range& alphas, const Range& betas, 
  const Range& sab, const Float& az, const Float& tevz, int lasym,int lat, 
  const Float& sb, const Float& sb2, int iinc ){

   while ( i < int(x.size())-1 ){
     std::cout << " --- 150 --- " << std::endl; 
     if ( i > 2 and 0.5*(y[i-1]+y[i])*(x[i-1]-x[i]) < tolmin ){ return; }

     Float xMid = sigfig(0.5*(x[i-1]+x[i]),8,0);

     if ( xMid <= x[i] or xMid >= x[i-1] ){ return; }

     Float yGuess = 0.5*(y[i-1]+y[i-0]),
           yTrue = sig( e, xMid, u, tev, alphas, betas, sab, az, tevz, lasym, lat, sb, sb2, teff, iinc );

     // Point passes
     if ( abs(yTrue-yGuess) <= tol*abs(yTrue) ){ return; }

     // We need to bisect again
     ++i;
     x[i] = x[i-1];
     y[i] = y[i-1];
     x[i-1] = xMid;
     y[i-1] = yTrue;
   }
}




template <typename Range, typename Float>
inline auto sigu( int nemax, const Float& e, const Float& u, const Float& tev, 
  const Range& alphas, const Range& betas, const Range& sab, 
  const Float& tolin, const Float& az, const int& iinc, const int& lat, 
  const int& lasym, const Float& sb, const Float& sb2, const Float& teff, Range& s1, Range& s2 ){

  /*-------------------------------------------------------------------
   * Compute the secondary energy distribution scattering for cosine u.
   * Uses linear reconstruction with the cross section from function sig.
   *-------------------------------------------------------------------
   */
   using std::abs;
   int j=0, jbeta;
   Float sum = 0.0, xl=0.0, yl=0.0, tolmin = 1.e-6, bmax = 20;
   Range x(20,0.0), y(20,0.0), s(2*nemax,0.0);
   Float tevz = 0.0253;

   y[0] = sig( e, x[0], u, tev, alphas, betas, sab, az, tevz, lasym, 
                 lat, sb, sb2, teff, iinc );

   jbeta = (lasym > 0) ? 1 : -betas.size();


  do {
    std::cout << " --- 111 --- " << std::endl;
    x[1] = x[0]; 
    y[1] = y[0];


    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alphas, betas, sab, az, 
                lasym, teff, sb, sb2, iinc );

    int i = 1;
    do {

      do_150( i, y, x, tolmin, tolin, teff, e, u, tev, alphas, 
              betas, sab, az, tevz, lasym, lat, sb, sb2, iinc );

      do {
        std::cout << " --- 160 --- " << std::endl;
        j += 1;
        s1[j] = x[i];
        s2[j] = y[i];
        if ( j > 1 ) { sum += (y[i]+yl)*(x[i]-xl); }
        xl = x[i];
        yl = y[i];

        if (( j >= nemax-1 ) or ( jbeta > 0 and betas[jbeta-1] > bmax )){ 
          std::cout << " --- 170 --- " << std::endl;
          s1[0] = sum; s2[0] = j;
          return;
        } 

        --i;
        if ( i > 0 ){ break; }
        ++jbeta;
        if ( jbeta <= int(betas.size())) { break; }
    
      } while (i == 0);

    } while(i > 0);

  } while( jbeta <= int(betas.size()));
  std::cout << " --- 170 --- " << std::endl;
  s1[0] = sum; s2[0] = j;
}
#endif

