#include "calcem/calcem_util/sig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/begin_sigu.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/do150.h"
#include <range/v3/all.hpp>

#ifndef CALCEM_E_MU_EP_SIGU
#define CALCEM_E_MU_EP_SIGU


template <typename Range, typename Float>
auto do_150(int& i, int imax, Range& y, Range& x, Float tolmin, const Float& tol, 
  const Float& teff,Float& xm, const Float& e, 
  const Float& u, const Float& tev, const Range& alphas, const Range& betas, 
  const Range& sab, const Float& az, const Float& tevz, int lasym,int lat, 
  const Float& sb, const Float& sb2, int iinc ){

   while ( i != imax-1 ){
     std::cout << " --- 150 --- " << std::endl; 
     if ( i > 3-1 and 0.5*(y[i-1]+y[i-0])*(x[i-1]-x[i-0]) < tolmin ){ return; }
     xm = 0.5*(x[i-1]+x[i-0]); xm = sigfig(xm,8,0);
     if ( xm <= x[i-0] or xm >= x[i-1] ){ return; }

     Float yGuess = 0.5*(y[i-1]+y[i-0]),
           yTrue = sig( e, xm, u, tev, alphas, betas, sab, az, tevz, lasym, lat, sb, sb2, teff, iinc );

     // Point passes
     if ( abs(yTrue-yGuess) <= tol*abs(yTrue) ){ return; }

     // We need to bisect again
     i += 1;
     x[i] = x[i-1];
     y[i] = y[i-1];
     x[i-1] = xm;
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
   int j, jbeta, imax=20;
   Float sum, xl, yl, xm, test, tol, tolmin = 1.e-6, bmax = 20;
   Range x(imax), y(imax), s(2*nemax);
   Float tevz = 0.0253;

   tol=tolin;

   // adaptive calculation of cross section
   sum  = 0;
   x[0] = 0;
   y[0] = sig( e, x[0], u, tev, alphas, betas, sab, az, tevz, lasym, 
                 lat, sb, sb2, teff, iinc );

   jbeta = -betas.size();
   if (lasym > 0) jbeta=1;
   j = 0;
   xl = 0;
   yl = 0;


   // set up next panel
   while (true){  
     std::cout << " --- 111 --- " << std::endl;
     x[1] = x[0]; 
     y[1] = y[0];


     do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alphas, betas, sab, az, 
                 lasym, teff, sb, sb2, iinc );

     int i = 1;
     while ( true ){

       do_150( i, imax, y, x, tolmin, tol, teff, xm, e, u, tev, alphas, 
               betas, sab, az, tevz, lasym, lat, sb, sb2, iinc );

       while (true) {
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

         i -= 1;
         if ( i > 0 or jbeta+1 <= int(betas.size())){ break; }

         jbeta += 1;

         if ( i == 0 ){ continue; }
    
         std::cout << " --- 170 --- " << std::endl;
         s1[0] = sum; s2[0] = j;

         return;
      }

      if ( i > 0 ){ 
        continue; 
      }
      else if ( jbeta <= int(betas.size()) ){
        jbeta +=1;
        break;
      }
      else { 
          std::cout << "something is wrong" << std::endl;
          throw std::exception();
      }



    }

  }

   return;
}
#endif

