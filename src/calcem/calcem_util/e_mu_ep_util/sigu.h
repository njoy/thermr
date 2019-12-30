#include "calcem/calcem_util/sig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/begin_sigu.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/do150.h"
#include <range/v3/all.hpp>

#ifndef CALCEM_E_MU_EP_SIGU
#define CALCEM_E_MU_EP_SIGU


template <typename Range, typename Float>
auto do_150(int& i, int imax, Range& y, Range& x, Float tolmin, 
            const Float& tol,
            const Float& teff,
            Float& xm, Float& ym, Float& yt, const Float& e, const Float& u, 
            const Float& tev, const Range& alphas, const Range& betas, 
            const Range& sab, const Float& az, const Float& tevz, int lasym,
            int lat, const Float& sb, const Float& sb2, int iinc ){
   while (true){
     std::cout << " --- 150 --- " << std::endl; 
     if ( i == imax ){ break; } 
     if ( i > 3 and 0.5*(y[i-2]+y[i-1])*(x[i-2]-x[i-1]) < tolmin ){ break; }
     xm = 0.5*(x[i-2]+x[i-1]);
     xm = sigfig(xm,8,0);
     if ( xm <= x[i-1] or xm >= x[i-2] ){ break; }
     ym = 0.5*(y[i-2]+y[i-1]);
     yt = sig( e, xm, u, tev, alphas, betas, sab, az, tevz, lasym, 
               lat, sb, sb2, teff, iinc );

     Float test = tol*abs(yt);
     if ( abs(yt-ym) <= test ){ break; }

     i += 1;
     x[i-1] = x[i-2];
     y[i-1] = y[i-2];
     x[i-2] = xm;
     y[i-2] = yt;
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
   int i, j, jbeta, imax=20;
   Float sum, xl, yl, xm, ym, test, yt, tol, tolmin = 1.e-6, bmax = 20;
   Range x(imax), y(imax), s(2*nemax);
   Float tevz = 0.0253;

   // constant factors
   tol=tolin;
   for ( int i = 0; i < 2*nemax; ++i ){
     if ( i >= int(s.size()) ) { break; }
     s[i] = 0;
   }

   double root1 = (u*sqrt(e)+sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);
   double root2 = (u*sqrt(e)-sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);

   // adaptive calculation of cross section
   sum=0;
   x[0]=0;
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
     x[1]=x[0]; 
     y[1]=y[0];


     do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alphas, betas, 
         sab, az, lasym, teff, sb, sb2, iinc );

     i = 2;
     //std::cout << x[0] << "    " << x[1] << "     " << x[2] << std::endl;
     //std::cout << y[0] << "    " << y[1] << "     " << y[2] << std::endl;
     
     while ( true ){

       do_150( i, imax, y, x, tolmin, tol, teff, xm, ym, yt, e, u, tev, alphas, 
               betas, sab, az, tevz, lasym, lat, sb, sb2, iinc );


       while (true) {
         std::cout << " --- 160 --- " << std::endl;
         //std::cout << (y|ranges::view::all) << std::endl;
         //return;
         j += 1;
         s[2*j+1-1] = x[i-1];
         s[2*j+2-1] = y[i-1];
         s1[j] = x[i-1];
         s2[j] = y[i-1];
         if ( j > 1 ) { sum += (y[i-1]+yl)*(x[i-1]-xl); }
         xl = x[i-1];
         yl = y[i-1];
         if (( j >= nemax-1 ) or 
             ( jbeta > 0 and betas[jbeta-1] > bmax )){ 
             std::cout << " --- 170 --- " << std::endl;
             s[0] = sum;
             s[1] = j;
             s1[0] = sum;
             s2[0] = j;
             return s;
         } 

         i = i - 1;
         if ( i > 1 ){
             //std::cout << " go to 150 " << std::endl;
             break;
         }

         jbeta += 1;
         if ( jbeta <= int(betas.size()) ){
             //std::cout << " go to 111" << std::endl;
             break;
         }
         if ( i == 1 ){ 
             //std::cout << " go to 160 " << std::endl; 
             continue;
         }
    
         std::cout << " --- 170 --- " << std::endl;
         s[0] = sum;
         s[1] = j;
         s1[0] = sum;
         s2[0] = j;

         return s;
      }

      if ( i > 1 ){
        //std::cout << " go to 150 " << std::endl;
        continue;
      }
      else if ( jbeta <= int(betas.size()) ){
        //std::cout << " go to 111" << std::endl;
        break;
      }



    }

  }

   return s;
     /*
     // compare linear approximation to true function
     // 150 continue
     bool goTo150 = true;
     while (true){
       if ( i != imax and goTo150 ){


       if ( do150( i, x, y, xm, ym, yt, test, tolmin, e, u, tev, alphas, betas, 
         sab, tevz, lasym, az, lat, sb, sb2, teff, tol, 
	 iinc) ){continue;}

       //if (do150( i, x, y, xm, ym, yt, test )){ continue; }

         std::cout << std::setprecision(20) << 150 << "     " << y[0] << std::endl;
         if (i <= 3 or 0.5*(y[i-2]+y[i-1])*(x[i-2]-x[i-1]) >= tolmin) {
           xm = 0.5*(x[i-2]+x[i-1]);
           xm = sigfig(xm,8,0);
           if (xm > x[i-1] and xm < x[i-2]){
             ym=0.5*(y[i-2]+y[i-1]);
             yt = sig( e, xm, u, tev, alpha, beta, sab, az, tevz, lasym, 
                 az2, teff2, lat, sb, sb2, teff, iinc );
             test = tol*abs(yt);
      
             if (abs(yt-ym) > test) {
               // point fails
               i=i+1;
               x[i-1]=x[i-2];
               y[i-1]=y[i-2];
               x[i-2]=xm;
               y[i-2]=yt;
               // go to 150
              // return;
               continue;
             }
           }
         }
       }
       goTo150 = true;
  
       // point passes
       // 160 continue
       //std::cout << 160 << "    " << i << "    " << j << "      " << sum << std::endl;
       //if ( j == 220 )return;
       j=j+1;
       s[2*j+1-1]=x[i-1];
       s[2*j+2-1]=y[i-1];
       if (j > 1) sum=sum+(y[i-1]+yl)*(x[i-1]-xl);
       xl=x[i-1];
       yl=y[i-1];

       // if (j >= nemax-1) go to 170
       if (j >= nemax-1 or (jbeta > 0 and beta[jbeta-1] > bmax )) {
         //std::cout << j << "    " << nemax-1 << "    " << jbeta << "     " << bmax <<std::endl;
         s[1-1]=sum;
         s[2-1]=j;
         return; 
      }

      // continue bin loop and linearization loop
      i=i-1;
      //std::cout << i << "     " << jbeta << "    " << beta.size() << std::endl;

      // if (i > 1) go to 150
      if (i > 1) {
        continue;
      }
      //std::cout << "no longer doing 150" << std::endl;
      //
      //
      if ( (unsigned) jbeta > beta.size() and i == 1 ){
        jbeta = jbeta + 1;
        goTo150 = false;
        continue;
      } 

      break;
    } 

    jbeta=jbeta+1;

    // if (jbeta <= nbeta) go to 111
    if ( (unsigned) jbeta <= beta.size()) continue;
 

    // 170 continue
    //std::cout << 170 << std::endl;
    s[1-1]=sum;
    s[2-1]=j;
    return; 
  }

*/
}
#endif

