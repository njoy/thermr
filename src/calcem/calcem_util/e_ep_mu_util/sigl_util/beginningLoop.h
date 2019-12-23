#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"
#include <cmath>

template <typename Float>
inline Float maxOf3Vals( const Float& a, const Float& b, const Float& c ){
  return (a < b) ? (b < c ? c : b) : (a < c ? c : a);
}


template <typename Range, typename Float>
inline auto adaptiveLinearization( Range& x, Range& y, const Float& e, 
  const Float& ep, const Float& tev, const Float& tevz, const Range& alpha, 
  const Range& beta,const Range& sab, const Float& az, 
  const int& lasym, const Float& teff, const int& lat, 
  const Float& sb, const Float& sb2, const int& iinc, const Float& eps, 
  const Float& seep, const Float& s1bb  ){
  /* So here we consider three angles - a cosine value of -1, a cosine value 
   * that corresponds to alpha = sqrt(1+beta^2) [According to Eq. 227], and a
   * cosine value of 1. We calculate the incoherent cross sections for all 
   * of these angles (using sig), and determine the maximum cross section of 
   * these three.
   */

  // prime stack for equally-probable angles
  //std::cout << 130 << std::endl;
  // adaptive linearization
  // Consider a cosine mu equal to -1. What's the cross section?
  x[2] = -1;
  y[2] = sig(e,ep,x[2],tev,alpha,beta,sab,az,tevz,lasym,/*az2,teff2,*/lat,sb,sb2,teff,iinc);

  // Consider a cosine mu that corresponds to an alpha value of sqrt(1+beta^2).
  // What's the cross section?
  x[1] = 0.5 * seep * ( e + ep - (s1bb-1) * az * tev );
  if (std::abs(x[1]) > 1-eps) x[1] = 0.99;
  x[1] = sigfig(x[1],8,0);
  y[1] = sig(e,ep,x[1],tev,alpha,beta,sab,az,tevz,lasym,/*az2,teff2,*/lat,sb,sb2,teff,iinc);

  // Consider a cosine mu equal to 1. What's the cross section?
  x[0] = 1;
  y[0] = sig(e,ep,x[0],tev,alpha,beta,sab,az,tevz,lasym,/*az2,teff2,*/lat,sb,sb2,teff,iinc);

  Float ymax = maxOf3Vals(y[0],y[1],y[2]);
  //if ( e >= 1.05 and e < 1.050001 and ep > 2.621273e-2 and ep < 2.621274e-2 )std::cout << "in sigl     " << y[1] << "     " << y[2] << std::endl;
  return ( ymax < eps ) ? eps : ymax;

}



template <typename Range, typename Float>
inline void shiftOver( int& i, Range& x, Range& y, Float& muMid, const Float& yt ){
  // x = [   mu1      mu2          mu3            0    0   0 ... ]
  // y = [ s(mu1)   s(mu2)       s(mu3)           0    0   0 ... ]
  //  becomes 
  // x = [   mu1      mu2      0.5*(mu2+mu3)     mu3   0   0 ... ]
  // y = [ s(mu1)   s(mu2)   s(0.5*(mu2+mu3))  s(mu3)  0   0 ... ]
  i++;
  x[i-1] = x[i-2]; x[i-2] = muMid;
  y[i-1] = y[i-2]; y[i-2] = yt;
}

template <typename Range, typename Float>
inline auto do_110(int& i, Range& muVec, Range& xsVec, const Float& e, const Float& ep, 
  const Float& tev, const Range& alpha, const Range& beta, const Range& sab, 
  const Float& az, const Float tevz, const int& lasym, const int& lat, 
  const Float& sb, const Float& sb2, const Float& teff, const int& iinc, 
  const Float& mu_tol, const Float& tol, const Float& ymax){

  // If the spacing isn't fine enough, we'll bisect the grid
  //
  // x = [   mu1      mu2          mu3            0    0   0 ... ]
  // y = [ s(mu1)   s(mu2)       s(mu3)           0    0   0 ... ]
  //
  //           will be turned into
  //
  // x = [   mu1      mu2      0.5*(mu2+mu3)     mu3   0   0 ... ]
  // y = [ s(mu1)   s(mu2)   s(0.5*(mu2+mu3))  s(mu3)  0   0 ... ]
  //
  // where s(x) is the incoherent cross section [Eq. 225] evaluated at 
  // cosine x
   
  using std::abs;
  Float muMid, xs_guess, xs_true;
  while ( (unsigned) i < muVec.size() ){ // std::cout << 110 << std::endl;
    
    muMid = 0.5*( muVec[i-2] + muVec[i-1] );
    muMid = sigfig(muMid,8,0);
    xs_guess = 0.5*( xsVec[i-2] + xsVec[i-1] );
    xs_true  = sig(e,ep,muMid,tev,alpha,beta,sab,az,tevz,lasym,lat,sb,sb2,teff,iinc);
    
    if ( ( abs(xs_guess-xs_true) <= tol*abs(xs_true)+tol*ymax/50.0 and 
           abs(xsVec[i-2]-xsVec[i-1]) <= xs_guess+ymax/100.0 and 
           (muVec[i-2]-muVec[i-1]) < 0.5 ) or
          ( muVec[i-2]-muVec[i-1] < mu_tol ) ) { 
      return; 
    }

 
    shiftOver( i, muVec, xsVec, muMid, xs_true );
  }  
} 




template <typename Range, typename Float>
inline auto do_things(Float ep, Float e, Float tev, Float tolin, int nL, 
  int lat, int iinc  ){


  Float beta = std::abs((ep-e)/tev);
  Float tol = 0.5*tolin;

  nL = std::abs(nL);
  
  Float tevz = 0.0253;
  if ( lat == 1 and iinc == 2 ){ beta *= tev/tevz; }

  int i = 3;
  Float sum = 0;

  Range muVec(20), xsVec(20);

  muVec[2] = -1.0; xsVec[2] = sig(e,ep,x[2],tev,alpha,beta,sab); 
  //muVec[1] =  1.0; xsVec[1] =  1.0;


  




  //Float xsLeft = xsVec[3];
  //do_110(i, muVec, xsVec,  e,  ep, tev, alpha, beta, sab, az, tevz, lasym, 
  //       lat, sb,  sb2,  teff, iinc, mu_tol,  tol,  ymax);
  
  //sum = sum + 0.5*(xsVec[i]);




}





/*
*/













template <typename Range, typename Float>
inline auto do_110_120_130_for_sigl( int& i, Range& x, Range& y, 
  const Float& e, const Float& ep, const Float& tev, const Float& tevz, 
  const Range& alpha, const Range& beta,
  const Range& sab, 
  const Float& az, const int& lasym, const Float& teff, 
  const int& lat, const Float& sb, 
  const Float& sb2, const int& iinc, const int& nl, 
  const Float& sigmin, Range& s, int& nbin, Float& fract, 
  Float& xl, int& j, Float& ymax, 
  Float& yl, const Float& tol, 
  const Float& xtol ){
  /* For this, we fill up mu values into x, and S(a,b,mu) values into y (with
   * a and b being fixed, and mu corresponding to the values in x). The grid
   * is chosen adaptively. So let's consider x = [2, -8, -16, 0, 0, 0, ... ]
   * (in reality 2-->1 and -16 --> -1, but this is just good for discussion now)
   */



  Float gral, sum = 0;

  while (true){

    // Fills up x, y with mu and S(a,b,mu) values (respectively) so that they
    // are close enough to be reasonably interpolated.
    // x = [   mu1      mu2      mu3     ...    mu_i     0    0   0 ... ]
    // y = [ s(mu1)   s(mu2)   s(mu3)    ...  s(mu_i)    0    0   0 ... ]

    do_110(i, x, y, e, ep, tev, alpha, beta, sab, az, tevz, lasym, /*az2, 
        teff2,*/ lat, sb, sb2, teff, iinc, xtol, tol, ymax);

    // When do_100 returns, we x and y both have i-many nonzero entries
    // On the first iteration, xl = -1, and yl = S(a,b,mu=-1)

    while ( true ){

      // std::cout << 120 << std::endl;
      
      // Sum is meant to be the integral of cross section across all values of
      // mu. 
      // x = [   -1  ..   mu_j  ..   1   0  0  0 ... ]
      // y = [ s(-1) .. s(mu_j) .. s(1)  0  0  0 ... ]
      
      sum = sum + 0.5*(y[i-1]+yl)*(x[i-1]-xl);
      // If the tolerance between my i-th point and my (i-1)th point is 
      // reasonably low, then go left to check the other side.
      xl  = x[i-1];
      yl  = y[i-1];
      i   = i - 1;

      if ( i > 1 ){ break; }
      // if (i == 1) go to 120
      
      if (i != 1) { // don't go to 120 - either go to 130 or return 
        
        // if (sum > sigmin) go to 130
        if (sum > sigmin) { 
          s[0] = sum; 
          nbin  = nl - 1;
          fract = sum/nbin;
          sum   = 0;
          i     = 3;
          j     = 0;
          xl    = -1;
          gral  = 0;
          for ( int il = 1; il < nl; ++il ){ s[il] = 0; } 
          return std::tuple<Float,Float,bool> { gral, sum, false };
        } 
        else { 
          s[0] = 0; 
          for ( int il = 1; il < nl; ++il ){ s[il] = 0; } 
          return std::tuple<Float,Float,bool> { gral, sum, true };
        }


      } // don't go to 120
    } // go to 120
  } // go to 110

}

