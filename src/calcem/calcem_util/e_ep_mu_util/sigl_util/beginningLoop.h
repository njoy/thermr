#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"
#include <range/v3/all.hpp>
#include <cmath>

template <typename Float>
inline Float maxOf4Vals( const Float a, const Float b, const Float c, const Float d ){
  return (a < b) ? (b < c ? (c < d ? d : c ) : (b < d ? d : b)) 
                 : (a < c ? (c < d ? d : c ) : (a < d ? d : a));
}




template <typename Range, typename Float>
auto do_190(Range xsVec, Float xn, Float xil, Float muLeft, Float xsLeft, Float fract, int i){
  std::cout << " --- 190 --- " << std::endl;
  Float yn = xsLeft + (xsVec[i-1]-xsLeft)*(xn-muLeft)*xil;
  Float gral = gral + (xn-muLeft)*(xsLeft*0.5*(xn+muLeft) + 
        (xsVec[i-1]-xsLeft)*xil*(-muLeft*0.5*(xn+muLeft)
            + (1.0/3.0)*(xn*xn+xn*muLeft+muLeft*muLeft)));
  Float xbar = gral / fract;
  return std::make_tuple<double,double,double>(yn,gral,xbar);
}


template <typename Range, typename Float>
auto do_170_175_180(Float& fract, Float& sum, Range& xsVec, Range& muVec, 
  Float& xsLeft, Float& muLeft, int& i, int& j, Float& xil ){

  std::cout << " --- 170 --- " << std::endl;
  j++; 
  Float gral=0.0, xn, yn, xbar;

  Float test = (fract-sum)*(xsVec[i-1]-xsLeft)/((muVec[i-1]-muLeft)*xsLeft*xsLeft);
  if ( xsLeft >= 1e-32 and abs(test) <= 1e-3 ){
      // Could potentially do 180, 190
      Float xn = muLeft + ( fract-sum)/xsLeft;
      if ( xn > muVec[i-1] ){
        std::cout << " --- 180 --- " << std::endl;
        xn = muVec[i-1];
        return xn;

      } 
      if ( xn >= muLeft and xn <= muVec[i-1] ){
        std::cout << " --- 190 --- " << std::endl;
        return xn;
      }

  }
  std::cout << " --- 175 --- " << std::endl;
  Float f = (xsVec[i-1]-xsLeft)*xil;
  Float disc = (xsLeft/f)*(xsLeft/f)+2.0*(fract-sum)/f;
  if ( disc < 0 ){ throw std::exception(); }
  if (f > 0){ xn = muLeft - (xsLeft/f) + pow(disc,0.5); }
  if (f < 0){ xn = muLeft - (xsLeft/f) - pow(disc,0.5); }

  if (xn > muLeft and xn <= muVec[i-1]){ 
      std::cout << "to to 190" << std::endl; 
      return xn;
  }
  else if (xn > muLeft and xn < (muVec[i-1]+1e-3*(muVec[i-1]-muLeft))){ 
      std::cout << " --- 180 --- " << std::endl;
      xn = muVec[i-1];
      return xn;
  }
  else { throw std::exception(); }


}



template <typename Range, typename Float>
inline auto adaptiveLinearization( Range& x, Range& y, const Float& e, 
  const Float& ep, const Float& tev, const Float& tevz, const Range& alphas, 
  const Range& betas, const Range& sab, const Float& az, 
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
  y[2] = sig(e,ep,x[2],tev,alphas,betas,sab,az,tevz,lasym,/*az2,teff2,*/lat,sb,sb2,teff,iinc);

  // Consider a cosine mu that corresponds to an alpha value of sqrt(1+beta^2).
  // What's the cross section?
  x[1] = 0.5 * seep * ( e + ep - (s1bb-1) * az * tev );
  if (std::abs(x[1]) > 1-eps) x[1] = 0.99;
  x[1] = sigfig(x[1],8,0);
  y[1] = sig(e,ep,x[1],tev,alphas,betas,sab,az,tevz,lasym,/*az2,teff2,*/lat,sb,sb2,teff,iinc);

  // Consider a cosine mu equal to 1. What's the cross section?
  x[0] = 1;
  y[0] = sig(e,ep,x[0],tev,alphas,betas,sab,az,tevz,lasym,/*az2,teff2,*/lat,sb,sb2,teff,iinc);

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
  const Float& tev, const Range& alphas, const Range& betas, const Range& sab, 
  const Float& az, const Float tevz, const int& lasym, const int& lat, 
  const Float& sb, const Float& sb2, const Float& teff, const int& iinc, 
  const Float& tol, const Float& ymax){

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
    xs_true  = sig(e,ep,muMid,tev,alphas,betas,sab,az,tevz,lasym,lat,sb,sb2,teff,iinc);
    
    if ( ( abs(xs_guess-xs_true) <= tol*abs(xs_true)+tol*ymax/50.0 and 
           abs(xsVec[i-2]-xsVec[i-1]) <= xs_guess+ymax/100.0 and 
           (muVec[i-2]-muVec[i-1]) < 0.5 ) or
          ( muVec[i-2]-muVec[i-1] < 1e-5 ) ) { 
      return; 
    }

 
    shiftOver( i, muVec, xsVec, muMid, xs_true );
  }  
} 




template <typename Range, typename Float>
inline auto do_things(Float ep, Float e, Float tev, Float tolin, int nL, 
  int lat, int iinc, Range& alphas, Range& betas, Range& sab, Float az,
  int lasym, Float sigma_b, Float sigma_b2, Float teff ){

  using std::abs;
  using std::min;
  using std::pow;

  Float beta = std::abs((ep-e)/tev);
  Float tol = 0.5*tolin;

  nL = std::abs(nL);
  
  Float tevz = 0.0253;
  if ( lat == 1 and iinc == 2 ){ beta *= tev/tevz; }

  int i = 3;
  Float sum = 0;

  Range muVec(20), xsVec(20);

  muVec[0] =  1.0; 
  muVec[2] = -1.0; 
  muVec[1] = (ep==0.0) ? 0.0 : 
                         0.5 * (e+ep-(pow(1.+beta*beta,0.5)-1.)*az*tev)
                             * pow(e*ep,-0.5);
  if (abs(muVec[1]) > 0.99){ muVec[1] = 0.99; }

  xsVec[0] = sig(e,ep,muVec[0],tev,alphas,betas,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  xsVec[1] = sig(e,ep,muVec[1],tev,alphas,betas,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  xsVec[2] = sig(e,ep,muVec[2],tev,alphas,betas,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc);


  Float muLeft = muVec[2],
        xsLeft = xsVec[2];

  Float xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);



  while ( true ){
    do_110(i, muVec, xsVec,  e,  ep, tev, alphas, betas, sab, az, tevz, lasym, 
           lat, sigma_b,  sigma_b2,  teff, iinc, tol,  xsMax);

    // 200                      
    while ( true ){
      sum += 0.5*( xsVec[i-1] + xsLeft )*(muVec[i-1]-muLeft);
      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      i -= 1;
    
      if ( i != 1 ){ break; }
    }

    if ( i <= 1 ){ break; }
  }



  std::vector<double> s(65,0.0);

  s[0] = sum;
  if ( sum <= 1e-32 ){
    for (int i = 0; i < nL; ++i){
        s[i] = 0.0;
    }
  }

  // 130
  int nbin = nL - 1;
  Float fract = sum/(1.0*nbin);
  sum = 0.0;
  Float gral = 0.0;
  for (int i = 1; i < nL; ++i){
    s[i] = 0.0;
  }
  int j = 0;



  //---------------------------------------------------------------------------
  i = 3;
  //---------------------------------------------------------------------------
  muVec[0] =  1.0; 
  muVec[2] = -1.0; 
  muVec[1] = (ep==0.0) ? 0.0 : 
                         0.5 * (e+ep-(pow(1.+beta*beta,0.5)-1.)*az*tev)
                             * pow(e*ep,-0.5);
  if (abs(muVec[1]) > 0.99){ muVec[1] = 0.99; }

  xsVec[0] = sig(e,ep,muVec[0],tev,alphas,betas,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  xsVec[1] = sig(e,ep,muVec[1],tev,alphas,betas,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  xsVec[2] = sig(e,ep,muVec[2],tev,alphas,betas,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc);


  muLeft = muVec[2],
  xsLeft = xsVec[2];

  xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);
  std::cout << (xsVec|ranges::view::all) << std::endl;
  std::cout << (muVec|ranges::view::all) << std::endl;

  Float xn=0.0;

  while (true){
    std::cout << " --- 150 ---  "<< std::endl;
    // 150
    do_110(i, muVec, xsVec,  e,  ep, tev, alphas, betas, sab, az, tevz, lasym, 
           lat, sigma_b,  sigma_b2,  teff, iinc, tol,  xsMax);


    while (true){
      std::cout << " --- 160 ---  "<< std::endl;
      // 160
      Float add = 0.5*(xsVec[i-1]+xsLeft)*(muVec[i-1]-muLeft);


      if (muVec[i-1] == muLeft) {  // 250
          std::cout << " --- 250 ---  "<< std::endl;
          //std::cout << "go to 250" << std::endl; 
          muLeft = muVec[i-1];
          xsLeft = xsVec[i-1];
          i -= 1;
          if ( i > 1 ){ break; }
          if ( i != 1 ){ return; }
      }

      Float xil = 1.0/(muVec[i-1]-muLeft);
      Float shade = 0.99999999;

      if ( i == 1 and j == nbin-1 ){
          std::cout << " --- 165 ---  "<< std::endl;
          xn = muVec[i-1];
          j++;
          std::cout << "go to 190" << std::endl;
          do_190(xsVec, xn, xil, muLeft, xsLeft, fract, i);
          return ;
      }
      else if ( sum + add >= fract * shade and j < nbin - 1 ){
          xn = do_170_175_180(fract, sum, xsVec, muVec, xsLeft, muLeft, i, j, xil );
          do_190(xsVec, xn, xil, muLeft, xsLeft, fract, i);
      }

      sum += add;
      gral = gral + 0.5 * (xsLeft*muVec[i-1]-xsVec[i-1]*muLeft)*(muVec[i-1]+muLeft)+0.333333*(xsVec[i-1]-xsLeft)*(muVec[i-1]*muVec[i-1]+muVec[i-1]*muLeft+muLeft*muLeft);
      std::cout << " --- 250 ---  "<< std::endl;
      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      i -= 1;
      if ( i > 1 ){ break; }
      if ( i != 1 ){ return; }



    }

  }


  //std::cout << (xsVec|ranges::view::all) << std::endl;
  //std::cout << (muVec|ranges::view::all) << std::endl;






  //std::cout << nbin << std::endl;

  //std::cout << muVec[0] << "   " << muVec[1] << "    " << muVec[2] << std::endl;


  //std::cout << i << std::endl;
  //std::cout << sum << std::endl;

  //std::cout << (muVec|ranges::view::all) << std::endl;
  //std::cout << std::endl;
  //std::cout << (xsVec|ranges::view::all) << std::endl;



  return;
  std::cout << sigma_b << "   " << sigma_b2<<std::endl;
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

