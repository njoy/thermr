#include "calcem/calcem_util/sig.h"
#include "general_util/sigfig.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include <range/v3/all.hpp>
#include <cmath>
#include <tuple>

template <typename Range, typename Float>
auto initialize_XS_MU_vecs( Range& muVec, Range& xsVec, const Float& e, 
  const Float& ep, const Float& az, const Float& tev, const Range& alphas, 
  const Range& betas, const Range& sab, int lasym, int lat, 
  const Float& sigma_b, const Float& sigma_b2, const Float& teff, int iinc){

  Float beta = std::abs((ep-e)/tev);
  if ( lat == 1 and iinc == 2 ){ beta *= tev/0.0253; }

  muVec[0] =  1.0; 
  muVec[2] = -1.0; 
  muVec[1] = (ep==0.0) ? 0.0 : 
                         0.5 * (e+ep-(pow(1.+beta*beta,0.5)-1.)*az*tev)
                             * pow(e*ep,-0.5);
  if (abs(muVec[1]) > 0.99){ muVec[1] = 0.99; }

  xsVec[0] = sig(e,ep,muVec[0],tev,alphas,betas,sab,az,0.0253,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  xsVec[1] = sig(e,ep,muVec[1],tev,alphas,betas,sab,az,0.0253,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  xsVec[2] = sig(e,ep,muVec[2],tev,alphas,betas,sab,az,0.0253,lasym,lat,sigma_b,sigma_b2,teff,iinc);
}



template <typename Float>
inline Float maxOf4Vals( const Float a, const Float b, const Float c, const Float d ){
  return (a < b) ? (b < c ? (c < d ? d : c ) : (b < d ? d : b)) 
                 : (a < c ? (c < d ? d : c ) : (a < d ? d : a));
}


template <typename Range, typename Float>
auto do_190(Range& muVec, Range& xsVec, Float& xn, Float& xil, Float& muLeft, Float& xsLeft, 
            Float& fract, int& i, Float& gral, int& nL, int& nbin, Range&s, int& j, Float& sum){
  std::cout << " --- 190 --- " << std::endl;
  Float yn = xsLeft + (xsVec[i-1]-xsLeft)*(xn-muLeft)*xil;
  gral = gral + (xn-muLeft)*(xsLeft*0.5*(xn+muLeft) + 
        (xsVec[i-1]-xsLeft)*xil*(-muLeft*0.5*(xn+muLeft)
            + (1.0/3.0)*(xn*xn+xn*muLeft+muLeft*muLeft)));
  Float xbar = gral / fract;

  if ( nL >= 0 ){
    Range p (nL,0.0);
    legndr(xbar,p,nL);
    for (int k = 1; k < nL; ++k){
      s[k] += p[k]/nbin;
    }
  }
  else {
    s[j] = xbar;
  }
  

  muLeft = xn;
  xsLeft = yn;
  sum = 0.0;
  gral = 0.0;
  if (j == nbin){ return 260; }
  if (muLeft < muVec[i-1] ){ return 160; }
  return 250;

}


template <typename Range, typename Float>
auto do_170_175_180(Float& fract, Float& sum, Range& xsVec, Range& muVec, 
  Float& xsLeft, Float& muLeft, int& i, int& j, Float& xil ){
  std::cout << " --- 170 --- " << std::endl;
  Float gral=0.0, xn, yn, xbar;
  j++; 

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
      //std::cout << "to to 190" << std::endl; 
      return xn;
  }
  else if (xn > muLeft and xn < (muVec[i-1]+1e-3*(muVec[i-1]-muLeft))){ 
      std::cout << " --- 180 --- " << std::endl;
      xn = muVec[i-1];
      return xn;
  }
  else { std::cout << "throw std::exception();" << std::endl; return 1.0; }

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
  const Float& az, const int& lasym, const int& lat, 
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
  while ( (unsigned) i < muVec.size() ){ // //std::cout << 110 << std::endl;
      std::cout << " --- 110 --- " << std::endl;
    
    muMid = 0.5*( muVec[i-2] + muVec[i-1] );
    muMid = sigfig(muMid,8,0);
    xs_guess = 0.5*( xsVec[i-2] + xsVec[i-1] );
    xs_true  = sig(e,ep,muMid,tev,alphas,betas,sab,az,0.0253,lasym,lat,sb,sb2,teff,iinc);

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
inline auto sigl(Float ep, Float e, Float tev, Float tolin, int nl, 
  int lat, int iinc, Range& alphas, Range& betas, Range& sab, Float az,
  int lasym, Float sigma_b, Float sigma_b2, Float teff ){
    std::cout.precision (15);

  Float yn = 0.0;
  Float xbar = 0.0;
  using std::abs;
  using std::min;
  using std::pow;

  Float tol = 0.5*tolin;

  int nL = std::abs(nl);
  int nbin = nL - 1;
  
  Float sum = 0, gral = 0;
  Range muVec(20), xsVec(20);


  int i = 3;
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  Float muLeft = muVec[2],
        xsLeft = xsVec[2];
  Float xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);



  do {
    do_110(i, muVec, xsVec,  e,  ep, tev, alphas, betas, sab, az, lasym, 
           lat, sigma_b,  sigma_b2,  teff, iinc, tol,  xsMax);

    do { // If i = 2, then do this action twice. Else do it once
      std::cout << " --- 120 --- " << std::endl;
      sum += 0.5*( xsVec[i-1] + xsLeft )*( muVec[i-1] - muLeft );
      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      i -= 1;
    } while ( i == 1 );

  } while ( i > 1 );



  std::vector<double> s(65,0.0);
  s[0] = (sum <= 1e-32) ? 0.0 : sum;


  std::cout << " --- 130 --- " << std::endl;
  Float fract = sum/(1.0*nbin);
  sum = 0.0;
  int j = 0;


  i = 3;
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  muLeft = muVec[2],
  xsLeft = xsVec[2];
  xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);


  Float xn = 0.0;

  while (true){
    std::cout << " --- 150 ---  "<< std::endl;
    do_110(i, muVec, xsVec,  e,  ep, tev, alphas, betas, sab, az, lasym, 
           lat, sigma_b,  sigma_b2,  teff, iinc, tol,  xsMax);


    while (true){
      std::cout << " --- 160 ---  "<< std::endl;
      Float add = 0.5*(xsVec[i-1]+xsLeft)*(muVec[i-1]-muLeft);

      if (muVec[i-1] == muLeft) { 

          std::cout << " --- 250 ---  "<< std::endl;
          muLeft = muVec[i-1];
          xsLeft = xsVec[i-1];
          i -= 1;
          if ( i > 1 ){ break; }
          if ( i < 1 ){ return s; }
      }

      Float xil = 1.0/(muVec[i-1]-muLeft);
      Float shade = 0.99999999;

      if ( i == 1 and j == nbin-1 ){
          std::cout << " --- 165 ---  "<< std::endl;
          xn = muVec[i-1];
          j++;
          int out = do_190(muVec, xsVec, xn, xil, muLeft, xsLeft, fract, i, gral, nl, nbin, s, j, sum);
          if (j == nbin){ return s; }
          else if ( out == 160 ){
              continue;
          }
          else { // 250
            std::cout << " --- 250 ---  "<< std::endl;
            muLeft = muVec[i-1];
            xsLeft = xsVec[i-1];
            i -= 1;
            if ( i > 1 ){ break; }
            if ( i < 1 ){ return s; }
          }

      }
      else if ( sum + add >= fract * shade and j < nbin - 1 ){
          xn = do_170_175_180(fract, sum, xsVec, muVec, xsLeft, muLeft, i, j, xil );
          int out = do_190(muVec, xsVec, xn, xil, muLeft, xsLeft, fract, i, gral, nl, nbin, s, j, sum);
          if (j == nbin){ return s; }
          else if ( out == 160 ){
              continue;
          }
          else {
            std::cout << " --- 250 ---  "<< std::endl;
            muLeft = muVec[i-1];
            xsLeft = xsVec[i-1];
            i -= 1;
            if ( i > 1 ){ break; }
            if ( i < 1 ){ return s; }
          }
          
      }


      sum  += add;
      gral += 0.5 * (xsLeft*muVec[i-1]-xsVec[i-1]*muLeft) * 
                    (muVec[i-1]+muLeft) + 0.333333333*(xsVec[i-1]-xsLeft) * 
                    (muVec[i-1]*muVec[i-1]+muVec[i-1]*muLeft+muLeft*muLeft);

      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      i -= 1;
      if ( i > 1 ){ break; }
      if ( i < 1 ){ return s; }


    }

  }




}










