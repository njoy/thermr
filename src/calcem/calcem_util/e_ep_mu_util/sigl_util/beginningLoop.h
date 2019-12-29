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
auto populateXSvector(Range& xsVec, Float& xn, Float& invDeltaMu, Float& muLeft, Float& xsLeft, 
  Float& fract, int& i, Float& gral, int& nL, int& nbin, Range&s, const int& j, Float& sum){
  //std::cout << " --- 190 --- " << std::endl;
  gral += (xn-muLeft) * 
          ( xsLeft*0.5*(xn+muLeft) + 
            (xsVec[i-1]-xsLeft)*invDeltaMu*(-muLeft*0.5*(xn+muLeft) + 
            0.333333333*(xn*xn+xn*muLeft+muLeft*muLeft))
          );

  Float xbar = gral / fract;

  if (nL < 0){
    s[j] = xbar;
  }
  else {
    Range p (nL,0.0);
    legndr(xbar,p,nL);
    for (int k = 1; k < nL; ++k){ s[k] += p[k]/nbin; }
  }

  xsLeft = xsLeft + (xsVec[i-1]-xsLeft)*(xn-muLeft)*invDeltaMu;
  muLeft = xn;
  sum  = 0.0;
  gral = 0.0;
}


template <typename Range, typename Float>
auto getNextMuValue(Float& fract, Float& sum, Range& xsVec, Range& muVec, 
  Float& xsLeft, Float& muLeft, int& i, Float& invDeltaMu ){
  //std::cout << " --- 170 --- " << std::endl;
  Float gral=0.0, xn;
  Float test = (fract-sum)*(xsVec[i-1]-xsLeft)/((muVec[i-1]-muLeft)*xsLeft*xsLeft);

  if ( xsLeft >= 1e-32 and abs(test) <= 1e-3 ){ // Could potentially do 180, 190
      xn = muLeft + ( fract-sum )/xsLeft;
      if ( xn > muVec[i-1] ){
        //std::cout << " --- 180 --- " << std::endl;
        return muVec[i-1];
      } 
      else if ( xn >= muLeft ){
        return xn;
      }
  }
  //std::cout << " --- 175 --- " << std::endl;
  Float f = (xsVec[i-1]-xsLeft)*invDeltaMu;
  Float sqrt_disc = pow(pow((xsLeft/f),2)+2.0*(fract-sum)/f,0.5);

  xn = ( f > 0 ) ? muLeft - (xsLeft/f) + sqrt_disc 
                 : muLeft - (xsLeft/f) - sqrt_disc;

  if ( xn > muLeft ){
    if ( xn <= muVec[i-1] ){ return xn; }
    if ( xn < (muVec[i-1]+1e-3*(muVec[i-1]-muLeft))){ return muVec[i-1]; }
  }
  throw std::exception();
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
  const Float& tol, const Float& xsMax){

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
    //std::cout << " --- 110 --- " << std::endl;
    
    muMid = 0.5*( muVec[i-2] + muVec[i-1] );
    muMid = sigfig(muMid,8,0);
    xs_guess = 0.5*( xsVec[i-2] + xsVec[i-1] );
    xs_true  = sig(e,ep,muMid,tev,alphas,betas,sab,az,0.0253,lasym,lat,sb,sb2,teff,iinc);

    if ( ( abs(xs_guess-xs_true) <= tol*abs(xs_true)+tol*xsMax/50.0 and 
           abs(xsVec[i-2]-xsVec[i-1]) <= xs_guess+xsMax/100.0 and 
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

  using std::abs; using std::min; using std::pow;

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
    do_110(i,muVec,xsVec, e, ep,tev,alphas,betas,sab,az,lasym,lat,sigma_b, sigma_b2, teff,iinc,tol, xsMax);
    do { // If i = 2, then do this action twice. Else do it once
      //std::cout << " --- 120 --- " << std::endl;
      sum += 0.5*( xsVec[i-1] + xsLeft )*( muVec[i-1] - muLeft );
      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      i -= 1;
    } while ( i == 1 );

  } while ( i > 1 );



  std::vector<double> s(65,0.0);
  s[0] = (sum <= 1e-32) ? 0.0 : sum;
  //std::cout << " --- 130 --- " << std::endl;
  Float fract = sum/(1.0*nbin);
  sum = 0.0;
  int j = 0;


  i = 3;
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,sigma_b,sigma_b2,teff,iinc);
  muLeft = muVec[2],
  xsLeft = xsVec[2];
  xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);


  Float xn = 0.0;

  do { 
    do_110(i,muVec,xsVec, e, ep,tev,alphas,betas,sab,az,lasym,lat,sigma_b, sigma_b2, teff,iinc,tol, xsMax);

    do {
      //std::cout << " --- 160 ---  "<< std::endl;
      Float add = 0.5*(xsVec[i-1]+xsLeft)*(muVec[i-1]-muLeft);

      if (muVec[i-1] == muLeft) { 
          //std::cout << " --- 250 ---  "<< std::endl;
          muLeft = muVec[i-1];
          xsLeft = xsVec[i-1];
          i -= 1;
          if ( i < 1 ){ return s; }
          if ( i > 1 ){ break; }
      }

      Float invDeltaMu = 1.0/(muVec[i-1]-muLeft);

      if ( i == 1 and j == nbin-1 ){ //std::cout << " --- 165 ---  "<< std::endl;
        j++;
        xn = muVec[i-1];
        populateXSvector(xsVec, xn, invDeltaMu, muLeft, xsLeft, fract, i, gral, nl, nbin, s, j, sum);
        return s;
      }
      else if ( sum + add >= fract * 0.99999999 and j < nbin-1 ){
        j++;
        xn = getNextMuValue(fract, sum, xsVec, muVec, xsLeft, muLeft, i, invDeltaMu );
        populateXSvector(xsVec, xn, invDeltaMu, muLeft, xsLeft, fract, i, gral, nl, nbin, s, j, sum);

        if (muLeft < muVec[i-1] ){ continue; }
        else { // 250
          //std::cout << " --- 250 ---  "<< std::endl;
          muLeft = muVec[i-1];
          xsLeft = xsVec[i-1];
          i -= 1;
          if ( i < 1 ){ return s; }
          if ( i > 1 ){ break; }
        }
      }

      sum  += add;
      gral += 0.5 * (xsLeft*muVec[i-1]-xsVec[i-1]*muLeft) * 
                    (muVec[i-1]+muLeft) + 0.333333333*(xsVec[i-1]-xsLeft) * 
                    (muVec[i-1]*muVec[i-1]+muVec[i-1]*muLeft+muLeft*muLeft);

      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      i -= 1;
      if ( i < 1 ){ return s; }
      if ( i > 1 ){ break; }


    } while ( i <= 1 );

  } while ( i > 1 );


  return s;


}










