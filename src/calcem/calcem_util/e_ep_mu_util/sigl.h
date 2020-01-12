#include "calcem/calcem_util/sig.h"
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

  muVec[2] = -1.0; 
  muVec[1] = (ep==0.0) ? 
      0.0 
    : std::min(0.5 * (e+ep-(pow(1.+beta*beta,0.5)-1.)*az*tev) * pow(e*ep,-0.5), 0.99);
  muVec[0] =  1.0; 
  muVec[1] = 0.0;

  xsVec[0] = sig(e,ep,muVec[0],tev,alphas,betas,sab,az,0.0253,lasym,lat,sigma_b,
                 sigma_b2,teff,iinc);
  xsVec[1] = sig(e,ep,muVec[1],tev,alphas,betas,sab,az,0.0253,lasym,lat,sigma_b,
                 sigma_b2,teff,iinc);
  xsVec[2] = sig(e,ep,muVec[2],tev,alphas,betas,sab,az,0.0253,lasym,lat,sigma_b,
                 sigma_b2,teff,iinc);
}

template <typename Float>
inline Float maxOf4Vals( const Float a, const Float b, const Float c, const Float d ){
  return (a < b) ? (b < c ? (c < d ? d : c ) : (b < d ? d : b)) 
                 : (a < c ? (c < d ? d : c ) : (a < d ? d : a));
}

template <typename Range, typename Float>
auto populateXSvector(Range& xsVec, Float& xn, Float& invDeltaMu, Float& muLeft, 
  Float& xsLeft, Float& fract, int& i, Float& gral, const int& nbin, Range& s, 
  const int& j, bool equiprobableBins){
  //std::cout << " --- 190 --- " << std::endl;
  gral += (xn-muLeft) * 
          ( xsLeft*0.5*(xn+muLeft) + 
            (xsVec[i-1]-xsLeft)*invDeltaMu*(-muLeft*0.5*(xn+muLeft) + 
            0.333333333*(xn*xn+xn*muLeft+muLeft*muLeft))
          );

  if (equiprobableBins){
    s[j-1] = gral / fract;  // xbar
  }
  else {        // Here we fill in the legendre components
    Range p (nbin+1,0.0);
    legndr(gral/fract,p,nbin+1);
    for (int k = 1; k < nbin+1; ++k){ s[k-1] += p[k]/nbin; }
  }

  xsLeft = xsLeft + (xsVec[i-1]-xsLeft)*(xn-muLeft)*invDeltaMu;
  muLeft = xn;
  gral = 0.0;
}


template <typename Range, typename Float>
auto getNextMuValue(Float& fract, const Float& sum, Range& xsVec, Range& muVec, 
  Float& xsLeft, Float& muLeft, int& i, Float& invDeltaMu ){
  //std::cout << " --- 170 --- " << std::endl;
  Float xn;
  Float test = (fract-sum)*(xsVec[i-1]-xsLeft)/((muVec[i-1]-muLeft)*xsLeft*xsLeft);

  if ( xsLeft >= 1e-32 and abs(test) <= 1e-3 ){ // Could potentially do 180, 190
      xn = muLeft + (fract-sum)/xsLeft;
      if ( muVec[i-1] < xn ){ //std::cout << " --- 180 --- " << std::endl;
        return muVec[i-1];
      } 
      else if ( muLeft <= xn ){
        return xn;
      }
  }
  //std::cout << " --- 175 --- " << std::endl;
  Float f = (xsVec[i-1]-xsLeft)*invDeltaMu;
  Float sqrt_disc = pow(pow((xsLeft/f),2)+2.0*(fract-sum)/f,0.5);

  xn = ( f > 0 ) ? muLeft - (xsLeft/f) + sqrt_disc 
                 : muLeft - (xsLeft/f) - sqrt_disc;

  if ( xn <= muLeft ){ throw std::exception(); }

  if ( xn <= muVec[i-1] ){ return xn; }
  if ( xn < (muVec[i-1]+1e-3*(muVec[i-1]-muLeft))){ return muVec[i-1]; }
  else { throw std::exception(); }

}



template <typename Range, typename Float>
inline void shiftOver( int& i, Range& muVec, Range& xsVec, Float& muMid, const Float& xsTrue ){
  // Inserts the xsTrue value adn muMid into their respective vectors. Increase
  // i by one so that now we can check between midpoint and i-1 to see if we 
  // need any additional points in there.
  i++;
  muVec[i-1] = muVec[i-2]; muVec[i-2] = muMid;
  xsVec[i-1] = xsVec[i-2]; xsVec[i-2] = xsTrue;
}

template <typename Range, typename Float>
inline auto addMidpointsRight(int& i, Range& muVec, Range& xsVec, const Float& e, 
  const Float& ep, const Float& tev, const Range& alphas, const Range& betas, 
  const Range& sab, const Float& az, const int& lasym, const int& lat, 
  const Float& sb, const Float& sb2, const Float& teff, const int& iinc, 
  const Float& tol){
  using std::abs;

  Float xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);
  Float muMid, xs_guess, xs_true;
  while ( (unsigned) i < muVec.size() ){ //std::cout << " --- 110 --- " << std::endl;
    muMid    = 0.5*( muVec[i-2] + muVec[i-1] ); //muMid = sigfig(muMid,8,0);
    xs_guess = 0.5*( xsVec[i-2] + xsVec[i-1] );
    xs_true  = sig(e,ep,muMid,tev,alphas,betas,sab,az,0.0253,lasym,lat,sb,sb2,teff,iinc);

    if ( ( abs(xs_guess-xs_true) <= tol*abs(xs_true)+tol*xsMax/50.0 and 
           abs(xsVec[i-2]-xsVec[i-1]) <= xs_guess+xsMax/100.0 and 
           (muVec[i-2]-muVec[i-1]) < 0.5 ) or
         ( muVec[i-2]-muVec[i-1] < 1e-5 ) ) { break; }
    shiftOver( i, muVec, xsVec, muMid, xs_true );
    std::cout << "       " << (muVec|ranges::view::all) << std::endl;
  }  
} 


template <typename Range, typename Float>
inline auto getPDF(Float ep, Float e, Float tev, Float tol, int lat, int iinc, 
  const Range& alphas, const Range& betas, const Range& sab, Float az,
  int lasym, Float sigma_b, Float sigma_b2, Float teff ){
  // This computes integral -1 -> +1 of sigma(E->E',mu) dmu
  // through trapezoidal integration

  using std::abs; using std::min; using std::pow;
  Float pdf = 0;
  Range muVec(20,0.0), xsVec(20,0.0);
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,
                        sigma_b,sigma_b2,teff,iinc);
  Float muLeft = muVec[2], xsLeft = xsVec[2];
  // The outer loop will check between i-2 and i-1 to see if we need a midpoint
  // in between there. if so, it'll put it in there and then check between 
  // midpoint and i-1. If we don't need a midpoint, we go to the inner loop
  // and move to the left and keep checking. 
  int i = 3;
  do {
    addMidpointsRight(i,muVec,xsVec,e,ep,tev,alphas,betas,sab,az,lasym,lat,
                      sigma_b,sigma_b2,teff,iinc,tol);
    do { // If i = 2, then do this action twice. Else do it once
      pdf += 0.5*( xsVec[i-1] + xsLeft )*( muVec[i-1] - muLeft );
      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      --i;
    } while ( i == 1 );
  } while ( i > 1 );
  return pdf;
}

template <typename Range, typename Float>
inline auto sigl(Float ep, Float e, Float tev, Float tol, int lat, int iinc, 
  const Range& alphas, const Range& betas, const Range& sab, Float az,
  int lasym, Float sigma_b, Float sigma_b2, Float teff, const int nbin, 
  bool equiprobableBins = true ){

  using std::abs; using std::min; using std::pow;

  Range muVec(20), xsVec(20), s(nbin,0.0);

  tol *= 0.5;

  Float pdf = getPDF(ep, e, tev, tol, lat, iinc, alphas, betas, sab, az, lasym, 
                     sigma_b, sigma_b2, teff);

  if ( pdf <= 1e-32 ){ return s; }

  //std::cout << " --- 130 --- " << std::endl;
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,
                        sigma_b,sigma_b2,teff,iinc);
  Float muLeft = muVec[2],
        xsLeft = xsVec[2],
        gral   = 0.0,
        sum    = 0.0,
        // This ``fract'' is how much the integrated xs should be (roughly) 
        // in each angle bin. 
        fract  = pdf/(1.0*nbin);

  int i = 3, j = 0;

  do { 
    addMidpointsRight(i,muVec,xsVec,e,ep,tev,alphas,betas,sab,az,lasym,lat,
                      sigma_b,sigma_b2,teff,iinc,tol);
    do { //std::cout << " --- 160 ---  "<< std::endl;
      if (muVec[i-1] == muLeft) {  
        // Is this our first time in the outer do loop? i.e. are we trying to
        // compare -1 with -1? If so, let's just move i down one and try again.
        //std::cout << " --- 250 ---  "<< std::endl;
        --i;
        if ( i < 1 ){ return s; }
        if ( i > 1 ){ break; }
      }

      std::cout << "now I'm confident with " << muLeft << "  and   " << muVec[i-1] << "      " << i <<  "   " << j << std::endl;
      Float invDeltaMu = 1.0/(muVec[i-1]-muLeft);

      // This is the integrated little block of xs (from muLeft -> muVec[i-1]).
      // We're going to want to see if adding this to the sum will make us reach
      // our desired xs from fract (recall that fract is the xs we will have 
      // in each bin if they're perfectly equiprobable).
      Float add = 0.5*(xsVec[i-1]+xsLeft)*(muVec[i-1]-muLeft);

      if ( i == 1 and j == nbin-1 ){ //std::cout << " --- 165 ---  "<< std::endl;
        j++;
        Float xn = muVec[i-1];
        populateXSvector(xsVec, xn, invDeltaMu, muLeft, xsLeft, fract, i, gral, 
                         nbin, s, j, equiprobableBins);
        return s;
      }
      else if ( sum + add >= fract * 0.99999999 and j < nbin-1 ){
        // Here, we're trying to get sum+add as close to possible as fract. 
        // Fract is equalt to the cross section (integrated over all mu) divided
        // by the number of bins. So we're going to keep on adding to this sum
        // tally and giving it more and more bin bits until it looks like its 
        // about to overflow.
        j++;
        Float xn = getNextMuValue(fract, sum, xsVec, muVec, xsLeft, muLeft, i, invDeltaMu );
        populateXSvector(xsVec, xn, invDeltaMu, muLeft, xsLeft, fract, i, gral, nbin, s, j, equiprobableBins);
        sum = 0.0;

        if (muLeft < muVec[i-1] ){ continue; }
        else {  //std::cout << " --- 250 ---  "<< std::endl;
          muLeft = muVec[i-1];
          xsLeft = xsVec[i-1];
          --i;
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
      --i;

    } while ( i == 1 );

  } while ( i > 1 );

  return s;
}











