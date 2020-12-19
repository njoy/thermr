#include "inelastic/sig.h"
#include "coherentElastic/coherentElastic_util/sigcoh_util/legndr.h"
#include <range/v3/all.hpp>
#include <cmath>
#include <tuple>



template <typename Range, typename Float>
auto initialize_XS_MU_vecs( Range& muVec, Range& xsVec, const Float& e, 
  const Float& ep, const Float& az, const Float& tev, const Range& alphas, 
  const Range& betas, const Range& sab, int lasym, int lat, 
  const Range& boundXsVec, const Float& teff, int iinc){

  Float beta = std::fabs((ep-e)/tev);
  if ( lat == 1 and iinc == 2 ){ beta *= tev/0.0253; }

  muVec[2] = -1.0; 
  muVec[1] = (ep==0.0) ? 
      0.0 
    : 0.5 * (e+ep-(std::pow(1.+beta*beta,0.5)-1.)*az*tev) * std::pow(e*ep,-0.5);
  if (std::fabs(muVec[1]) > 0.999){ muVec[1] = 0.99; }
  muVec[0] =  1.0; 

  xsVec[0] = sig(e,ep,muVec[0],tev,alphas,betas,sab,az,0.0253,lasym,lat,boundXsVec,
                 teff,iinc);
  xsVec[1] = sig(e,ep,muVec[1],tev,alphas,betas,sab,az,0.0253,lasym,lat,boundXsVec,
                 teff,iinc);
  xsVec[2] = sig(e,ep,muVec[2],tev,alphas,betas,sab,az,0.0253,lasym,lat,boundXsVec,
                 teff,iinc);
}

template <typename Float>
inline Float maxOf4Vals( const Float a, const Float b, const Float c, const Float d ){
  return (a < b) ? (b < c ? (c < d ? d : c ) : (b < d ? d : b)) 
                 : (a < c ? (c < d ? d : c ) : (a < d ? d : a));
}

template <typename Range, typename Float>
auto populateXSvector(Range& xsVec, Float& mu, Float& invDeltaMu, Float& muLeft, 
  Float& xsLeft, Float& xs_per_bin, int& i, Float& gral, const int& nbin, Range& s, 
  const int& j, bool equiprobableBins){ //std::cout << " --- 190 --- " << std::endl;
  gral += (mu-muLeft) * 
          ( xsLeft*0.5*(mu+muLeft) + 
            (xsVec[i-1]-xsLeft)*invDeltaMu*(-muLeft*0.5*(mu+muLeft) + 
            0.333333333*(mu*mu+mu*muLeft+muLeft*muLeft))
          );

  if (equiprobableBins){
    s[j-1] = gral / xs_per_bin;  // xbar
  }
  else {        // Here we fill in the legendre components
    Range p (nbin+1,0.0);
    legndr(gral/xs_per_bin,p,nbin+1);
    for (int k = 1; k < nbin+1; ++k){ s[k-1] += p[k]/nbin; }
  }

  xsLeft = xsLeft + (xsVec[i-1]-xsLeft)*(mu-muLeft)*invDeltaMu;
  muLeft = mu;
  gral = 0.0;
}


template <typename Range, typename Float>
auto getNextMuValue(Float& xs_per_bin, const Float& sum, Range& xsVec, 
  Float& xsLeft, Float& muLeft, Float& muRight, int& i, Float& invDeltaMu ){
  Float mu, //std::cout << " --- 170 --- " << std::endl;
    test = (xs_per_bin-sum)*(xsVec[i-1]-xsLeft)/((muRight-muLeft)*xsLeft*xsLeft);

  if ( xsLeft >= 1e-32 and std::fabs(test) <= 1e-3 ){ // Could potentially do 180, 190
      mu = muLeft + (xs_per_bin-sum)/xsLeft; //std::cout << " --- 180 --- " << std::endl;
      if      ( muRight <  mu ){ return muRight; } 
      else if ( muLeft  <= mu ){ return mu;      }
  }
  //std::cout << " --- 175 --- " << std::endl;
  Float f = (xsVec[i-1]-xsLeft)*invDeltaMu;
  Float sqrt_disc = std::pow(std::pow((xsLeft/f),2)+2.0*(xs_per_bin-sum)/f,0.5);

  mu = ( f > 0 ) ? muLeft - (xsLeft/f) + sqrt_disc 
                 : muLeft - (xsLeft/f) - sqrt_disc;

  if ( mu <= muLeft ){ throw std::exception(); }

  if ( mu <= muRight ){ return mu; }
  if ( mu < (muRight+1e-3*(muRight-muLeft))){ return muRight; }
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
  const Range& boundXsVec, const Float& teff, const int& iinc, 
  const Float& tol, const Float& xsMax ){

  Float muMid, xs_guess, xs_true;
  while ( (unsigned) i < muVec.size() ){ //std::cout << " --- 110 --- " << std::endl;
    muMid    = 0.5*( muVec[i-2] + muVec[i-1] ); muMid = sigfig(muMid,8,0);
    xs_guess = 0.5*( xsVec[i-2] + xsVec[i-1] );
    xs_true  = sig(e,ep,muMid,tev,alphas,betas,sab,az,0.0253,lasym,lat,boundXsVec,teff,iinc);

    if ( ( std::fabs(xs_guess-xs_true) <= tol*std::fabs(xs_true)+tol*xsMax/50.0 and 
           std::fabs(xsVec[i-2]-xsVec[i-1]) <= xs_guess+xsMax/100.0 and 
           (muVec[i-2]-muVec[i-1]) < 0.5 ) or
         ( muVec[i-2]-muVec[i-1] < 1e-5 ) ) { break; }
    shiftOver( i, muVec, xsVec, muMid, xs_true );
  }  
} 



template <typename Range, typename Float>
inline auto getPDF(Float ep, Float e, Float tev, Float tol, int lat, int iinc, 
  const Range& alphas, const Range& betas, const Range& sab, Float az,
  int lasym, const Range& boundXsVec, Float teff ){
  // This computes integral -1 -> +1 of sigma(E->E',mu) d mu
  // through trapezoidal integration

  Float pdf = 0;
  Range muVec(20,0.0), xsVec(20,0.0);
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,
                        boundXsVec,teff,iinc);
  Float xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001);
  Float muLeft = muVec[2], xsLeft = xsVec[2];
  // The outer loop will check between i-2 and i-1 to see if we need a midpoint
  // in between there. if so, it'll put it in there and then check between 
  // midpoint and i-1. If we don't need a midpoint, we go to the inner loop
  // and move to the left and keep checking. 
  int i = 3;
  Float muRight = muVec[i-1];
  do {
    addMidpointsRight(i,muVec,xsVec,e,ep,tev,alphas,betas,sab,az,lasym,lat,
                      boundXsVec,teff,iinc,tol,xsMax);
    do { // If i = 2, then do this action twice. Else do it once
      pdf += 0.5*( xsVec[i-1] + xsLeft )*( muRight - muLeft );
      muLeft = muVec[i-1];
      xsLeft = xsVec[i-1];
      --i;
      if (i > 0){ muRight = muVec[i-1]; }
    } while ( i == 1 );
  } while ( i > 1 );
  return pdf;
}






template <typename Range, typename Float>
inline auto sigl(Float ep, Float e, Float tev, Float tol, int lat, int iinc, 
  const Range& alphas, const Range& betas, const Range& sab, Float az,
  int lasym, const Range& boundXsVec, Float teff, Range& s, 
  bool equiprobableBins = true ){

  int nbin = s.size();
  for ( auto& sVal : s ){ sVal = 0.0; }
  tol *= 0.5;

  // Integrate the incoherent inelastic scattering xs over all mu values, so 
  // that we can later divide by the requested # of angle bins so that we know
  // how much to fill up each bin
  Float pdf = getPDF(ep, e, tev, tol, lat, iinc, alphas, betas, sab, az, lasym, 
                     boundXsVec, teff);


  //Range s(nbin,0.0);
  if ( pdf <= 1e-32 ){ return pdf; }

  Range muVec(20), xsVec(20); 
  initialize_XS_MU_vecs(muVec,xsVec,e,ep,az,tev,alphas,betas,sab,lasym,lat,
                        boundXsVec,teff,iinc); // 130

  Float xsMax  = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001),
        muLeft = muVec[2],
        muRight,
        xsLeft = xsVec[2],
        gral   = 0.0,
        sum    = 0.0,
        xs_per_bin  = pdf/(1.0*nbin); // about how much integrated xs we want per bin 

  int i = 3, j = 0;

  do { 
    addMidpointsRight(i,muVec,xsVec,e,ep,tev,alphas,betas,sab,az,lasym,lat,
                      boundXsVec,teff,iinc,tol,xsMax);
    do { //std::cout << " --- 160 ---  "<< std::endl;
      muRight = muVec[i-1];
      if (muRight == muLeft) {  
        // Is this our first time in the outer do loop? i.e. are we trying to
        // compare -1 with -1? If so, let's just move i down one and try again.
        --i; //std::cout << " --- 250 ---  "<< std::endl;
        muRight = muVec[i-1];
        if ( i < 1 ){ return pdf; }
        if ( i > 1 ){ break; }
      }

      Float invDeltaMu = 1.0/(muRight-muLeft);

      // This is the integrated little block of xs (from muLeft -> muVec[i-1]).
      // We're going to want to see if adding this to the sum will make us reach
      // our desired xs from xs_per_bin (recall that xs_per_bin is the xs we will have 
      // in each bin if they're perfectly equiprobable).
      Float add = 0.5*(xsVec[i-1]+xsLeft)*(muRight-muLeft);

      if ( i == 1 and j == nbin-1 ){ //std::cout << " --- 165 ---  "<< std::endl;
        j++;
        Float mu = muRight;
        populateXSvector(xsVec, mu, invDeltaMu, muLeft, xsLeft, xs_per_bin, i, 
                         gral, nbin, s, j, equiprobableBins);
        return pdf;
      }
      else if ( sum + add >= xs_per_bin * 0.99999999 and j < nbin-1 ){
        // Here, we're trying to get sum+add as close to possible as xs_per_bin. 
        // Fract is equalt to the cross section (integrated over all mu) divided
        // by the number of bins. So we're going to keep on adding to this sum
        // tally and giving it more and more bin bits until it looks like its 
        // about to overflow.
        j++;
        Float mu = getNextMuValue(xs_per_bin, sum, xsVec, xsLeft, muLeft, 
                                  muRight, i, invDeltaMu );
        populateXSvector(xsVec, mu, invDeltaMu, muLeft, xsLeft, xs_per_bin, i, 
                         gral, nbin, s, j, equiprobableBins);
        sum = 0.0;

        if (muLeft < muRight ){ continue; }
        else {  //std::cout << " --- 250 ---  "<< std::endl;
          muLeft = muRight;
          xsLeft = xsVec[i-1];
          --i;
          muRight = muVec[i-1];
          if ( i < 1 ){ return pdf; }
          if ( i > 1 ){ break; }
        }
      }

      sum  += add;

      gral += 0.5 * (xsLeft*muRight-xsVec[i-1]*muLeft) * 
                    (muRight+muLeft) + 0.333333333*(xsVec[i-1]-xsLeft) * 
                    (muRight*muRight+muRight*muLeft+muLeft*muLeft);

      muLeft = muRight;
      xsLeft = xsVec[i-1];
      --i;
      muRight = muVec[i-1];

    } while ( i == 1 );

  } while ( i > 1 );

  return pdf;
}











