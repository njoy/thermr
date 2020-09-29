#ifndef THERMR_SIG
#define THERMR_SIG

#include "calcem/calcem_util/sig_util/terpq.h"
#include "general_util/sigfig.h"
#include <cmath>

template <typename Range>
auto get(const Range& sab, int a, int b, int betas_size){
  return sab[a*betas_size+b];
}

template <typename Float>
auto cutoff( const Float proposedAnswer, const Float& cutoff=1e-10 ){
    return (proposedAnswer < cutoff) ? 0.0 : proposedAnswer;
}

template <typename Float>
Float freeGas( Float alpha, Float beta, Float sab_to_xs_consts ){ 
  // Follows Eq. 225, 228, 229 in LEAPR manual.
  Float sab = pow(4.0*M_PI*alpha,-0.5) * exp(-(alpha*alpha+beta*beta)/(4.0*alpha));
  return cutoff(sab_to_xs_consts*exp(-beta*0.5)*sab);
} 


template <typename Float, typename Range>
auto getIndices( const Range& vector, const Float& value ){
  for ( size_t i = 0; i < vector.size()-1; ++i ){
    if (value < vector[i+1]){ return i+1; }
  }
  return vector.size()-1;
//  int ia = distance(alphas.cbegin(),lower_bound(alphas.cbegin(),alphas.cend(),  alpha)),
//      ib = distance(betas.cbegin(), lower_bound(betas.cbegin(), betas.cend(),abs(beta)));
//  ia = std::max(ia,1); ib = std::max(ib,1);
//  return std::make_tuple(ia,ib);
}

template <typename Float, typename Range>
inline auto interpolateSAB( int a, int b, const Range& alphas, 
  const Range& betas, const Range& sab, Float alpha, Float absBetaScaled, 
  Float beta ){

  a = ( a+1 == int(alphas.size()) ) ? a-1 : a;
  b = ( b+1 == int( betas.size()) ) ? b-1 : b;

  int n = betas.size();
  Range sVec(3);
  for ( int i = -1; i < 2; ++i ){
    sVec[i+1] = terpq( alphas[a-1], alphas[a], alphas[a+1], alpha,
                get(sab,a-1,b+i,n), get(sab,a,b+i,n), get(sab,a+1,b+i,n) );
  }
  Float s  = terpq( betas[b-1], betas[b], betas[b+1], absBetaScaled, 
                     sVec[0], sVec[1], sVec[2] );
  return exp(s-beta*0.5);
}

template <typename Float>
inline auto doSCT( Float alpha, Float teff, Float sigc, 
  Float sigma_b2, Float tev, Float sabflg, Float beta, Float sigma_b ){
 /* The SCT approximation is calculated, according to Eq. 230. Note that since
  * some evaluations give S(a,b) for a molecule or compund (e.g. C6H6 or BeO), 
  * the corresponding SCT approximation must contin terms for both atoms.
  * So in this case (when s2 > 0) we need to recalculate Eq. 230 taking into 
  * account the secondary atom's parameters.
  */
  using std::abs; using std::pow;
  Float arg, sigVal = 0, s;
  // If sigma_b2 > 0, then we add the second atom's contribution to sig. 
  for ( const Float& sigma_b_val : {sigma_b, sigma_b2} ){ // Eq. 230 
    arg = pow((alpha-abs(beta)),2)*tev/(4.0*alpha*teff) + (abs(beta)+beta)/2.0;
    s   = (-arg > sabflg) ? exp(-arg)/(sqrt(4.0*M_PI*alpha*teff/tev)) : 0;
    sigVal += sigc * sigma_b_val * s;
  }
  return sigVal;
}

template <typename Float, typename Range>
inline auto sig( const Float& e, const Float& ep, const Float& u, 
  const Float& tev, const Range& alphas, const Range& betas, const Range& sab, 
  const Float az, const Float tevz, const int lasym, const int lat, 
  const Range& boundXsVec, const Float teff, 
  const int iinc ){
 /*-------------------------------------------------------------------
  * Compute the differential scattering cross section from e to
  * ep through the angle with the cosine u from endf tabulated
  * data or an analytic law.
  *-------------------------------------------------------------------
  */
  using std::abs; using std::distance; using std::pow; using std::max;

  Float sabflg=-225.e0, amin=1.e-6, test1=0.2e0, test2=30.e0,
        sigc = sqrt(ep/e)/(2.0*tev);

  if (not (iinc == 1 or iinc == 2) ){ throw std::exception(); }

  Float alpha = max((e+ep-2*u*sqrt(e*ep))/(az*tev),amin),
        beta  = (ep-e)/tev;

  Float maxAlpha = alphas[alphas.size()-1],
        maxBeta  = betas [betas.size() -1];

  Float sigma_b  = boundXsVec[0];
  Float sigma_b2 = boundXsVec[1];
 
  // ----------------------------------------------------------------------
  // ------------------------ FREE GAS SCATTERING ------------------------- 
  // ----------------------------------------------------------------------
  if ( iinc == 1 ){ 
    return freeGas( alpha, beta, sigma_b*sigc ); 
  }

  // ----------------------------------------------------------------------
  // -------------------------- BOUND SCATTERING -------------------------- 
  // ----------------------------------------------------------------------
  Float alphaScaled = ( lat == 1 ) ? alpha*tev/tevz : alpha;
  Float  betaScaled = ( lat == 1 ) ? beta*tev/tevz  : beta;

  // ........................................................................//
  // Check to see if requested alpha or beta values are outside the grids    //
  // provided. If so, just do the SCT approximation.                         //
  // ........................................................................//
  bool betaOutOfGrid = (lasym==1 and (betaScaled<betas[0] or betaScaled >maxBeta)) 
                    or (lasym!=1 and                     abs(betaScaled)>maxBeta);
  if (alphaScaled > maxAlpha or betaOutOfGrid) {  // go to 170
    return cutoff(doSCT(alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b));
  }

  // ........................................................................//
  //  Diffusion in liquids can create a sort of singularity when alpha is    //
  //  small and beta = 0. If it looks like we have a singularity, then we    //
  //  extrapolate using a beta^2/alpha law.
  // ........................................................................//
  Float cliq = (sab[0] - sab[1])*alphas[0]/(betas[1]*betas[1]);
  if (cliq != 0.0 and alphaScaled < alphas[0] and lasym != 1 and abs(betaScaled) < test1){
    Float s = sab[0] + log(alphas[0]/alphaScaled)*0.5 - cliq*pow(betaScaled,2)/alphaScaled;
    return cutoff( sigc * sigma_b * exp(s-beta*0.5) );
  }

  // ........................................................................//
  // Find where in the alphas, betas vectors our desired values are at.      // 
  // ........................................................................//
  Float absBeta = (lasym == 1 and beta < 0) ? -abs(betaScaled) : abs(betaScaled);
  int ia = getIndices(alphas,alphaScaled);
  int ib = getIndices(betas, absBeta);

  // ........................................................................//
  // Not really sure why this is here or who it's helping. Sure isn't me.    //
  // ........................................................................//
  if ( ( alphaScaled*az >= test2 or abs(betaScaled) >= test2 ) ){
    for (int a : {ia-1,ia}){
      for (int b : {ib-1,ib}){
        if ( get(sab,a,b,betas.size()) <= sabflg ){
          return cutoff( doSCT(alpha,teff,sigc,sigma_b2,tev,sabflg,beta,sigma_b) );
        }
      }
    }
    //return cutoff( doSCT(alpha,teff,sigc,sigma_b2,tev,sabflg,beta,sigma_b) );
  }


  // ........................................................................//
  // Just interpolate on the S(a,b) grid for the right value.                //
  // ........................................................................//

  auto sabVal = interpolateSAB( ia, ib, alphas, betas, sab, alphaScaled, abs(betaScaled), beta);
  return cutoff( sigc * sigma_b * sabVal);


}


#endif
