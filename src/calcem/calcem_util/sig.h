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
  //std::cout << exp(-(alpha*alpha+beta*beta)/(4.0*alpha)) << std::endl;
  //std::cout << exp(-(alpha*alpha+beta*beta)/(4.0*alpha))*exp(beta*0.5)<< std::endl;
  //std::cout << pow(4.0*M_PI*alpha,-0.5) << std::endl;
  return cutoff(sab_to_xs_consts*exp(-beta*0.5)*sab);
} 


template <typename Float, typename Range>
inline auto interpolateSAB( int a, int b, const Range& alphas, 
  const Range& betas, const Range& sab, Float alpha, Float abs_beta_scaled, 
  Float beta ){
  
  a = ( a+1 == int(alphas.size()) ) ? a-1 : a;
  b = ( b+1 == int( betas.size()) ) ? b-1 : b;

  int n = betas.size();
  std::vector<double> sVec(3);
  for ( int i = -1; i < 2; ++i ){
      sVec[i+1] = terpq( alphas[a-1], alphas[a], alphas[a+1], alpha,
                  get(sab,a-1,b+i,n), get(sab,a,b+i,n), get(sab,a+1,b+i,n) );
  }
  Float s  = terpq( betas[b-1], betas[b], betas[b+1], abs_beta_scaled, sVec[0], sVec[1], sVec[2] );
  return exp(s-beta*0.5);
}


template <typename Float>
inline auto doSCT( Float a, Float teff, Float sigc, 
  Float sigma_b2, Float tev, Float sabflg, Float bb, Float sigma_b ){
 /* The SCT approximation is calculated, according to Eq. 230. Note that since
  * some evaluations give S(a,b) for a molecule or compund (e.g. C6H6 or BeO), 
  * the corresponding SCT approximation must contin terms for both atoms.
  * So in this case (when s2 > 0) we need to recalculate Eq. 230 taking into 
  * account the secondary atom's parameters.
  */
  using std::abs;
  Float arg, sigVal = 0, s;
  // If sigma_b2 > 0, then we're going to have to add the second atom's contribution to sig. 
  for ( const Float& sigma_b_val : {sigma_b, sigma_b2} ){ // Eq. 230 in the NJOY manual
    arg = (a-abs(bb))*(a-abs(bb))*tev/(4.0*a*teff) + (abs(bb)+bb)/2.0;
    s   = (-arg > sabflg) ? exp(-arg)/(sqrt(4.0*M_PI*a*teff/tev)) : 0;
    sigVal += sigc * sigma_b_val * s;
  }
  return sigVal;
}


template <typename Float, typename Range>
inline auto sig( const Float& e, const Float& ep, const Float& u, 
  const Float& tev, const Range& alphas, const Range& betas, const Range& sab, 
  const Float az, const Float tevz, const int lasym, const int lat, 
  const Float cliq, const Float sigma_b, const Float sigma_b2, const Float teff, 
  const int iinc ){

 /*-------------------------------------------------------------------
  * Compute the differential scattering cross section from e to
  * ep through the angle with the cosine u from endf tabulated
  * data or an analytic law.
  *-------------------------------------------------------------------
  */
  using std::abs; using std::distance; using std::lower_bound; using std::pow;

  int i,ib,ia;
  Float bb,sigc,s,s1,s2,s3,bbm,sctResult;
  Float sabflg=-225.e0, amin=1.e-6, test1=0.2e0, test2=30.e0;
  Float sigVal;

  // common factors.
  sigc = sqrt(ep/e) / (2.0*tev);

  if (iinc != 1 and iinc != 2){ throw std::exception(); }


  Float alpha = (e+ep-2*u*sqrt(e*ep))/(az*tev);
  Float beta  = (ep-e)/tev;
  if ( alpha < amin ){ alpha = amin; }

  Float maxAlpha = alphas[alphas.size()-1],
        maxBeta  = betas [betas.size() -1];
 
  // ----------------------------------------------------------------------
  // ------------------------ FREE GAS SCATTERING ------------------------- 
  // ----------------------------------------------------------------------
  if ( iinc == 1 ){  return freeGas( alpha, beta, sigma_b*sigc ); }



  // ----------------------------------------------------------------------
  // -------------------------- BOUND SCATTERING -------------------------- 
  // ----------------------------------------------------------------------
  Float alpha_scaled = ( lat == 1 ) ? alpha*tev/tevz : alpha;
  Float  beta_scaled = ( lat == 1 ) ? beta*tev/tevz  : beta;


  // ........................................................................//
  // Check to see if requested alpha or beta values are outside the grids    //
  // provided. If so, just do the SCT approximation.                         //
  // ........................................................................//
  
  bool beta_out_of_bounds = (lasym == 1 and (beta_scaled<betas[0] or beta_scaled >maxBeta)) 
                         or (lasym != 1 and                      abs(beta_scaled)>maxBeta);
  if ( alpha_scaled > maxAlpha or beta_out_of_bounds) {  // go to 170
    return cutoff( doSCT( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
  }


  // ........................................................................//
  // Not really sure why this is here or who it's helping. Sure isn't me.    //
  // ........................................................................//
  if (cliq != 0.0 and alpha_scaled < alphas[0] and lasym != 1 and abs(beta_scaled) < test1){
    s = sab[0] + log(alphas[0]/alpha_scaled)*0.5 - cliq*pow(beta_scaled,2)/alpha_scaled;
    return cutoff( sigc * sigma_b * exp(s-beta*0.5) );
  }


  // ........................................................................//
  // Not really sure why this is here or who it's helping. Sure isn't me.    //
  // ........................................................................//
   int ia1 = distance(alphas.cbegin(),lower_bound(alphas.cbegin(),alphas.cend(),  alpha_scaled)),
      ib1 = distance(betas.cbegin(), lower_bound(betas.cbegin(), betas.cend(),abs(beta_scaled))),
      nbeta = betas.size();
  ia = std::max(ia1,1); ib = std::max(ib1,1);
  if ( ( alpha_scaled*az >= test2 or abs(beta_scaled) >= test2 ) ){
    for (int a : {ia-1,ia}){
      for (int b : {ib-1,ib}){
        if ( get(sab,a,b,nbeta) <= sabflg ){
          return cutoff( doSCT(alpha,teff,sigc,sigma_b2,tev,sabflg,beta,sigma_b) );
        }
      }
    }
    return cutoff( doSCT(alpha,teff,sigc,sigma_b2,tev,sabflg,beta,sigma_b) );
  }


  auto sabVal = interpolateSAB( ia, ib, alphas, betas, sab, alpha_scaled, abs(beta_scaled), beta);
  return cutoff( sigc * sigma_b * sabVal);


}


#endif
