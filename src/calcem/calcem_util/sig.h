#ifndef THERMR_SIG
#define THERMR_SIG

#include "calcem/calcem_util/sig_util/terpq.h"
#include "general_util/sigfig.h"
#include <cmath>

template <typename Float>
auto cutoff( const Float proposedAnswer, const Float& cutoff=1e-10 ){
    return (proposedAnswer < cutoff) ? 0.0 : proposedAnswer;
}

template <typename Float>
inline auto do160(Float sigc, Float sigma_b, Float s, Float bb, Float sabflg ){
  // Based off of Eq. 225, it appears that s = S(a,b). 
  // sigc = sqrt(E'/E)/(2 kb T), sigma_b = sigma_b, and bb = beta. 
  // Note that calling do160 at the very and of sig seems mysterious, because
  // S(a,b) = S(a0,b0) + 1/2 ln(a0/a) - a0/a * b^2/b0^2 ( S(a0,b0) - S(a0,b1) )
  // and I really don't know where this comes from. (note the last term is b/c
  // of variable cliq that is calculated elsewhere and fed into sig)
  if (s-bb*0.5 <= sabflg ){ return 0.0; }
  return sigc * sigma_b * exp(s-bb/2);
}

template <typename Float, typename Range>
inline auto do155( int ia, int ib, const Range& alpha,
    const Range& beta, const std::vector<std::vector<Float>>& sab,
    Float a, Float bbb, Float sigc, Float sigma_b, Float bb, Float sabflg ){
   if ( (unsigned) ia+1 == alpha.size()) ia -= 1;
   if ( (unsigned) ib+1 == beta.size())  ib -= 1;
   ia -= 1; ib -= 1;
   Float s1 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+0], sab[ia+1][ib+0], sab[ia+2][ib+0] ),
          s2 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+1], sab[ia+1][ib+1], sab[ia+2][ib+1] ),
          s3 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+2], sab[ia+1][ib+2], sab[ia+2][ib+2] );
   Float s  = terpq( beta[ib], beta[ib+1], beta[ib+2], bbb, s1, s2, s3 );
   return do160(sigc, sigma_b, s, bb, sabflg );
}


template <typename Float>
inline auto doSCTApproximation( Float a, Float teff, Float sigc, 
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
  const Float& tev,  
  const Range& alpha, const Range& beta,
  const std::vector<std::vector<Float>>& sab, const Float az, 
  const Float tevz, const int lasym,// const Float az2, const Float teff2, 
  const int lat, const Float cliq, const Float sigma_b, const Float sigma_b2, 
  const Float teff, const int iinc ){

 /*-------------------------------------------------------------------
  * Compute the differential scattering cross section from e to
  * ep through the angle with the cosine u from endf tabulated
  * data or an analytic law.
  *-------------------------------------------------------------------
  */
  int i,ib,ia;
  Float bb,a,sigc,b,s,s1,s2,s3,bbm,sctResult;
  Float sabflg=-225.e0, amin=1.e-6, test1=0.2e0, test2=30.e0;
  Float sigVal;

  // common factors.
  sigc=sqrt(ep/e)/(2.0*tev);


  if (iinc != 1 and iinc != 2){ throw std::exception(); }
  // -------------------------------------------------------------------------
  // So the whole source of confusion here is that leapr has this option with
  // lat, where if you set lat = 1 it will scale all alpha, beta values by
  // 0.0253/tev, and if lat = 0 then no such operation will be performed. Here,
  // we're checking to see if lat = 1 or not and undoing this scaling. 
  // -------------------------------------------------------------------------
  Float a_tev   = (e+ep-2*u*sqrt(e*ep))/(az*tev);
  Float a_tevz  = (e+ep-2*u*sqrt(e*ep))/(az*tevz);
  Float bb_tev  = (ep-e)/tev;     
  Float bb_tevz = (ep-e)/tevz;

  if ( a_tev < amin ){ a_tev = amin; a_tevz = amin * tev / tevz; }

  if ( lat == 1 ){ b = std::abs(bb_tevz); a = a_tevz; bbm = bb_tevz; }
  else           { b = std::abs(bb_tev ); a = a_tev ; bbm = bb_tev ; }
  b = sigfig(b,8,0);


  // ----------------------------------------------------------------------
  // ------------------------ FREE GAS SCATTERING ------------------------ 
  // ----------------------------------------------------------------------
  // free-gas scattering. Plug Eq. 229 into Eq. 225, solve.
  if ( iinc == 1 ){ 
    // S(a,b) = 1/sqrt(4pi*a) * exp( - ( a^2 + b^2 ) / 4a )          [Eq. 229]
    Float sab = 1.0/sqrt(4.0*M_PI*a_tev) * exp( -(a_tev*a_tev + bb_tev*bb_tev)/(4.0*a_tev) );
    // sigma = sigma_b / 2 kb T  * sqrt(E'/E) * exp( -b/2 ) * S(a,b) [Eq. 225]
    // Note that sigma_b is the bound scattering cross section, which is 
    // defined by Eq. 228 to be 
    //             sigma_b = sigma_f * ( A + 1 )^2 / A^2             [Eq. 228]
    // where sigma_f is the characteristic free scattering cross section
    sigVal = sigma_b/(2.0*tev) * sqrt(ep/e) * exp(-bb_tev/2.0) * sab;
    return cutoff( sigVal );
  } // end free gas scattering option


  // ----------------------------------------------------------------------
  // -------------------------- BOUND SCATTERING -------------------------- 
  // ----------------------------------------------------------------------
  //bool outOfRange = false;

  if (a > alpha[alpha.size()-1]) {  // go to 170
    return cutoff( doSCTApproximation( a_tev, teff, sigc, sigma_b2, tev, sabflg, bb_tev, sigma_b ) );
  }

  Float maxBeta = beta[beta.size()-1];

  // lasym = 1 means Symmetric S(a,b) and ortho/para hydrogen (need +/- beta sides)
  if ( (lasym == 1 and (bbm > maxBeta or bbm < beta[0])) or 
                       (b   > maxBeta                 ) ){ 
      return cutoff( doSCTApproximation( a_tev, teff, sigc, sigma_b2, tev, sabflg, bb_tev, sigma_b ) );
  }

  if (cliq == 0.0 or a >= alpha[0] or lasym == 1 or b >= test1 ){ 
    // Find index ia such that alpha[ia-1] < a < alpha[ia]. Same for ib
    auto ia1 = std::distance(alpha.cbegin(), std::lower_bound(alpha.cbegin(), alpha.cend(), a));
    auto ib1 = std::distance(beta.cbegin(),  std::lower_bound(beta.cbegin(),  beta.cend(), b));
    ia = ia1 > 1 ? ia1 : 1;
    ib = ib1 > 1 ? ib1 : 1;

    // 150 continue
    if (a*az < test2 and b < test2) { 
      return cutoff( do155( ia, ib, alpha, beta, sab, a, b, sigc, sigma_b, bb_tev, sabflg ) );
    }
    if ( a*az >= test2 or b >= test2 ){
      if ( sab[ia-1][ib-1]   <= sabflg or sab[ia+1-1][ib-1]   <= sabflg or 
           sab[ia-1][ib+1-1] <= sabflg or sab[ia+1-1][ib+1-1] <= sabflg ) {
        return cutoff( doSCTApproximation( a_tev, teff, sigc, sigma_b2, tev, sabflg, bb_tev, sigma_b ) );
      } 
      else { 
        return cutoff( do155( ia, ib, alpha, beta, sab, a, b, sigc, sigma_b, bb_tev, sabflg ) );
      }
    }
  }


  s=sab[0][0]+log(alpha[0]/a)/2-cliq*b*b/a;
  if (s < sabflg) s=sabflg;
  return cutoff( do160(sigc,sigma_b,s,bb_tev,sabflg) );

}



#endif
