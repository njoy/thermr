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
Float freeGas( Float alpha, Float beta, Float sab_to_xs_consts ){
  // S(a,b) = 1/sqrt(4pi*a) * exp( - ( a^2 + b^2 ) / 4a )          [Eq. 229]
  Float sab = 1.0/sqrt(4.0*M_PI*alpha) * exp( -(alpha*alpha + beta*beta)/(4.0*alpha) );
  // sigma = sigma_b / 2 kb T  * sqrt(E'/E) * exp( -b/2 ) * S(a,b) [Eq. 225]
  // Note that sigma_b is the bound scattering cross section, which is 
  // defined by Eq. 228 to be 
  //             sigma_b = sigma_f * ( A + 1 )^2 / A^2             [Eq. 228]
  // where sigma_f is the characteristic free scattering cross section
  return cutoff( sab_to_xs_consts * exp(-beta/2.0) * sab );
} 




template <typename Float>
inline auto do160(Float sigc, Float sigma_b, Float s, Float bb, Float sabflg ){
  // Based off of Eq. 225, it appears that s = S(a,b). 
  // sigc = sqrt(E'/E)/(2 kb T), sigma_b = sigma_b, and bb = betas. 
  // Note that calling do160 at the very and of sig seems mysterious, because
  // S(a,b) = S(a0,b0) + 1/2 ln(a0/a) - a0/a * b^2/b0^2 ( S(a0,b0) - S(a0,b1) )
  // and I really don't know where this comes from. (note the last term is b/c
  // of variable cliq that is calculated elsewhere and fed into sig)
  if (s-bb*0.5 <= sabflg ){ return 0.0; }
  return sigc * sigma_b * exp(s-bb/2);
}

template <typename Float, typename Range>
inline auto do155( int ia, int ib, const Range& alphas,
    const Range& betas, const std::vector<std::vector<Float>>& sab,
    Float a, Float bbb, Float sigc, Float sigma_b, Float bb, Float sabflg ){
   if ( (unsigned) ia+1 == alphas.size()) ia -= 1;
   if ( (unsigned) ib+1 == betas.size())  ib -= 1;
   ia -= 1; ib -= 1;
   Float s1 = terpq( alphas[ia], alphas[ia+1], alphas[ia+2], a, sab[ia][ib+0], sab[ia+1][ib+0], sab[ia+2][ib+0] ),
          s2 = terpq( alphas[ia], alphas[ia+1], alphas[ia+2], a, sab[ia][ib+1], sab[ia+1][ib+1], sab[ia+2][ib+1] ),
          s3 = terpq( alphas[ia], alphas[ia+1], alphas[ia+2], a, sab[ia][ib+2], sab[ia+1][ib+2], sab[ia+2][ib+2] );
   Float s  = terpq( betas[ib], betas[ib+1], betas[ib+2], bbb, s1, s2, s3 );
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
  const Range& alphas, const Range& betas,
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
  using std::abs; using std::distance; using std::lower_bound;

  int i,ib,ia;
  Float bb,sigc,s,s1,s2,s3,bbm,sctResult;
  Float sabflg=-225.e0, amin=1.e-6, test1=0.2e0, test2=30.e0;
  Float sigVal;

  // common factors.
  sigc=sqrt(ep/e)/(2.0*tev);


  if (iinc != 1 and iinc != 2){ throw std::exception(); }
  // -------------------------------------------------------------------------
  // So the whole source of confusion here is that leapr has this option with
  // lat, where if you set lat = 1 it will scale all alphas, betas values by
  // 0.0253/tev, and if lat = 0 then no such operation will be performed. Here,
  // we're checking to see if lat = 1 or not and undoing this scaling. 
  // -------------------------------------------------------------------------

  Float alpha = (e+ep-2*u*sqrt(e*ep))/(az*tev);
  Float beta  = (ep-e)/tev;
  if ( alpha < amin ){ alpha = amin; }

  Float maxAlpha = alphas[alphas.size()-1],
        maxBeta  = betas [betas.size() -1];
 
  // ----------------------------------------------------------------------
  // ------------------------ FREE GAS SCATTERING ------------------------ 
  // ----------------------------------------------------------------------
  if ( iinc == 1 ){ 
    return freeGas( alpha, beta, sigma_b*sigc );
  }

  // ----------------------------------------------------------------------
  // -------------------------- BOUND SCATTERING -------------------------- 
  // ----------------------------------------------------------------------
  Float alpha_scaled = alpha;
  Float beta_scaled = beta;
  if ( lat == 1 ){ alpha_scaled *= tev/tevz; }
  if ( lat == 1 ){ beta_scaled  *= tev/tevz; }

  bool beta_out_of_bounds = (lasym == 1 and (beta_scaled<betas[0] or beta_scaled>maxBeta)) 
                         or (lasym != 1 and              abs(beta_scaled)>maxBeta) ;
  if ( alpha_scaled > maxAlpha or beta_out_of_bounds) {  // go to 170
    return cutoff( doSCTApproximation( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
  }

  if (cliq != 0.0 and alpha_scaled < alphas[0] and lasym != 1 and abs(beta_scaled) < test1){
    s=sab[0][0]+log(alphas[0]/alpha_scaled)/2-cliq*(beta_scaled*beta_scaled)/alpha_scaled;
    return cutoff( sigc * sigma_b * exp(s-beta*0.5) );
  }

  int ia1 = distance(alphas.cbegin(),lower_bound(alphas.cbegin(),alphas.cend(),  alpha_scaled));
  int ib1 = distance(betas.cbegin(), lower_bound(betas.cbegin(), betas.cend(),abs(beta_scaled)));
  ia = std::max(ia1,1); ib = std::max(ib1,1);

  if ( ( alpha_scaled*az >= test2 or abs(beta_scaled) >= test2 ) and 
       ( sab[ia-1][ib-1]   <= sabflg or sab[ia+1-1][ib-1]   <= sabflg or 
         sab[ia-1][ib+1-1] <= sabflg or sab[ia+1-1][ib+1-1] <= sabflg ) ) {
    return cutoff( doSCTApproximation( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
  }
  else {
    return cutoff( do155( ia, ib, alphas, betas, sab, alpha_scaled, abs(beta_scaled), sigc, sigma_b, beta, sabflg ) );
  }


}


template <typename Float, typename Range>
inline auto sigDIFFERENT( const Float& e, const Float& ep, const Float& u, 
  const Float& tev,  
  const Range& alphas, const Range& betas,
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
  // lat, where if you set lat = 1 it will scale all alphas, betas values by
  // 0.0253/tev, and if lat = 0 then no such operation will be performed. Here,
  // we're checking to see if lat = 1 or not and undoing this scaling. 
  // -------------------------------------------------------------------------
  Float a_tev  = (e+ep-2*u*sqrt(e*ep))/(az*tev);
  Float b_tev  = (ep-e)/tev;     
  Float a_tevz = (e+ep-2*u*sqrt(e*ep))/(az*tevz);
  Float b_tevz = (ep-e)/tevz;

  Float alpha = (e+ep-2*u*sqrt(e*ep))/(az*tev);
  Float beta  = (ep-e)/tev;
  if ( alpha < amin ){ alpha = amin; }

  if ( a_tev < amin ){ a_tev = amin; a_tevz = amin * tev / tevz; }
  
  
  // ----------------------------------------------------------------------
  // ------------------------ FREE GAS SCATTERING ------------------------ 
  // ----------------------------------------------------------------------
  if ( iinc == 1 ){ 
    return freeGas( alpha, beta, sigma_b*sigc );
  }

  // ----------------------------------------------------------------------
  // -------------------------- BOUND SCATTERING -------------------------- 
  // ----------------------------------------------------------------------
  // IINC must equal 2 here
 
  if ( lat == 1 ){ b = std::abs(b_tevz); a = a_tevz; bbm = b_tevz; }
  else           { b = std::abs(b_tev ); a = a_tev ; bbm = b_tev ; }
  b = sigfig(b,8,0);

  std::vector<double> alphas2 (alphas.size(),0.0);
  for (size_t i = 0; i < alphas.size(); ++i){
    alphas2[i] = alphas[i];
    if (lat == 1){ alphas2[i] *= tevz/tev; }
  }
  std::vector<double> betas2 (betas.size(),0.0);
  for (size_t i = 0; i < betas.size(); ++i){
    betas2[i] = betas[i];
    if (lat == 1){ betas2[i] *= tevz/tev; }
  }

  if ( lat == 1 ){ 
      test1 *= tevz/tev; 
      test2 *= tevz/tev; 
  }

  Float maxBeta = betas[betas.size()-1];
  Float maxBeta2 = betas2[betas2.size()-1];
  Float minAlpha = alphas2[0];
  Float maxAlpha= alphas2[alphas2.size()-1];

  if (alpha > maxAlpha) {  // go to 170
    return cutoff( doSCTApproximation( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
  }



  // lasym = 1 means Symmetric S(a,b) and ortho/para hydrogen (need +/- betas sides)
  //if ( (lasym == 1 and (         bbm  > maxBeta or bbm < betas[0])) or 
  //                     (std::abs(bbm) > maxBeta                 ) ){ 

  if ( (lasym == 1 and (         beta  > maxBeta2 or beta < betas[0])) or 
                       (std::abs(beta) > maxBeta2                 ) ){ 
    return cutoff( doSCTApproximation( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
  }

  if (cliq == 0.0 or alpha >= minAlpha or lasym == 1 or beta >= test1 ){ 
    // Find index ia such that alphas[ia-1] < a < alphas[ia]. Same for ib
    auto ia1 = std::distance(alphas.cbegin(), std::lower_bound(alphas.cbegin(), alphas.cend(), a));
    auto ib1 = std::distance(betas.cbegin(),  std::lower_bound(betas.cbegin(),  betas.cend(), b));


    //std::cout << (betas|ranges::view::all) << std::endl;
    //std::cout << (betas2|ranges::view::all) << std::endl;
    //std::cout << std::endl;

    auto ia2 = std::distance(alphas2.cbegin(), std::lower_bound(alphas2.cbegin(), alphas2.cend(), alpha));
    auto ib2 = std::distance(betas2.cbegin(), std::lower_bound(betas2.cbegin(), betas2.cend(), beta));



    //std::cout << lat << "      " <<ib1 << "    " << ib2 << std::endl;
    //std::cout << lat << "      " <<ib1 << "    " << betas[ib1-1] << "   " << b << "   " << betas[ib1] << std::endl;
    //std::cout << lat << "      " <<ib2 << "    " << betas2[ib2-1] << "   " << beta << "   " << betas2[ib2] << std::endl;
    //std::cout << lat << "  " << ib2 << "    " << betas2[ib2-1] << "   " << beta << "   " << betas2[ib2] << std::endl;
    ia = ia1 > 1 ? ia1 : 1;
    ib = ib1 > 1 ? ib1 : 1;
    auto ia_2 = ia2 > 1 ? ia2 : 1;
    auto ib_2 = ib2 > 1 ? ib2 : 1;


    if ( (alpha*az >= test2 or beta >= test2) and
         (sab[ia-1][ib-1] <= sabflg or sab[ia][ib-1] <= sabflg or 
          sab[ia-1][ib  ] <= sabflg or sab[ia][ib]   <= sabflg ) ) {
      return cutoff( doSCTApproximation( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
    }
    else { 
      return cutoff( do155( ia_2, ib_2, alphas2, betas2, sab, alpha, std::abs(beta), sigc, sigma_b, b_tev ) );
    }



    /*
    if (a*az < test2 and b < test2) { 
      //return cutoff( do155( ia, ib, alphas, betas, sab, a, b, sigc, sigma_b, b_tev, sabflg ) );
      return cutoff( do155( ia_2, ib_2, alphas2, betas2, sab, alpha, std::abs(beta), sigc, sigma_b, b_tev, sabflg ) );
    }
    else {
    if ( a*az >= test2 or b >= test2 ){
      if ( sab[ia-1][ib-1] <= sabflg or sab[ia][ib-1] <= sabflg or 
           sab[ia-1][ib  ] <= sabflg or sab[ia][ib]   <= sabflg ) {
        return cutoff( doSCTApproximation( alpha, teff, sigc, sigma_b2, tev, sabflg, beta, sigma_b ) );
      } 
      else { 
        return cutoff( do155( ia_2, ib_2, alphas2, betas2, sab, alpha, std::abs(beta), sigc, sigma_b, b_tev, sabflg ) );
      }
    }
    */
  }


  s = std::max(sab[0][0]+log(alphas[0]/a)/2-cliq*b*b/a , sabflg);
  if (s-b_tev*0.5 <= sabflg ){ return 0.0; }
  return cutoff( sigc * sigma_b * exp(s-b_tev/2));


}



#endif
