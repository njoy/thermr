#include <iostream>
#include <vector>
#include "terpq.h"


auto do170( int lat, double a, double b, double tevz, double tev, double teff, 
  double teff2, double az, double az2, double sigc, double sb, double sb2, 
  double arg, double e, double sigmin, double u, double c, double bb, 
  int sabflg ){

  // short collision time for large beta or alpha
  // 170 continue
  
  double sigVal, tfff, tfff2, s, s2 = 0;
  if (lat == 1) { 
    a = a * tevz / tev; 
    b = b * tevz / tev; 
  }
  tfff  = teff;
  tfff2 = teff2;

  arg = (a-b) * (a-b) * tev / (4*a*tfff) + (b+bb) / 2;
  s = -arg > sabflg ? exp(-arg) / (c*sqrt(a*tfff/tev)) : 0;
  sigVal = sigc * sb * s;
  if (sb2 > 0) {
    a *= az / az2;
    arg = (a-b)*(a-b)*tev/(4*a*tfff2)+(b+bb)/2;
    s2 = -arg > sabflg ? exp(-arg)/(c*sqrt(a*tfff2/tev)) : 0;
    sigVal = sigVal + sigc * sb2 * s2;
  }
  return sigVal < sigmin ? 0 : sigVal;
}



auto do155( int ia, int ib, int nalpha, int nbeta, double s, double bb,
  double sigc, double sb, double sigmin, double sabflg, 
  std::vector<double> alpha, std::vector<double> beta, 
  std::vector<std::vector<double>> sab, double a, double bbb ) {
  //155 continue
  if (ia+2 == nalpha) { ia = ia - 1; }
  if (ib+2 == nbeta)  { ib = ib - 1; }

  double s1, s2, s3;

  s1 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib], 
    sab[ia+1][ib], sab[ia+2][ib] );
  s2 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+1], 
    sab[ia+1][ib+1], sab[ia+2][ib+1] );
  s3 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+2], 
    sab[ia+1][ib+2], sab[ia+2][ib+2] );

  s = terpq( beta[ib], beta[ib+1], beta[ib+2], bbb, s1, s2, s3 );

  double sigVal = s-bb/2 > sabflg ? exp(s-bb/2) : 0.0;
  return sigVal < sigmin ? 0 : sigc * sb * sigVal;

}


auto sig( double e, double ep, double u, double tev, double tevz,
  std::vector<double>& alpha, std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, double az, double az2, int lat, 
  int iinc, int lasym, double cliq, double sb, double sb2, double teff,
  double teff2 ){
  /*-------------------------------------------------------------------
   * Compute the differential scattering cross section from e to
   * ep through the angle with the cosine u from endf tabulated
   * data or an analytic law.
   *-------------------------------------------------------------------
   */
  int nb1, na1,i,ib,ia,nbeta = beta.size(), nalpha = alpha.size();
  double bb,a,sigc,b,c,bbb,s,s1,s2,s3,arg,tfff,tfff2,rat,bbm,sigmin=1.e-10, 
         sabflg=-225.e0,amin=1.e-6,test1=0.2e0,test2=30.e0,sigVal;

  // common factors.
  bb = (ep-e) / tev;
  a = (e+ep-2*u*sqrt(e*ep)) / (az*tev);
  if (a < amin) { a = amin; }
  sigc = sqrt(ep/e) / ( 2 * tev );
  b = std::abs(bb);
  c = sqrt(4*M_PI);

  // tabulated s(alpha,beta).
  if (iinc != 2) { 
    // free-gas scattering.
    // 200 continue
    if (iinc != 1) {
      std::cout << "illegal option in sig! oh no!" << std::endl;
      return 0.0;
    }

    s = 0;
    arg = (a+bb)*(a+bb)/(4*a);
    if (-arg > sabflg) s = exp(-arg)/(c*sqrt(a));
    sigVal = sigc * sb * s;
    return sigVal < sigmin ? 0.0 : sigVal;

  }

  if (lat == 1){
    b = b * tev / tevz;
    a = a * tev / tevz;
  }

  if (a > alpha[nalpha-1]) { 
    return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
      sb, sb2, arg, e, sigmin, u, c, bb, sabflg );
  }

  if (lasym == 1) {
    bbm = lat == 1 ? bb * tev / tevz : bb;
    if ( bbm > beta[nbeta-1] or bbm < beta[0] ) {
      return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
        sb, sb2, arg, e, sigmin, u, c, bb, sabflg );
    }
  }
  else {
    if (b > beta[nbeta-1]) {
      return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
        sb, sb2, arg, e, sigmin, u, c, bb, sabflg );
    }
  } 
  nb1 = nbeta - 1;
  na1 = nalpha - 1;

  bbb = lasym == 1 and bb < 0 ? -b : b;

  for ( size_t i = 0; i < nb1; ++i ){
    ib = i;
    if (bbb < beta[i+1]) break;
  }
  for ( size_t i = 0; i < na1; ++i ) {
    ia = i;
    if ( a < alpha[i+1] ){ break; }
  }

  if (cliq == 0 or a >= alpha[0] or lasym == 1 or b > test1 ) {
    // This is 150
    if (a*az < test2 and b < test2) { 
      return do155( ia, ib, nalpha, nbeta, s, bb, sigc, sb, sigmin, sabflg, 
        alpha, beta, sab, a, bbb );
    } 

    if ( sab[ia][ib]   <= sabflg or sab[ia+1][ib]   <= sabflg or 
         sab[ia][ib+1] <= sabflg or sab[ia+1][ib+1] <= sabflg ){
      return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
        sb, sb2, arg, e, sigmin, u, c, bb, sabflg );
    }

  }

  s = sab[0][0] + log(alpha[0]/a)/2 - cliq*b*b/a;

  if (s < sabflg) { s = sabflg; }

  // go to 160
  sigVal = s-bb/2 > sabflg ? exp(s-bb/2) : 0;
  return sigVal < sigmin ? 0.0 : sigc * sb * sigVal;


  // free-gas scattering.
  // 200 continue
  // if (iinc != 1) go to 300
  arg = (a+bb)*(a+bb)/(4*a);
  s = (-arg > sabflg) ? exp(-arg)/(c*sqrt(a)) : 0;
  return sigVal < sigmin ? 0.0 : sigc * sb * s;

}


