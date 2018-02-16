#include <iostream>
#include <vector>
#include "terpq.h"

auto do170( int lat, double a, double b, double tevz, double tev, 
    double teff, double teff2, double az, double az2, 
    double sigc, double sb, double sb2, double arg, double e, 
    double tfff, double tfff2, double sigmin, double u, 
    double c, double bb, int sabflg ){
   // short collision time for large beta or alpha
   // 170 continue
   double sigVal;
   if (lat == 1) b = b * tevz / tev;
   if (lat == 1) a = a * tevz / tev;
   tfff = teff;
   tfff2 = teff2;

   double s = 0;
   arg = (a-b)*(a-b)*tev/(4*a*tfff)+(b+bb)/2;
   if (-arg > sabflg) s = exp(-arg) / (c*sqrt(a*tfff/tev));
   sigVal = sigc * sb * s;
   if (sb2 > 0) {
      a = a * az / az2;
      arg = (a-b)*(a-b)*tev/(4*a*tfff2)+(b+bb)/2;
      double s2 = 0;
      if (-arg > sabflg) s2 = exp(-arg)/(c*sqrt(a*tfff2/tev));
      sigVal = sigVal + sigc * sb2 * s2;
    }
    if (std::abs(e-10) < .01 and std::abs(u-.99219) < .0001) {
    }
   if (sigVal < sigmin) sigVal = 0;
   std::cout << sigVal << std::endl;
   return sigVal; 
}



auto do155( int ia, int ib, int nalpha, int nbeta, double s, double bb,
  double sigc, double sb, double sigmin, double sabflg, 
  std::vector<double> alpha, std::vector<double> beta, 
  std::vector<std::vector<double>> sab, double a, double bbb ) {
  //155 continue
  if (ia+2 == nalpha) ia = ia - 1;
  if (ib+2 == nbeta) ib = ib - 1;

  double s1, s2, s3;
  s1 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib], sab[ia+1][ib], sab[ia+2][ib] );
  std::cout << " ------------------- " << std::endl;
  std::cout << s1 << std::endl;
  s2 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+1], sab[ia+1][ib+1], sab[ia+2][ib+1] );
  std::cout << s2 << std::endl;
  s3 = terpq( alpha[ia], alpha[ia+1], alpha[ia+2], a, sab[ia][ib+2], sab[ia+1][ib+2], sab[ia+2][ib+2] );
  std::cout << s3 << std::endl;
  s = terpq( beta[ib], beta[ib+1], beta[ib+2], bbb, s1, s2, s3 );
  std::cout << s << std::endl;
  std::cout << bb << std::endl;
  std::cout << " ------------------- " << std::endl;

  std::cout << "160" << std::endl;
  double sigVal = 0;
  if (s-bb/2 > sabflg) sigVal = exp(s-bb/2);
  sigVal = sigc * sb * sigVal;
  //std::cout << "sigVal  " << sigVal << std::endl;
  if (sigVal < sigmin) sigVal = 0;
  //std::cout << "sigVal  " << sigVal << std::endl;
  return sigVal;

}





auto do150( double a, double az, double test2, double b, int ia, int ib,
  int sabflg, int lat, double tevz, double tev, 
  std::vector<std::vector<double>> sab, double teff, double teff2,
  double az2, double sigc, double sb, double sb2, double e, double tfff, 
  double tfff2, double arg, double sigmin, double u, double c, double bb, 
  int nalpha, int nbeta, double s, std::vector<double> alpha, 
  std::vector<double> beta, double bbb ){  //150 continue

  if (a*az < test2 and b < test2) { 
    std::cout << "155" << std::endl; 
    return do155( ia, ib, nalpha, nbeta, s, bb, sigc, sb, sigmin, sabflg, 
        alpha, beta, sab, a, bbb );
  } // go to 155

  if ( sab[ia][ib]   <= sabflg or sab[ia+1][ib]   <= sabflg or 
       sab[ia][ib+1] <= sabflg or sab[ia+1][ib+1] <= sabflg ){
    return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
       sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
  }

  return 0.0;
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
   if (a < amin) a = amin;
   sigc = sqrt(ep/e) / tev / 2;
   b = std::abs(bb);
   c = sqrt(4*M_PI);

   // tabulated s(alpha,beta).
   if (iinc != 2) { 
     std::cout << "200" << std::endl;
     // free-gas scattering.
     // 200 continue
      if (iinc != 1) {
        std::cout << "300" << std::endl;
        std::cout << "illegal option in sig! oh no!" << std::endl;
        sigVal = 0;
        return sigVal;
      }

      s = 0;
      arg = (a+bb)*(a+bb)/(4*a);
      if (-arg > sabflg) s = exp(-arg)/(c*sqrt(a));
      sigVal = sigc * sb * s;
      if (sigVal < sigmin) sigVal = 0;
      return sigVal;

    }

   if (lat == 1) b = b * tev / tevz;
   if (lat == 1) a = a * tev / tevz;
   if (a > alpha[nalpha-1]) std::cout << "go to 170 A" << std::endl;
   if (a > alpha[nalpha-1]) { 
     return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
       sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
   }

   if (lasym == 1) {
      bbm = bb;
      if (lat == 1) bbm = bb * tev / tevz;
      if (bbm > beta[nbeta-1]) { std::cout << "go to 170 B" << std::endl; }
      if (bbm > beta[nbeta-1]) {
        return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
          sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
      }
      if (bbm < beta[0]) { std::cout << "go to 170 C" << std::endl; }
      if (bbm < beta[0]) {
        return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
          sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
      }
   }
   else {
     if (b > beta[nbeta-1]) { std::cout << "go to 170 D" << std::endl; }
     if (b > beta[nbeta-1]) {
        return do170( lat, a, b, tevz, tev, teff, teff2, az, az2, sigc, 
          sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
     }
   } 
   nb1 = nbeta - 1;
   na1 = nalpha - 1;
   bbb = b;
   if (lasym == 1 and bb < 0) bbb = -b;
   for ( size_t i = 0; i < nb1; ++i ){
      ib = i;
      if (bbb < beta[i+1]) break;
   }
   for ( size_t i = 0; i < na1; ++i ) {
     ia = i;
     if ( a < alpha[i+1] ){ break; }
   }

   if (cliq == 0 or a >= alpha[0]) {
     std::cout << "150 A" << std::endl;

     std::cout << " bb:   " << bb << std::endl;
     return do150( a, az, test2, b, ia, ib, sabflg, lat, tevz, tev, sab, 
       teff, teff2, az2, sigc, sb, sb2, e, tfff, tfff2, arg, sigmin, u, c, bb, 
       nalpha, nbeta, s, alpha, beta, bbb );

   }

   if ( lasym == 1 ) {
     std::cout << "150 B" << std::endl;
     return do150( a, az, test2, b, ia, ib, sabflg, lat, tevz, tev, sab, 
       teff, teff2, az2, sigc, sb, sb2, e, tfff, tfff2, arg, sigmin, u, c, bb, 
       nalpha, nbeta, s, alpha, beta, bbb );
    }



   if (b > test1) {
     std::cout << "150 C" << std::endl;
     return do150( a, az, test2, b, ia, ib, sabflg, lat, tevz, tev, sab, 
       teff, teff2, az2, sigc, sb, sb2, e, tfff, tfff2, arg, sigmin, u, c, bb, 
       nalpha, nbeta, s, alpha, beta, bbb );
    }



   s = sab[0][0] + log(alpha[0]/a)/2 - cliq*b*b/a;

   if (s < sabflg) s = sabflg;

   // go to 160
   std::cout << "160" << std::endl;
  sigVal = 0;
  if (s-bb/2 > sabflg) sigVal = exp(s-bb/2);
  sigVal = sigc * sb * sigVal;
  if (sigVal < sigmin) sigVal = 0;
  return sigVal;


   // free-gas scattering.
  // 200 continue
  // if (iinc != 1) go to 300
   s = 0;
   arg = (a+bb)*(a+bb)/(4*a);
   if (-arg > sabflg) s = exp(-arg)/(c*sqrt(a));
   sigVal = sigc * sb * s;
   if (sigVal < sigmin) sigVal = 0;
   return sigVal;

   // other options not yet implemented.
  // 300 continue
  // call error('sig','illegal option.',' ')
   sigVal = 0;
   return sigVal;
  }


