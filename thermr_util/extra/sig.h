#include <iostream>
#include <vector>

auto do170( int lat, double a, double b, double tevz, double tev, 
    double rtev, double teff, double teff2, double az, double az2, 
    double sigc, double sb, double sb2, double arg, double e, 
    double tfff, double tfff2, double sigmin, double u, 
    double c, double bb, int sabflg ){
   // short collision time for large beta or alpha
   // 170 continue
   double sigVal;
   if (lat == 1) b = b * tevz * rtev;
   if (lat == 1) a = a * tevz * rtev;
   tfff = teff;
   tfff2 = teff2;

   double s = 0;
   arg = (a-b)*(a-b)*tev/(4*a*tfff)+(b+bb)/2;
   if (-arg > sabflg) s = exp(-arg) / (c*sqrt(a*tfff*rtev));
   sigVal = sigc * sb * s;
   if (sb2 > 0) {
      a = a * az / az2;
      arg = (a-b)*(a-b)*tev/(4*a*tfff2)+(b+bb)/2;
      double s2 = 0;
      if (-arg > sabflg) s2 = exp(-arg)/(c*sqrt(a*tfff2*rtev));
      sigVal = sigVal + sigc * sb2 * s2;
    }
    if (std::abs(e-10) < .01 and std::abs(u-.99219) < .0001) {
    }
   if (sigVal < sigmin) sigVal = 0;
   std::cout << sigVal << std::endl;
   return sigVal; 
}



auto do150( double a, double az, double test2, double b, int ia, int ib,
  int sabflg, int lat, double tevz, double tev, 
  std::vector<std::vector<double>> sab, double rtev, double teff, double teff2,
  double az2, double sigc, double sb, double sb2, double e, double tfff, 
  double tfff2, double arg, double sigmin, double u, double c, int bb ){  //150 continue
  //if (a*az < test2 and b < test2) go to 155
  if (a*az < test2 and b < test2) { 
    std::cout << "155" << std::endl; 
  } // go to 155

  //if (sab(ia,ib) <= sabflg) go to 170
  if (sab[ia][ib] <= sabflg){ // go to 170
    return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
       sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
  }

  //if (sab(ia+1,ib) <= sabflg) go to 170
  if (sab[ia+1][ib] <= sabflg){// go to 170
    return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
       sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
  }

  //if (sab(ia,ib+1) <= sabflg) go to 170
  if (sab[ia][ib+1] <= sabflg){// go to 170
    return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
       sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
  }

  //if (sab(ia+1,ib+1) <= sabflg) go to 170
  if (sab[ia+1][ib+1] <= sabflg){// go to 170
    return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
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
   int nb1, na1,i,ib,ia;
   double rtev,bb,a,sigc,b,c,bbb,s,s1,s2,s3,arg,tfff,tfff2,rat,bbm;
   double sigmin=1.e-10, sabflg=-225.e0,amin=1.e-6,test1=0.2e0,test2=30.e0;
   double sigVal;

   int nbeta = beta.size(), nalpha = alpha.size();

   // common factors.
   rtev = 1 / tev;
   bb = (ep-e) * rtev;
   a = (e+ep-2*u*sqrt(e*ep)) / (az*tev);
   if (a < amin) a = amin;
   sigc = sqrt(ep/e) * rtev / 2;
   b = std::abs(bb);
   c = sqrt(4*M_PI);

   // tabulated s(alpha,beta).
   //if (iinc != 2) go to 200
   if (lat == 1) b = b * tev / tevz;
   if (lat == 1) a = a * tev / tevz;
   if (a > alpha[nalpha-1]) std::cout << "go to 170 A" << std::endl;
   if (a > alpha[nalpha-1]) { 
     return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
       sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
   }

   if (lasym == 1) {
      bbm = bb;
      if (lat == 1) bbm = bb * tev / tevz;
      if (bbm > beta[nbeta-1]) { std::cout << "go to 170 B" << std::endl; }
      if (bbm > beta[nbeta-1]) {
        return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
          sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
      }
      if (bbm < beta[0]) { std::cout << "go to 170 C" << std::endl; }
      if (bbm < beta[0]) {
        return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
          sb, sb2, arg, e, tfff, tfff2, sigmin, u, c, bb, sabflg );
      }
   }
   else {
     if (b > beta[nbeta-1]) { std::cout << "go to 170 D" << std::endl; }
     if (b > beta[nbeta-1]) {
        return do170( lat, a, b, tevz, tev, rtev, teff, teff2, az, az2, sigc, 
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
      //if (bbb < beta[i+1]) return sigVal;
   }
   for ( size_t i = 0; i < na1; ++i ) {
     ia = i;
     if ( a < alpha[i+1] ){ break; }
     //if ( a < alpha[i+1] ){ return sigVal; }
   }

   if (cliq == 0 or a >= alpha[0]) std::cout << "150 A" << std::endl;
   // if (cliq == 0 or a >= alpha[0]) go to 150
   if ( lasym == 1 ) std::cout << "150 B" << std::endl;
   // if ( lasym == 1 ) go to 150
   if (b > test1) std::cout << "150 C" << std::endl;
   // if (b > test1) go to 150
   s = sab[0][0] + log(alpha[0]/a)/2 - cliq*b*b/a;
   if (s < sabflg) s = sabflg;
   // go to 160
  
  //150 continue
   //if (a*az < test2 and b < test2) go to 155
   //if (sab(ia,ib) <= sabflg) go to 170
   //if (sab(ia+1,ib) <= sabflg) go to 170
   //if (sab(ia,ib+1) <= sabflg) go to 170
   //if (sab(ia+1,ib+1) <= sabflg) go to 170
  //155 continue
   if (ia+1 == nalpha) ia = ia - 1;
   if (ib+1 == nbeta) ib = ib - 1;
   //call terpq(alpha(ia),sab(ia,ib),alpha(ia+1),sab(ia+1,ib),&
   //  alpha(ia+2),sab(ia+2,ib),a,s1)
   //call terpq(alpha(ia),sab(ia,ib+1),alpha(ia+1),sab(ia+1,ib+1),&
   //  alpha(ia+2),sab(ia+2,ib+1),a,s2)
   //call terpq(alpha(ia),sab(ia,ib+2),alpha(ia+1),sab(ia+1,ib+2),&
   //  alpha(ia+2),sab(ia+2,ib+2),a,s3)
   //call terpq(beta(ib),s1,beta(ib+1),s2,beta(ib+2),s3,bbb,s)
  //160 continue
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


