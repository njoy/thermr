#include "sig.h"


auto do_113_116( int& jbeta, const int& lat, std::vector<double>& x, 
  std::vector<double>& y, const double& e, const double& tev, const double& tevz,
  const double& root1, const double& u, double& xl, double& yl, 
  std::vector<double>& alpha, std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, const double& bbm, const double& az, 
  const int lasym, const double& az2, const double& teff, const double& teff2, 
  const double& cliq, const double& sb, const double& sb2, int& i, const int& iinc){
  double zero = 0;
  while (true){
    // 113 continue
    std::cout << 113 << std::endl;
    if (jbeta == 0) jbeta=1;
    if (jbeta <= 0) {
      if (lat == 1) {
        x[1-1]=e-beta[-jbeta-1]*tevz;
      }
      else { 
        x[1-1]=e-beta[-jbeta-1]*tev;
      }
    }
    else {
      if (lat == 1) {
        x[1-1]=e+beta[jbeta-1]*tevz;
      }
      else {
        x[1-1]=e+beta[jbeta-1]*tev;
      }
    }
    // if (x[1-1] > x[2-1]) go to 116
    if (x[1-1] > x[2-1]) break;

    jbeta=jbeta+1;
    // go to 113
  } // This is doing 113
   
  // 116 continue
  std::cout << 116 << std::endl;
  if (u < zero and root1*root1 > 1.01*x[2-1] and root1*root1 < x[1-1]) {
    x[1-1]=root1*root1;
  }
  y[1-1] = sig(e,x[1-1],xl,tev,alpha,beta,sab,bbm,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  i=2;
} 




auto sigu( int nlin, int nemax, double& e, double& u, double& ep, double& tev, 
  std::vector<double>& alpha, std::vector<double>& beta, 
  std::vector<std::vector<double>>&sab, std::vector<double>&s, double& tolin, 
  double& az, double& tevz, int iinc, int lat, double bbm, int lasym, 
  double& az2, double& teff2, double& cliq, double& sb, double& sb2, 
  double& teff ){

  /*-------------------------------------------------------------------
   * Compute the secondary energy distribution scattering for cosine u.
   * Uses linear reconstruction with the cross section from function sig.
   *-------------------------------------------------------------------
   */
   int i,j,jbeta,imax=20;
   double sum,xl,yl,xm,ym,test,yt,tol,root1,root2;
   std::vector<double> x(imax), y(imax);
   double zero=0, half=.5, tolmin=1.e-6, bmax=20;

   // constant factors
   tol=tolin;
   for ( int i = 0; i < 2*nemax; ++i ){
     s[i] = 0;
   }

   root1=(u*sqrt(e)+sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);
   root2=(u*sqrt(e)-sqrt(u*u*e+(az-1)*(az+1)*e))/(az+1);

   // adaptive calculation of cross section
   sum=0;
   x[1-1]=0;
   y[1-1] = sig(e,x[1-1],xl,tev,alpha,beta,sab,bbm,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   jbeta=-beta.size();
   if (lasym > 0) jbeta=1;
   j=0;
   xl=0;
   yl=0;

   // set up next panel
   while (true){  
     // 111 continue
     std::cout << 111 << std::endl;
     x[2-1]=x[1-1]; 
     y[2-1]=y[1-1];




     do_113_116( jbeta, lat, x, y, e, tev, tevz, root1, u, xl, yl, alpha, beta, sab, bbm, az, lasym, az2, teff, teff2, cliq, sb, sb2, i, iinc);




     // compare linear approximation to true function
     // 150 continue
     bool do150 = true;
     while (do150){
       std::cout << 150 << std::endl;
       // if (i == imax) go to 160
       if (i != imax){
         // if (i > 3 and half*(y[i-1-1]+y[i-1])*(x[i-1-1]-x[i-1]) < tolmin) go to 160
         if (i <= 3 or half*(y[i-1-1]+y[i-1])*(x[i-1-1]-x[i-1]) >= tolmin) {
           xm=half*(x[i-1-1]+x[i-1]);
           // if (xm <= x[i-1] or xm.ge.x[i-1-1]) go to 160
           if (xm > x[i-1] and xm > x[i-1-1]){
             ym=half*(y[i-1-1]+y[i-1]);
             // yt=sig(e,xm,u,tev,nalpha,alpha,nbeta,beta,sab)
             yt = sig(e,xm,xl,tev,alpha,beta,sab,bbm,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
             test=tol*abs(yt);
      
             if (abs(yt-ym) <= test) {
               std::cout << "go to 160" << std::endl;
             }
             if (abs(yt-ym) > test) {
               // point fails
               i=i+1;
               x[i-1]=x[i-1-1];
               y[i-1]=y[i-1-1];
               x[i-1-1]=xm;
               y[i-1-1]=yt;
               // go to 150
               continue;
             }
           }
         }
       }
  
       // point passes
       // 160 continue
       std::cout << 160 << std::endl;
       j=j+1;
       s[2*j+1-1]=x[i-1];
       s[2*j+2-1]=y[i-1];
       if (j > 1) sum=sum+(y[i-1]+yl)*(x[i-1]-xl);
       xl=x[i-1];
       yl=y[i-1];

       // if (j >= nemax-1) go to 170
       if (j >= nemax-1 or (jbeta > 0 and beta[jbeta-1] > bmax )) {
         s[1-1]=sum;
         s[2-1]=j;
         return; 
      }

      // continue bin loop and linearization loop
      i=i-1;

      // if (i > 1) go to 150
      if (i > 1) {
        std::cout << "got to 150" << std::endl;
        continue;
      }
      std::cout << "no longer doing 150" << std::endl;
      break;
    } 

    jbeta=jbeta+1;

    // if (jbeta <= nbeta) go to 111
    if (jbeta <= beta.size()) continue;
 
    // if (i == 1) go to 160
    if (i == 1) {
      std::cout << "go to 160" << std::endl;
    }

    // 170 continue
    std::cout << 170 << std::endl;
    s[1-1]=sum;
    s[2-1]=j;
    return; 
  }
}


