#include "sig.h"

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
   int i,j,jbeta;
   double sum,xl,yl,xm,ym,test,yt,tol,root1,root2;
   int imax=20;
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
   //y[1-1]=sig(e,x[1-1],u,tev,nalpha,alpha,nbeta,beta,sab)
   y[1-1] = sig(e,x[1-1],xl,tev,alpha,beta,sab,bbm,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   jbeta=-beta.size();
   if (lasym > 0) jbeta=1;
   j=0;
   xl=0;
   yl=0;

   // set up next panel
   // 111 continue
   x[2-1]=x[1-1]; 
   y[2-1]=y[1-1];
  // 113 continue
   if (jbeta == 0) jbeta=1;
   if (jbeta <= 0) {
      if (lat == 1) {
         x[1-1]=e-beta[-jbeta-1]*tevz;
      }
      else { 
         x[1-1]=e-beta[-jbeta-1]*tev;
      }
      // x[1-1]=sigfig(x[1-1],8,0)
      // if (x[1-1] == e) x[1-1]=sigfig(e,8,-1)
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
   jbeta=jbeta+1;
   // go to 113
   
   // 116 continue
   if (u < zero and root1*root1 > 1.01*x[2-1] and root1*root1 < x[1-1]) {
      x[1-1]=root1*root1;
   }
   // x[1-1]=sigfig(x[1-1],8,0)
   // y[1-1]=sig(e,x[1-1],u,tev,nalpha,alpha,nbeta,beta,sab)
   y[1-1] = sig(e,x[1-1],xl,tev,alpha,beta,sab,bbm,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   i=2;

   // compare linear approximation to true function
   // 150 continue
   // if (i == imax) go to 160
   // if (i > 3 and half*(y[i-1-1]+y[i-1])*(x[i-1-1]-x[i-1]) < tolmin) go to 160
   xm=half*(x[i-1-1]+x[i-1]);
   // xm=sigfig(xm,8,0);
   // if (xm <= x[i-1] or xm.ge.x[i-1-1]) go to 160
   ym=half*(y[i-1-1]+y[i-1]);
   // yt=sig(e,xm,u,tev,nalpha,alpha,nbeta,beta,sab)
   yt = sig(e,xm,xl,tev,alpha,beta,sab,bbm,az,tevz,lasym,az2,teff2,lat,cliq,sb,sb2,teff,iinc);
    if (abs(u-.99219) < .0001) {
       if (abs(e-10) < .01) {
       }
    }
   test=tol*abs(yt);
   // if (abs(yt-ym) <= test) go to 160

   // point fails
   i=i+1;
   x[i-1]=x[i-1-1];
   y[i-1]=y[i-1-1];
   x[i-1-1]=xm;
   y[i-1-1]=yt;
   // go to 150

   // point passes
   // 160 continue
   j=j+1;
   s[2*j+1-1]=x[i-1];
   s[2*j+2-1]=y[i-1];
   if (j > 1) sum=sum+(y[i-1]+yl)*(x[i-1]-xl);
   xl=x[i-1];
   yl=y[i-1];
   // if (j >= nemax-1) go to 170
   if (jbeta > 0) {
     //if (beta(jbeta) > bmax) go to 170
   }

   // continue bin loop and linearization loop
   i=i-1;
   // if (i > 1) go to 150
   jbeta=jbeta+1;
   // if (jbeta <= nbeta) go to 111
   // if (i == 1) go to 160
   // 170 continue
   s[1-1]=sum;
   s[2-1]=j;

   return; 
}


