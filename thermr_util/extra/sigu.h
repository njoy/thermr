#include <iostream>
#include <vector>



void do111( std::vector<double>& x, std::vector<double>& y ){
  x[1] = x[0];
  y[1] = y[0];
}

void do113( int& jbeta, int& lat, double e, std::vector<double>& x, 
    std::vector<double>& beta, double tev, double tevz ){
   if ( jbeta == 0 ){ jbeta = 1; } 
   if ( jbeta <= 0 ){
      if (lat == 1){
         x[0] = e-beta[-jbeta-1]*tevz;
      }
      else {
         x[0] = e-beta[-jbeta-1]*tev;
      }
      //x[0] = sigfig(x(1),8,0);
      //if (x[0] == e) x[0] = sigfig(e,8,-1);
   } 
   else {
      if (lat == 1){
         x[0] = e + beta[jbeta-1]*tevz;
      }
      else {
         x[0] = e + beta[jbeta-1]*tev;
      } 
   }
}

void do170( std::vector<double>& s, double sum, int j ){
  s[0] = sum;
  s[1] = j;
}

void do160( int i, int j, std::vector<double>& s, double& sum, 
   std::vector<double>& x, std::vector<double>& y, int jbeta, 
   double& xl, double& yl, std::vector<double>& beta, int nemax, int bmax ) {
   // point passes
   j = j + 1;
   s[2*j+1] = x[i];
   s[2*j+2] = y[i];
   if (j > 1) sum = sum + (y[i]+yl) * (x[i]-xl);
   xl = x[i];
   yl = y[i];
   if (j >= nemax-1) { do170( s, sum, j ); }         // go to 170
   if (jbeta > 0) {
     if (beta[jbeta] > bmax) { do170( s, sum, j ); } // go to 170
   }

   // continue bin loop and linearization loop
   i = i - 1;
   //if (i > 1) go to 150
   jbeta = jbeta + 1;
   if (jbeta <= beta.size() ) { do111( x, y ); }     // go to 111
   // if (i.eq.1) go to 160

}



auto sigu( double e, double u, double tev, double tevz,
    std::vector<double> alpha, std::vector<double> beta, 
    std::vector<std::vector<double>> sab, double tolin, double az, 
    int nemax, int lasym, int lat ){

  /*-------------------------------------------------------------------
   * Compute the secondary energy distribution scattering for cosine u.
   * Uses linear reconstruction with the cross section from function sig.
   *-------------------------------------------------------------------
   */
   int i, j, jbeta, imax = 20;
   double sum, xl, yl, xm, ym, test, yt, tol, root1, root2;
   std::vector<double> x ( imax ), y ( imax );
   double tolmin = 1.e-6;
   double bmax = 20;

   // constant factors
   tol = tolin;
   std::vector<double> s ( 2 * nemax, 0.0 );

   root1 = ( u*sqrt(e) + sqrt( u*u*e + (az-1) * (az+1) * e ) ) / (az+1);
   root2 = ( u*sqrt(e) - sqrt( u*u*e + (az-1) * (az+1) * e ) ) / (az+1);
   //std::cout << root1 << "    " << root2 << std::endl;

   // adaptive calculation of cross section
   sum = 0;
   x[0] = 0;

   // FIX THIS YOU SHOULD REALLY IMPLEMENT SIG
   //y[0] = sig(e,x(1),u,tev,nalpha,alpha,nbeta,beta,sab);
   y[0] = 0.0;

   jbeta = -beta.size();
   if (lasym > 0) { jbeta = 1; }
   j  = 0; 
   xl = 0;
   yl = 0;

   // set up next panel
   // 111 
   std::cout << "111" << std::endl;
   do111( x, y );

   // 113 
   std::cout << "113" << std::endl;
   do113( jbeta, lat, e, x, beta, tev, tevz );


   if (x[0] > x[1] ) {   // go to 116
      // 116 continue
      std::cout << "116" << std::endl;
      if (u < 0 and root1*root1 > 1.01*x[1] and root1*root1 < x[0]) {
        x[0] = root1*root1;
      }
      //x(1)=sigfig(x(1),8,0)
      //y(1)=sig(e,x(1),u,tev,nalpha,alpha,nbeta,beta,sab)
      i = 2;

   }
   else { 
     jbeta = jbeta + 1;
     // go to 113
   } 
   /*

   !--compare linear approximation to true function
  150 continue
   if (i.eq.imax) go to 160
   if (i.gt.3.and.0.5*(y(i-1)+y(i))*(x(i-1)-x(i)).lt.tolmin) go to 160
   xm=0.5*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   if (xm.le.x(i).or.xm.ge.x(i-1)) go to 160
   ym=0.5*(y(i-1)+y(i))
   yt=sig(e,xm,u,tev,nalpha,alpha,nbeta,beta,sab)
    if (abs(u-.99219).lt..0001) then
       if (abs(e-10).lt..01) then
       endif
    endif
   test=tol*abs(yt)
   if (abs(yt-ym).le.test) go to 160

   !--point fails
   i=i+1
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=yt
   go to 150

   !--point passes
  160 continue
   j=j+1
   s(2*j+1)=x(i)
   s(2*j+2)=y(i)
   if (j.gt.1) sum=sum+(y(i)+yl)*(x(i)-xl)
   xl=x(i)
   yl=y(i)
   if (j.ge.nemax-1) go to 170
   if (jbeta.gt.0) then
     if (beta(jbeta).gt.bmax) go to 170
   endif

   !--continue bin loop and linearization loop
   i=i-1
   if (i.gt.1) go to 150
   jbeta=jbeta+1
   if (jbeta.le.nbeta) go to 111
   if (i.eq.1) go to 160
  170 continue
   s(1)=sum
   s(2)=j

   return
   end subroutine sigu

   */
}

