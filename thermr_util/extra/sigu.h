#include <iostream>
#include <vector>


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
   x[1] = x[0];
   y[1] = y[0];

   // 113 
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

   if (x[0] > x[1] ) {   // go to 116
      // 116 continue
      if (u < 0 and root1*root1 > 1.01*x[1] and root1*root1 < x[0]) {
        x[0] = root1*root1;
      }
      //x(1)=sigfig(x(1),8,0)
      //y(1)=sig(e,x(1),u,tev,nalpha,alpha,nbeta,beta,sab)
      i = 2;
      std::cout << "HERE   " << x[0] << std::endl;

   }
}

