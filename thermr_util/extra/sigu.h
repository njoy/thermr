#include <iostream>
#include <vector>


bool do150( std::vector<double>& x, std::vector<double>& y, int i, double yt,
   double tol ){
   bool failed;
   // compare linear approximation to true function
   
   //if (i == imax) go to 160
   if (i == 20) { failed = false; return failed; }

   //if (i == 3 and 0.5*(y[i-1]+y[i])*(x[i-1]-x[i]) < tolmin) go to 160
   if (i == 3 and .5*(y[i-1]+y[i])*(x[i-1]-x[i]) < 1e-6 ) { failed = false; return failed; } 

   double xm = 0.5*(x[i-1]+x[i]);
   //xm = sigfig(xm,8,0)
   //if (xm <= x[i] or xm >= x[i-1]) go to 160
   double ym = 0.5*(y[i-1]+y[i]);
   //yt=sig(e,xm,u,tev,nalpha,alpha,nbeta,beta,sab)
   double test = tol*std::abs(yt);
   if (std::abs(yt-ym) <= test) { failed = false; return failed; }

   // point fails
   i = i + 1;
   x[i] = x[i-1];
   y[i] = y[i-1];
   x[i-1] = xm;
   y[i-1] = yt;
   //go to 150
   failed = true;
   return failed;

}



int do160( int& i, int& j, std::vector<double>& s, double& sum, 
   std::vector<double>& x, std::vector<double>& y, int& jbeta, 
   double& xl, double& yl, std::vector<double>& beta, int& nemax, double bmax ) {
   // point passes
   j = j + 1;
   s[2*j+1-1] = x[i-1];
   s[2*j+2-1] = y[i-1];
   if (j > 1) sum = sum + (y[i-1]+yl) * (x[i-1]-xl);
   xl = x[i-1];
   yl = y[i-1];
   //if (j >= nemax-1) { do170( s, sum, j ); }         // go to 170
   if (j >= nemax-1) { return 170; }         // go to 170
   if (jbeta > 0) {
  //   if (beta[jbeta] > bmax) { do170( s, sum, j ); } // go to 170
     if (beta[jbeta-1] > bmax) { return 170; } // go to 170
   }

   // continue bin loop and linearization loop
   i = i - 1;
   //if (i > 1) go to 150
   jbeta = jbeta + 1;
   //if (jbeta <= beta.size() ) { do111( x, y ); }     // go to 111
   if (jbeta <= int(beta.size()) ) {return 111; }     // go to 111
   // if (i.eq.1) go to 160
   return 160;
   

}



std::tuple<std::vector<double>,std::vector<double>,std::vector<double>>
  sigu( double e, double u, double tev, double tevz,
    std::vector<double> alpha, std::vector<double> beta, 
    std::vector<std::vector<double>> sab, double tolin, double az, 
    int nemax, int lasym, int lat ){

  /*-------------------------------------------------------------------
   * Compute the secondary energy distribution scattering for cosine u.
   * Uses linear reconstruction with the cross section from function sig.
   *-------------------------------------------------------------------
   */
   int i, j, jbeta, imax = 20;
   double sum, xl, yl, xm, ym, test, yt = 0.0, tol, root1, root2;
   std::vector<double> x ( imax ), y ( imax );
   double tolmin = 1.e-6;
   double bmax = 20;

   // constant factors
   tol = tolin;
   std::vector<double> s ( 2 * nemax, 0.0 );

   root1 = ( u*sqrt(e) + sqrt( u*u*e + (az-1) * (az+1) * e ) ) / (az+1);
   root2 = ( u*sqrt(e) - sqrt( u*u*e + (az-1) * (az+1) * e ) ) / (az+1);

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

   int counter = 0;
   // set up next panel
   while ( counter < 4 ){
       // 111 
       std::cout << "111" << std::endl;
       x[1] = x[0];
       y[1] = y[0];

    do {
      // 113 
      std::cout << "113" << std::endl;
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

      if ( x[0] <= x[1] ){
        jbeta = jbeta + 1;
      }
    } while ( x[0] <= x[1] ); // if x[0] > x[1] we leave and go to 116


     // 116 continue
     std::cout << "116" << std::endl;
     if (u < 0 and root1*root1 > 1.01*x[1] and root1*root1 < x[0]) {
       x[0] = root1*root1;
     }
     //x(1)=sigfig(x(1),8,0)
     //y(1)=sig(e,x(1),u,tev,nalpha,alpha,nbeta,beta,sab)
     i = 2;


     // 150
     bool failed;
     do {
       std::cout << "150" << std::endl;
       failed = do150( x, y, i, yt, tol );
     } while ( failed );


     // 160
     int whereTo;
     do {
       std::cout << "160" << std::endl;
       whereTo = do160( i, j, s, sum, x, y, jbeta, xl, yl, beta, nemax, bmax );
       std::cout << "bring me to " << whereTo << std::endl;
  //     std::cout << s[0] << "   " << s[1] << "   " << s[2] << std::endl;
  //     std::cout << s[3] << "   " << s[4] << "   " << s[5] << std::endl;
  //     std::cout << s[6] << "   " << s[7] << "   " << s[8] << std::endl;
  //     std::cout << s[9] << std::endl;
       std::cout << "   " << std::endl;
     } while ( whereTo == 160 );

     if ( whereTo == 170 ){ 
       std::cout << "170" << std::endl;
       s[0] = sum;
       s[1] = j;
       //for ( auto entry : s ) { std::cout << "---------- " << entry << std::endl; }
       return { x, y, s };
     }
     counter += 1;
  }
  if ( counter > 4) {std::cout << "HAD TO BREAK" << std::endl; }

   std::cout << "OH NO DONT END LIKE THIS" << std::endl;
  return { x, y, s };
}

