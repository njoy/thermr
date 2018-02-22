#include <iostream>
#include <vector>
#include "sig.h"
#include "legndr.h"

auto do190( double& yl, std::vector<double>& y, std::vector<double>& x,
  double& gral, int nlin, double& xl, double& fract,
  std::vector<double>& s, double xn, int& i, double& xil, int& j, int nbin,
  std::vector<double>& p, double& sum, int nl){

  std::cout << "190" << std::endl;
  // 190 continue
  double yn = yl + (y[i-1]-yl) * (xn-xl) * xil;
  gral += (xn-xl) * ( yl * 0.5* (xn+xl) + (y[i-1]-yl) * 
    xil * (-xl*0.5*(xn + xl) + (1.0/3.0) * ( xn*xn + xn*xl + xl*xl ) ) );
  double xbar = gral / fract;

  // compute legendre components
  if (nlin >= 0){
    legndr(xbar,p);
    for ( int il = 1; il < nl; ++il ){
      s[il]=s[il]+p[il]/nbin;
    }
  } 
  // output equally probable angles
  else{
     s[j+1-1]=xbar;
  }

  // continue bin loop and linearization loop
  xl=xn;
  yl=yn;
  sum=0;
  gral=0;
  if (j == nbin) {
    //   std::cout << "go to 260" << std::endl;
    return 260;
  }
  if (xl < x[i-1]){
    //    std::cout << "go to 160" << std::endl;
    return 160;
  }
  return 250;
}


auto do170( double& xn, double& xl, int& j, double& yl, double& sigmin, double& fract,
  double& sum, double& ytol, std::vector<double>& x, std::vector<double>& y,
  double& xil, std::vector<double>& s, int& nl, int& nlin, 
  double& gral, int& nbin, int& i, std::vector<double>& p){

  // 170 continue
   j += 1;
   if (yl < sigmin) {
     std::cout << "go to 175" << std::endl;
   }
   double test=(fract-sum)*(y[i-1]-yl)/((x[i-1]-xl)*yl*yl);
   if (std::abs(test) > ytol){
     std::cout << "go to 175" << std::endl;
   }
   xn=xl+(fract-sum)/yl;
   if (xn > x[i-1]) { 
     std::cout << "go to 180" << std::endl;
   }
   if (xn >= xl and xn <=  x[i-1]){
     std::cout << "go to 190" << std::endl;
   }
   std::cout << "call error('sigl','no legal solution.',' ')" << std::endl;

  //175 continue
   double f=(y[i-1]-yl)*xil;
   double rf=1/f;
   double disc=(yl*rf)*(yl*rf)+2*(fract-sum)*rf;
   if (disc < 0) {
     std::cout << " write(strng,'(''disc='',1p,e12.4)') disc" << std::endl;
     std::cout << "call mess('sigl',strng,'set to abs value and continue')" << std::endl;
      disc=std::abs(disc);
   }
   if (f > 0 ) xn=xl-(yl*rf)+sqrt(disc);
   if (f < 0 ) xn=xl-(yl*rf)-sqrt(disc);
   if (xn > xl and xn <= x[i-1]){
     std::cout << "go to 190" << std::endl;
   }
   if (xn > xl and xn < (x[i-1]+ytol*(x[i-1]-xl))){
     std::cout << "go to 180" << std::endl;
   }
   std::cout << "call error('sigl','no legal solution (quadratic path).',' ')" << std::endl;
  // 180 continue
   xn=x[i-1];
  // 190 continue
   double yn=yl+(y[i-1]-yl)*(xn-xl)*xil;
   gral=gral+(xn-xl)*(yl*0.5*(xn+xl)
     +(y[i-1]-yl)*xil*(-xl*0.5*(xn+xl)
     +(1/3)*(xn*xn+xn*xl+xl*xl)));
   double xbar=gral/ fract;

   // compute legendre components
   if (nlin >= 0){
      // call legndr(xbar,p,nl);
      for ( int il = 2; il < nl; ++il ){
        s[il-1]=s[il-1]+p[il-1]/nbin;
      }
   } 
   // output equally probable angles
   else{
      s[j+1-1]=xbar;
  }

   // continue bin loop and linearization loop
   xl=xn;
   yl=yn;
   sum=0;
   gral=0;
   if (j == nbin) {
     std::cout << "go to 260" << std::endl;
   }
   if (xl < x[i-1]){
     //std::cout << "go to 160" << std::endl;
   }
   return 0;
}



auto do160( double& sum, double& xl, double& yl, std::vector<double>& x,
    std::vector<double>& y, int& i, int& j, double& gral, double& fract, 
    std::vector<double>& s, std::vector<double>& p, int& nl, 
    int& nlin, int& nbin, double& shade, double& sigmin, double& ytol){
  // 160 --> either 250, 165, or 170. These mostly somehow get to 190, and then
  // either return completely or go to 250. 250 will go back to 150 or 160. So
  // hopefully getting this done will help with the logic web
 
  // check bins for this panel
  //160 continue
  
  // ----- So chances are that we're going to have to go to 250 at some point
  // ----- during this evaluation, so let's separate out the change that we 
  // ----- have something else to do before then
   double add=0.5*(y[i-1]+yl)*(x[i-1]-xl);
   double xil, xn;

   if ( x[i-1] != xl ){
     xil = 1 / (x[i-1]-xl);
     if (i == 1 and j == nbin-1){
       // 165 continue
       std::cout << "165" << std::endl;
       xn = x[i-1];
       j += 1;
       return do190( yl, y, x, gral, nlin, xl, fract, s, xn, i, xil, j, 
           nbin, p, sum, nl);

      }
      if (sum+add >= fract*shade and j < nbin-1) {
        std::cout << "go to 170" << std::endl;
        return do170( xn, xl, j, yl, sigmin, fract, sum, ytol, x, y, xil, s, 
            nl, nlin, gral, nbin, i, p);

     }
     sum += add;
     gral+=0.5*(yl*x[i-1]-y[i-1]*xl)*(x[i-1]+xl)
       +(1.0/3)*(y[i-1]-yl)*(x[i-1]*x[i-1]+x[i-1]*xl+xl*xl);
   }

   // 250 continue
   //std::cout << "250" << "     " << x[i-1] << "    " << i  << std::endl;
   xl = x[i-1];
   yl = y[i-1];
   i = i - 1;
   if (i  > 1) { return 150; }
   if (i == 1) { return 160; }

   return 0;  // 260 continue
 
  
}






auto sigl( double e, double ep, double tev, std::vector<double> alpha,
  std::vector<double> beta, std::vector<double>& s, 
  std::vector<std::vector<double>> sab, 
  double tolin, int nlin, double az, double az2,  
  double teff, double teff2, int lat, double tevz, int lasym, int iinc, 
  double cliq, double sb, double sb2, int nbin ){

 /*-------------------------------------------------------------------
  * Compute the cross section and legendre components or equally-
  * probable angles for the scattering from e to ep.  Uses linear
  * reconstruction of the angular distribution computed by sig.
  *-------------------------------------------------------------------
  */
  
  int nl,i,j,il,nbeta = beta.size(), nalpha = alpha.size();
  double b,seep,sum,xl,yl,ymax,xm,ym,test,test2,fract,gral,add,xil,xn,f,rf,disc,yn,xbar,yt,tol,s1bb,xtol=.00001e0, ytol=.001e0, sigmin=1.e-32, eps=1.e-3, shade=.99999999e0;
  int imax=20;
  std::vector<double> x(imax), y(imax),p(nlin);

   // constant factors
   b=(ep-e)/tev;
   tol=0.5*tolin;
   nl=nlin;
   if (nl < 0) nl=-nl;
   b=std::abs(b);
   if (lat == 1 and iinc == 2) b=b*tev/tevz;
   s1bb=sqrt(1+b*b);
   if (ep != 0) seep=1/sqrt(e*ep);

   // adaptive calculation of cross section
   i=3;
   sum=0;
   x[2]=-1;
   xl=x[2];
   y[2]=sig(e,ep,x[2],tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
   yl=y[2];
   if (ep == 0) x[1]=0;
   if (ep != 0) x[1]=0.5*(e+ep-(s1bb-1)*az*tev)*seep;
   if (std::abs(x[1]) > 1-eps) x[1]=0.99e0;
   y[1]=sig(e,ep,x[1],tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
   x[0]=1;
   y[0]=sig(e,ep,x[0],tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
   ymax=y[1];
   if (y[0] > ymax) ymax=y[0];
   if (y[2] > ymax) ymax=y[2];
   if (ymax < eps) { ymax=eps; }

  //110 continue
  // if (i == imax) go to 120
   while (true){ 
     while ( true ){
     std::cout << "110" << std::endl;
     xm=0.5*(x[i-1]+x[i-2]);
     ym=0.5*(y[i-1]+y[i-2]);
     yt=sig(e,ep,xm,tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
     test=tol*std::abs(yt)+tol*ymax/50;
     test2=ym+ymax/100;
     if (std::abs(yt-ym) <= test and std::abs(y[i-2]-y[i-1]) <= test2 and 
       (x[i-2]-x[i-1]) < 0.5){
       break;
     }
     if (x[i-2]-x[i-1] < xtol) {
       break;
     }
     i = i + 1;
     x[i-1]=x[i-2];
     y[i-1]=y[i-2];
     x[i-2]=xm;
     y[i-2]=yt;
   }
   while (true) {
     std::cout << "120" << std::endl;
     // 120 continue
     sum = sum+0.5*(y[i-1]+yl)*(x[i-1]-xl);
     xl=x[i-1];
     yl=y[i-1];
     i=i-1;

     if ( i != 1 ){ break; }
   }
   if ( i < 1 ){ break; }
  }
   s[0]=sum;
   if (sum <= sigmin) {
     for ( int il = 0; il < nl-1; ++il ){ s[il] = 0.0; }
     return;
   } 

   // prime stack for equally-probable angles
   // 130 continue
   std::cout << "130" << std::endl;
   nbin=nl-1;
   fract=sum/nbin;

   sum = 0;
   gral = 0;
   for ( int il = 1; il < nl; ++il ){ 
     s[il] = 0; 
   }
   j=0;

   // adaptive linearization
   i = 3;
   x[2]=-1;
   xl=x[2];
   y[2]=sig(e,ep,x[2],tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
   if (ep == 0 ) x[1]=0;
   if (ep != 0 ) x[1]=0.5*(e+ep-(s1bb-1)*az*tev)*seep;
   if (abs(x[1]) > 1-eps) x[1]=0.99e0;
   y[1]=sig(e,ep,x[1],tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
   x[0]=+1;
   y[0]=sig(e,ep,x[0],tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
   ymax=y[0];
   if (y[1] > ymax) ymax=y[1];
   if (y[2] > ymax) ymax=y[2];
   if (ymax < eps) ymax=eps;
   
   int counter = 0;
   while (true){
     while (true){
       //150 continue
       std::cout << "150" << std::endl;
       if (i == imax) { std::cout << "go to 160!!!" << std::endl; break; }
       else {
         //if (i == imax) go to 160
         xm=0.5*(x[i-1]+x[i-2]);
         ym=0.5*(y[i-1]+y[i-2]);
         //std::cout << "----------------- " << i << "   " << xm << "    " << ym  << std::endl;
         yt=sig(e,ep,xm,tev,tevz,alpha,beta,sab,az,az2,lat,iinc,lasym,cliq,sb,sb2,teff,teff2);
     //std::cout << x[i-2]  << "   " << x[i-1] << std::endl;
         test=tol*std::abs(yt)+tol*ymax/50;
         test2=ym+ymax/100;
         if (abs(yt-ym) <= test and 
             std::abs(y[i-1]-y[i-2]) <= test2 and 
            (x[i-2]-x[i-1]) < 0.5) {
            break; 
         }
         if (x[i-2]-x[i-1] < xtol){
            break; 
         }
         i=i+1;
         x[i-1]=x[i-2];
         y[i-1]=y[i-2];
         x[i-2]=xm;
         y[i-2]=yt;
      //std::cout << "    " << x[0] << "     " << x[1] << "      " << x[2] << "    " << i << "   " << j << std::endl;
      counter +=1 ;
      if (counter > 3 ){return;}
         //go to 150
       }
     }

    // check bins for this panel
    //160 continue
    //while (true ){
    int whatToDo; 
     do {
      std::cout << "160" << std::endl;

       whatToDo = do160( sum, xl, yl, x, y, i, j, gral, fract, s, p, nl,  
          nlin, nbin, shade, sigmin, ytol);
        
       if ( whatToDo == 260 ){ std::cout << "260" << std::endl; return; }
 
    } while ( whatToDo == 160 );

    if ( whatToDo != 150 ){ break; }

   } // when near 150
}











