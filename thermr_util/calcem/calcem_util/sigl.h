#include "sig.h"

auto sigl( int nlin, int nalpha, int nbeta, int nlmax, double e, double ep,
    double tev, std::vector<double>& alpha, std::vector<double>& beta,
    std::vector<std::vector<double>>& sab, std::vector<double>& s, double tolin,
    double az, double tevz, int iinc, int lat, double bbm, 
  int lasym, double az2, double teff2, double cliq, double sb,
  double sb2, double teff ){

  /*-------------------------------------------------------------------
   * This is called by calcem, and uses sig.
   * * * Attempted description: * * *
   * Calcem has to turn sig(E->E',mu) --> sig(E->E') by integrating 
   * over mu. So we need to subdivide the cosine range until the 
   * actual angular function (given here) is within a nice tolerance
   *-------------------------------------------------------------------
   * * * Actual description: * * *
   * Compute the cross section and legendre components or equally-
   * probable angles for the scattering from e to ep.  Uses linear
   * reconstruction of the angular distribution computed by sig.
   *-------------------------------------------------------------------
   */
   int nl,i,j,il,nbin;
   double b,seep,sum,xl,yl,ymax,xm,ym,test,test2;
   double rnbin,fract,gral,add,xil,xn,f,rf,disc,yn,xbar,rfract;
   double yt,tol,s1bb;
   int imax=20;
   std::vector<double> x(imax), y(imax), p(nlin);
   // character(60)::strng
   double zero = 0, one = 1, half = 0.5, third = (1.0/3.0);
   double xtol = 0.00001;
   double ytol = 0.001;
   double sigmin = 1.0e-32, eps = 1.0e-3;
   double shade = 0.99999999; 

   // constant factors
   b=(ep-e)/tev;
   tol=half*tolin;
   nl=nlin;
   if (nl < 0) nl=-nl;
   b=abs(b);
   if (lat == 1 and iinc == 2) b=b*tev/tevz;
   s1bb=sqrt(1+b*b);
   if (ep != zero) seep=1/sqrt(e*ep);

   // adaptive calculation of cross section
   i=3;
   sum=0;
   x[3-1]=-1;
   xl=x[3-1];
   y[3-1]=sig(e,ep,x[3-1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   yl=y[3-1];
   if (ep == zero) x[2-1]=0;
   if (ep != zero) x[2-1]=half*(e+ep-(s1bb-1)*az*tev)*seep;
   if (abs(x[2-1]) > 1-eps) x[2-1]=0.99;
   y[2-1]=sig(e,ep,x[2-1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   x[1-1]=+1;
   y[1-1]=sig(e,ep,x[1-1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   ymax=y[2-1];
   if (y[1-1] > ymax) ymax=y[1-1];
   if (y[3-1] > ymax) ymax=y[3-1];
   if (ymax < eps) ymax=eps;

   bool do110 = true;
   bool do110_inner = true;
   bool do120 = true;
   while ( do110 ){
     while ( do110_inner ){ 
       std::cout << "110" << std::endl;
       // 110 continue
       
       // if (i == imax) go to 120
       if (i == imax) break;

       xm=half*(x[i-1-1]+x[i-1]);
       ym=half*(y[i-1-1]+y[i-1]);

       yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
         az2,teff2,lat,cliq,sb,sb2,teff,iinc);
      
       test=tol*abs(yt)+tol*ymax/50;
       test2=ym+ymax/100;
       if (abs(yt-ym) <= test and abs(y[i-1-1]-y[i-1]) <= test2 and
         (x[i-1-1]-x[i-1]) < half) { break; } // go to 120
       if (x[i-1-1]-x[i-1] < xtol) { break; } // go to 120
       i=i+1;
       x[i-1]=x[i-1-1];
       y[i-1]=y[i-1-1];
       x[i-1-1]=xm;
       y[i-1-1]=yt;
     }  // do 110 inner. This corresponds to    // go to 110 

     while ( do120 ){
       std::cout << "120" << std::endl;

       // 120 continue
       sum=sum+half*(y[i-1]+yl)*(x[i-1]-xl);
       xl=x[i-1];
       yl=y[i-1];
       i=i-1;

       // if (i > 1) go to 110
       if (i > 1) { break; }

       // if (i == 1) go to 120
       if (i != 1) { // don't go to 120 - either go to 130 or return 

         s[1-1]=sum;
         // if (sum > sigmin) go to 130
         if (sum > sigmin) { do110 = false; break;} 
         for ( int il = 0; il < nl; ++il ){
           s[il] = 0;
         } // end do 
         std::cout << "return from 120" << std::endl;
         return;
       } // don't go to 120

      } // go to 120
   } // go to 110

   // prime stack for equally-probable angles
   // 130 continue
   std::cout << "130" << std::endl;
   nbin=nl-1;
   rnbin=one/nbin;
   fract=sum*rnbin;
   rfract=1/fract;
   sum=0;
   gral=0;
   for ( int il = 1; il < nl; ++il ){
     s[il] = 0;
   }
   j=0;

   // adaptive linearization
   i=3;
   x[3-1]=-1;
   xl=x[3-1];
   y[3-1]=sig(e,ep,x[3-1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);
   if (ep == zero) x[2-1]=0;
   if (ep != zero) x[2-1]=half*(e+ep-(s1bb-1)*az*tev)*seep;
   if (abs(x[2-1]) > 1-eps) x[2-1]=0.99e0;
   y[2-1]=sig(e,ep,x[2-1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);

   x[1-1]=+1;
   y[1-1]=sig(e,ep,x[1-1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);


   ymax=y[1-1];
   if (y[2-1] > ymax) ymax=y[2-1];
   if (y[3-1] > ymax) ymax=y[3-1];
   if (ymax < eps) ymax=eps;
   return;

   bool do150 = true;
   while ( do150 ){
   // 150 continue
     std::cout << "150" << std::endl;
     //if (i == imax) go to 160
     xm=half*(x[i-1-1]+x[i-1]);
     ym=half*(y[i-1-1]+y[i-1]);
     yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
       az2,teff2,lat,cliq,sb,sb2,teff,iinc);

     test=tol*abs(yt)+tol*ymax/50;
     test2=ym+ymax/100;
     //if (abs(yt-ym) <= test and abs(y(i-1)-y(i)) <= test2 and &
     //  (x(i-1)-x(i)) < half) go to 160
     //if (x(i-1)-x(i) < xtol) go to 160
     i=i+1;
     x[i-1]=x[i-1-1];
     y[i-1]=y[i-1-1];
     x[i-1-1]=xm;
     y[i-1-1]=yt;
     // go to 150
   } 

   /*
   // check bins for this panel
  160 continue
   add=half*(y(i)+yl)*(x(i)-xl)
   if (x(i) == xl) go to 250
   xil=1/(x(i)-xl)
   if (i == 1 and j == nbin-1) go to 165
   if (sum+add >= fract*shade and j < nbin-1) go to 170
   sum=sum+add
   gral=gral+half*(yl*x(i)-y(i)*xl)*(x(i)+xl)&
     +third*(y(i)-yl)*(x(i)**2+x(i)*xl+xl**2)
   go to 250
  165 continue
   xn=x(i)
   j=j+1
   go to 190
  170 continue
   j=j+1
   if (yl < sigmin) go to 175
   test=(fract-sum)*(y(i)-yl)/((x(i)-xl)*yl**2)
   if (abs(test) > ytol) go to 175
   xn=xl+(fract-sum)/yl
   if (xn > x(i)) go to 180
   if (xn >= xl and xn <= x(i)) go to 190
   call error('sigl','no legal solution.',' ')
  175 continue
   f=(y(i)-yl)*xil
   rf=1/f
   disc=(yl*rf)**2+2*(fract-sum)*rf
   if (disc < zero) {
      write(strng,'(''disc='',1p,e12.4)') disc
      call mess('sigl',strng,'set to abs value and continue')
      disc=abs(disc)
   } // end if
   if (f > zero) xn=xl-(yl*rf)+sqrt(disc)
   if (f < zero) xn=xl-(yl*rf)-sqrt(disc)
   if (xn > xl and xn <= x(i)) go to 190
   if (xn > xl and xn < (x(i)+ytol*(x(i)-xl))) go to 180
   call error('sigl','no legal solution (quadratic path).',' ')
  180 continue
   xn=x(i)
  190 continue
   yn=yl+(y(i)-yl)*(xn-xl)*xil
   gral=gral+(xn-xl)*(yl*half*(xn+xl)&
     +(y(i)-yl)*xil*(-xl*half*(xn+xl)&
     +third*(xn**2+xn*xl+xl**2)))
   xbar=gral*rfract

   // compute legendre components
   if (nlin >= 0) {
      call legndr(xbar,p,nl)
      do il=2,nl
         s(il)=s(il)+p(il)*rnbin
      } // end do

   // output equally probable angles
   else {
      s(j+1)=xbar
   } // end if

   // continue bin loop and linearization loop
   xl=xn
   yl=yn
   sum=0
   gral=0
   if (j == nbin) go to 260
   if (xl < x(i)) go to 160
  250 continue
   xl=x(i)
   yl=y(i)
   i=i-1
   if (i > 1) go to 150
   if (i == 1) go to 160
  260 continue
   return
   end subroutine sigl

   */


}
       //counter += 1; if (counter > 20 ){std::cout << "return from counter" << std::endl; return;}
