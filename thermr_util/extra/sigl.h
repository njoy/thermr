#include <iostream>
#include <vector>
#include "sig.h"

auto sigl( double e, double ep, double tev, std::vector<double> alpha,
  std::vector<double> beta, std::vector<std::vector<double>> sab, 
  std::vector<double> s, double tolin, int nlin ){

 /*-------------------------------------------------------------------
  * Compute the cross section and legendre components or equally-
  * probable angles for the scattering from e to ep.  Uses linear
  * reconstruction of the angular distribution computed by sig.
  *-------------------------------------------------------------------
  */
  
  int nl,i,j,il;
  double b,seep,sum,xl,yl,ymax,xm,ym,test,test2,rnbin,fract,gral,add,xil,xn,f,rf,disc,yn,xbar,rfract,yt,tol,s1bb
  int imax=20;
  std::vector<double x(imax), y(imax),p(nlin);
  double xtol=.00001e0;
  double ytol=.001e0;
  double sigmin=1.e-32;
  double eps=1.e-3;
  double shade=.99999999e0;

   // constant factors
   b=(ep-e)/tev;
   tol=half*tolin;
   nl=nlin:
   if (nl < 0) nl=-nl;
   b=std::abs(b);
   if (lat == 1 and iinc == 2) b=b*tev/tevz;
   s1bb=sqrt(1+b*b);
   if (ep != 0) seep=1/sqrt(e*ep);

   // adaptive calculation of cross section
   i=3;
   sum=0;
   x[3-1]=-1;
   xl=x[3-1];
   y[3-1]=sig(e,ep,x(3),tev,nalpha,alpha,nbeta,beta,sab);
   yl=y[3-1];
   if (ep == 0) x(2)=0;
   if (ep != 0) x(2)=0.5*(e+ep-(s1bb-1)*az*tev)*seep;
   if (abs(x[2-1]) > 1-eps) x(2)=0.99e0;
   y[2-1]=sig(e,ep,x[2-1],tev,nalpha,alpha,nbeta,beta,sab);
   x[1-1]=1;
   y[1-1]=sig(e,ep,x[1-1],tev,nalpha,alpha,nbeta,beta,sab);
   ymax=y[2-1];
   if (y[1-1] > ymax) ymax=y[1-1];
   if (y[3-1] > ymax) ymax=y[3-1];
   if (ymax < eps) { ymax=eps; }
   /*
  //110 continue
  // if (i == imax) go to 120
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   ym=half*(y(i-1)+y(i))
   yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab)
   test=tol*abs(yt)+tol*ymax/50
   test2=ym+ymax/100
   if (abs(yt-ym).le.test.and.abs(y(i-1)-y(i)).le.test2.and.&
     (x(i-1)-x(i)) < half) go to 120
   if (x(i-1)-x(i) < xtol) go to 120
   i=i+1
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=yt
   go to 110
  120 continue
   sum=sum+half*(y(i)+yl)*(x(i)-xl)
   xl=x(i)
   yl=y(i)
   i=i-1
   if (i > 1) go to 110
   if (i == 1) go to 120
   s(1)=sum
   if (sum > sigmin) go to 130
   do il=1,nl
      s(il)=0
   enddo
   return

   !--prime stack for equally-probable angles
  130 continue
   nbin=nl-1
   rnbin=one/nbin
   fract=sum*rnbin
   rfract=1/fract
   sum=0
   gral=0
   do il=2,nl
      s(il)=0
   enddo
   j=0

   !--adaptive linearization
   i=3
   x(3)=-1
   xl=x(3)
   y(3)=sig(e,ep,x(3),tev,nalpha,alpha,nbeta,beta,sab)
   if (ep == zero) x(2)=0
   if (ep.ne.zero) x(2)=half*(e+ep-(s1bb-1)*az*tev)*seep
   if (abs(x(2)) > 1-eps) x(2)=0.99e0_kr
   x(2)=sigfig(x(2),8,0)
   y(2)=sig(e,ep,x(2),tev,nalpha,alpha,nbeta,beta,sab)
   x(1)=+1
   y(1)=sig(e,ep,x(1),tev,nalpha,alpha,nbeta,beta,sab)
   ymax=y(1)
   if (y(2) > ymax) ymax=y(2)
   if (y(3) > ymax) ymax=y(3)
   if (ymax < eps) ymax=eps
  150 continue
   if (i == imax) go to 160
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,8,0)
   ym=half*(y(i-1)+y(i))
   yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab)
   test=tol*abs(yt)+tol*ymax/50
   test2=ym+ymax/100
   if (abs(yt-ym).le.test.and.abs(y(i-1)-y(i)).le.test2.and.&
     (x(i-1)-x(i)) < half) go to 160
   if (x(i-1)-x(i) < xtol) go to 160
   i=i+1
   x(i)=x(i-1)
   y(i)=y(i-1)
   x(i-1)=xm
   y(i-1)=yt
   go to 150

   !--check bins for this panel
  160 continue
   add=half*(y(i)+yl)*(x(i)-xl)
   if (x(i) == xl) go to 250
   xil=1/(x(i)-xl)
   if (i == 1.and.j == nbin-1) go to 165
   if (sum+add.ge.fract*shade.and.j < nbin-1) go to 170
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
   if (xn.ge.xl.and.xn.le.x(i)) go to 190
   call error('sigl','no legal solution.',' ')
  175 continue
   f=(y(i)-yl)*xil
   rf=1/f
   disc=(yl*rf)**2+2*(fract-sum)*rf
   if (disc < zero) then
      write(strng,'(''disc='',1p,e12.4)') disc
      call mess('sigl',strng,'set to abs value and continue')
      disc=abs(disc)
   endif
   if (f > zero) xn=xl-(yl*rf)+sqrt(disc)
   if (f < zero) xn=xl-(yl*rf)-sqrt(disc)
   if (xn > xl.and.xn.le.x(i)) go to 190
   if (xn > xl.and.xn < (x(i)+ytol*(x(i)-xl))) go to 180
   call error('sigl','no legal solution (quadratic path).',' ')
  180 continue
   xn=x(i)
  190 continue
   yn=yl+(y(i)-yl)*(xn-xl)*xil
   gral=gral+(xn-xl)*(yl*half*(xn+xl)&
     +(y(i)-yl)*xil*(-xl*half*(xn+xl)&
     +third*(xn**2+xn*xl+xl**2)))
   xbar=gral*rfract

   !--compute legendre components
   if (nlin.ge.0) then
      call legndr(xbar,p,nl)
      do il=2,nl
         s(il)=s(il)+p(il)*rnbin
      enddo

   !--output equally probable angles
   else
      s(j+1)=xbar
   endif

   !--continue bin loop and linearization loop
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

