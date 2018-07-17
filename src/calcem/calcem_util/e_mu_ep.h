#include "iel/iel_util/terp1.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"


auto e_mu_ep( int matdp, int mtref, double bk, double t, double& teff, 
  double& teff2, std::vector<double>& scr, double za, double awr, int ncds,
  int nw, int nne, double cliq, int iinc, double emax, std::vector<double> egrid,
  double temp, double breakVal, std::vector<double>& esi ){

   // compute kernel and write in mf6/law7 angle-energy format
   // 510 continue
   int ltt=4;
   int math=matdp;
   int mfh=6;
   int mth=mtref;
   double tev=t*bk;
   teff=teff*bk;
   teff2=teff2*bk;
   scr[1-1]=za;
   scr[2-1]=awr;
   scr[3-1]=0;
   scr[4-1]=ltt; // temporary flag for this format
   scr[5-1]=1;
   scr[6-1]=0;
   //call contio(0,0,nscr,scr,nb,nw);
   ncds=ncds+1;
   scr[1-1]=1;
   scr[2-1]=1;
   scr[3-1]=0;
   scr[4-1]=7;
   scr[5-1]=1;
   scr[6-1]=2;
   scr[7-1]=2;
   scr[8-1]=2;
   scr[9-1]=1.e-5;
   scr[10-1]=1;
   scr[11-1]=emax;
   scr[12-1]=1;
   nw=12;
   //call tab1io(0,0,nscr,scr,nb,nw)
   ncds=ncds+2;
   scr[1-1]=0;
   scr[2-1]=0;
   scr[3-1]=0;
   scr[4-1]=0;
   scr[5-1]=1;
   scr[6-1]=nne;
   scr[7-1]=nne;
   scr[8-1]=2;
   nw=8;
   //call tab2io(0,0,nscr,scr,nb,nw)
   ncds=ncds+2;
   cliq=0;
   if (iinc == 1) std::cout << "go to 515" << std::endl;
   //if (sab(1,1) > sab(2,1)) cliq=(sab(1,1)-sab(1,2))*alpha[1-1]/(beta[2-1]*beta[2-1]);

   // loop over given incident energy grid.
   // the first pass computes the cross section and mu grid
  //515 continue
   int ie=0;
  //520 continue
   ie=ie+1;
   double enow=egrid[ie-1];
   if (ie > 1 and temp > breakVal) enow=enow*temp/breakVal;
   //enow=sigfig(enow,8,0)
   esi[ie-1]=enow;
   int j=0;
   double sum=0;

   /*
   // adaptive reconstruction of angular cross section
   u=-1
   x[2-1]=u
   call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   yy[2-1]=yu[1-1]
   xl=x[2-1]
   yl=yy[2-1]
   u=1
   x[1-1]=u
   call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   yy[1-1]=yu[1-1]
   i=2

   // adaptive reconstruction
 530 continue
   if (i == imax) go to 560
   xm=half*(x(i-1)+x(i))
   xm=sigfig(xm,7,0)
   if (xm <= x(i) or xm >= x(i-1)) go to 560
   call sigu(enow,xm,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   if (x(i-1)-x(i) > .25) go to 575
   ym=yy(i)+(xm-x(i))*(yy(i-1)-yy(i))/(x(i-1)-x(i))
   if (abs(yu[1-1]-ym) > 2*tol*ym+tolmin) go to 575
   ! point passes.  save top point in stack and continue.
  560 continue
   j=j+1
   if (j > mumax-1) call error('calcem','too many angles','see mumax')
   uj(j)=x(i)
   sj(j)=yy(i)
   if (j > 1) then
      sum=sum+half*(yy(i)+yl)*(x(i)-xl)
      xl=x(i)
      yl=yy(i)
   endif
   i=i-1
   if (i >= 2) go to 530
   go to 580
   ! test fails.  add point to stack and continue.
  575 continue
   i=i+1
   x(i)=x(i-1)
   x(i-1)=xm
   yy(i)=yy(i-1)
   yy(i-1)=yu[1-1]
   go to 530
   ! linearization complete.  write out result.
  580 continue
   j=j+1
   uj(j)=x[1-1]
   sj(j)=yy[1-1]
   nmu=j
   ubar[ie-1]=0
   sum=sum+half*(yy[1-1]+yl)*(x[1-1]-xl)
   xsi[ie-1]=sum/2
   do i=2,nmu
      ubar[ie-1]=ubar[ie-1]+half*(uj(i)-uj(i-1))*(sj(i)+sj(i-1))*(uj(i)+uj(i-1))
   enddo
   ubar[ie-1]=half*ubar[ie-1]/sum
   //if (iprint == 2) then
   //   write(nsyso,'(/i5,'' enow '',1p,e13.6,''   xsec '',e13.6,&
   //     &''   mubar  '',e13.6)') ie,enow,xsi[ie-1],ubar[ie-1]
   //   write(nsyso,'(''      num of mu '', i5)') nmu
   //   write(nsyso,'(/''            mu            theta      dsigma/dmu'')')
   //   do i=1,nmu
   //      write(nsyso,'(i5,1x,f15.8,1x,f12.4,1x,1p,e14.7)') i,uj(i),&
   //            acos(uj(i))*180.0/3.14159265359,sj(i)/2.0
   //   enddo
   //endif

   // now loop through the mu grid to write out the distributions
   mth=mtref
   scr[1-1]=0
   scr[2-1]=enow
   scr[3-1]=0
   scr[4-1]=0
   scr[5-1]=1
   scr[6-1]=nmu
   scr[7-1]=nmu
   scr[8-1]=2
   nw=8
   call tab2io(0,0,nscr,scr,nb,nw)
   ncds=ncds+1
   do il=1,nmu
      u=uj(il)
      call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
      nep=nint(yu[2-1])
      j=0
      do i=1,nep
        j=nep-i
        if (yu(2*(nep-i)+4)/sum > yumin) exit
      enddo
      nep=j
     if (iprint == 2) then
         write(nsyso,'(/'' mu = '',f15.8)') u
         write(nsyso,'('' (e-prime, pdf);  num of e-prime '', i5)') nep
         write(nsyso,*)
         // test yu()/sum below: is this pdf normalized to 1.0 ?
         write(nsyso,'(1p,3(1x,e14.7,1x,e14.7,1x))')&
               (yu(2*i+1),yu(2*i+2)/sum,i=1,nep)
      endif
      scr[1-1]=0
      scr[2-1]=u
      scr[3-1]=0
      scr[4-1]=0
      scr[5-1]=1
      scr[6-1]=nep
      scr[7-1]=nep
      scr[8-1]=2
      k=8;
      istart=1
     595 continue
      iend=nep
      if ((iend-istart) >= npage/2) iend=istart+npage/2-1
      j=k-1
      ib=istart-1
     596 continue
      j=j+2
      ib=ib+1
      scr(j)=yu(1+2*ib)
      scr(j+1)=yu(2+2*ib)*2/sum
      if (ib < iend) go to 596
      nw=j+1
      if (k == 0) go to 597
      k=0
      call tab1io(0,0,nscr,scr,nb,nw)
      if (nb == 0) go to 598
      istart=iend+1
      go to 595
    597 continue
      call moreio(0,0,nscr,scr,nb,nw)
      if (nb == 0) go to 598
      istart=iend+1
      go to 595
     598 continue
      ncds=ncds+1+(j*(nep+1)+5)/6
   enddo

   if (ie < nne) go to 520


     */

} 
