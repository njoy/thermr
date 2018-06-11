#include "../extra/terp1.h"
#include "../extra/legndr.h"
#include "calcem_util/sigl.h"
#include "calcem_util/sigfig.h"


void e_mu_ep( int ltt, int math, int matdp, int mfh, int mtref, double& tev,
  double& teff, double& teff2, const double& bk, 
  
  std::vector<double>& scr, double& za, 
  double& awr, int& ncds, double& emax, double& cliq, int iinc, int lat, 
  std::vector<double>& esi, std::vector<double>& xsi, const int& lasym, 
  std::vector<double>& alpha, std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, const double t, const double& tol, 
  const double& az, const double& az2, const double& sb, const double& sb2, 
  int& nnl, const int& nl, const int& jmax, const int& nne, int iprint,
  double& enow, double& break_val, double& temp, std::vector<double>& x, 
  int inew, int iold, int imax ){

  int nalpha = alpha.size(), nbeta = beta.size(), i, j, mumax, ie, nmu, nep, k, istart, iend, isave, npage, ib, nw, nexn, nb, ien, nex, ne, mth;
  double sum, u, xl, yl, xm, ym, tolmin, xs, yumin;
  std::vector<double> yy, yu, sj, uj, ubar, ex,egrid;


   // compute kernel and write in mf6/law7 angle-energy format
   // 510 continue
   ltt=4;
   math=matdp;
   mfh=6;
   mth=mtref;
   tev=t*bk;
   teff=teff*bk;
   teff2=teff2*bk;
   scr[1-1]=za;
   scr[2-1]=awr;
   scr[3-1]=0;
   scr[4-1]=ltt; // temporary flag for this format
   scr[5-1]=1;
   scr[6-1]=0;
   // call contio(0,0,nscr,scr,nb,nw)
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
   // call tab1io(0,0,nscr,scr,nb,nw)
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
   // call tab2io(0,0,nscr,scr,nb,nw)
   ncds=ncds+2;
   cliq=0;
   if (iinc == 1) std::cout << "go to 515" << std::endl;
   if (sab[1-1][1-1] > sab[2-1][1-1]) cliq=(sab[1-1][1-1]-sab[1-1][2-1])*alpha[1-1]/(beta[2-1]*beta[2-1]);

   // loop over given incident energy grid.
   // the first pass computes the cross section and mu grid
  
   // 515 continue
   ie=0;

  // 520 continue
   ie=ie+1;
   enow=egrid[ie-1];
   if (ie > 1 and temp > break_val) enow=enow*temp/break_val;
   enow=sigfig(enow,8,0);
   j=0;
   sum=0;

   // adaptive reconstruction of angular cross section
   u=-1;
   x[2-1]=u;
   // call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol);
   yy[2-1]=yu[1-1];
   xl=x[2-1];
   yl=yy[2-1];
   u=1;
   x[1-1]=u;
   // call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   yy[1-1]=yu[1-1];
   i=2;

   // adaptive reconstruction
 // 530 continue
   if (i == imax) std::cout << "go to 560" << std::endl;
   xm=0.5*(x[i-1-1]+x[i-1]);
   xm=sigfig(xm,7,0);
   if (xm <= x[i-1] or xm >= x[i-1-1]) std::cout << "go to 560" << std::endl;
   // call sigu(enow,xm,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
   if (x[i-1-1]-x[i-1] > 0.25) std::cout << "go to 575" << std::endl;
   ym=yy[i-1]+(xm-x[i-1])*(yy[i-1-1]-yy[i-1])/(x[i-1-1]-x[i-1]);
   if (abs(yu[1-1]-ym) > 2*tol*ym+tolmin) std::cout << "go to 575" << std::endl;
   // point passes.  save top point in stack and continue.
   // 560 continue
   j=j+1;
   if (j > mumax-1) std::cout << "call error('calcem','too many angles','see mumax')" << std::endl;
   uj[j-1]=x[i-1];
   sj[j-1]=yy[i-1];
   if (j > 1) {
      sum=sum+0.5*(yy[i-1]+yl)*(x[i-1]-xl);
      xl=x[i-1];
      yl=yy[i-1];
   } // endif
   i=i-1;
   if (i >= 2) std::cout << "go to 530" << std::endl;
   std::cout << "go to 580" << std::endl;
   // test fails.  add point to stack and continue.
   // 575 continue
   i=i+1;
   x[i-1]=x[i-1-1];
   x[i-1-1]=xm;
   yy[i-1]=yy[i-1-1];
   yy[i-1-1]=yu[1-1];
   std::cout << "go to 530" << std::endl;
   // linearization complete.  write out result.
   // 580 continue
   j=j+1;
   uj[j-1]=x[1-1];
   sj[j-1]=yy[1-1];
   nmu=j;
   ubar[ie-1]=0;
   sum=sum+0.5*(yy[1-1]+yl)*(x[1-1]-xl);
   xsi[ie-1]=sum/2;
   for ( int i = 1; i < nmu; ++i ){
      ubar[ie-1]=ubar[ie-1]+0.5*(uj[i]-uj[i-1])*(sj[i]+sj[i-1])*(uj[i]+uj[i-1]);
   } // end do
   ubar[ie-1]=0.5*ubar[ie-1]/sum;
   if (iprint == 2) {
      //write(nsyso,'(/i5,'' enow '',1p,e13.6,''   xsec '',e13.6,&
      //  &''   mubar  '',e13.6)') ie,enow,xsi(ie),ubar(ie)
      //write(nsyso,'(''      num of mu '', i5)') nmu
      //write(nsyso,'(/''            mu            theta      dsigma/dmu'')')
      for ( int i = 0; i < nmu; ++i ){
         //write(nsyso,'(i5,1x,f15.8,1x,f12.4,1x,1p,e14.7)') i,uj(i),&
         //      acos(uj(i))*180.0/3.14159265359,sj(i)/2.0
      } // end do
   } // endif

   // now loop through the mu grid to write out the distributions
   mth=mtref;
   scr[1-1]=0;
   scr[2-1]=enow;
   scr[3-1]=0;
   scr[4-1]=0;
   scr[5-1]=1;
   scr[6-1]=nmu;
   scr[7-1]=nmu;
   scr[8-1]=2;
   nw=8;
   // call tab2io(0,0,nscr,scr,nb,nw)
   ncds=ncds+1;
   for ( int il = 0; il < nmu; ++il ){
      u=uj[il];
      //call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
      nep=int(yu[2-1]);
      j=0;
      for ( int i = 1; i <= nep; ++i ){
        j=nep-i;
        if (yu[2*(nep-i)+4-1]/sum > yumin) return; //exit
      } // end do
      nep=j;
     if (iprint == 2) {
         //write(nsyso,'(/'' mu = '',f15.8)') u
         //write(nsyso,'('' (e-prime, pdf);  num of e-prime '', i5)') nep
         //write(nsyso,*)
         // test yu()/sum below: is this pdf normalized to 1.0 ?
         //write(nsyso,'(1p,3(1x,e14.7,1x,e14.7,1x))')&
         //      (yu(2*i+1),yu(2*i+2)/sum,i=1,nep)
      } // endif
      scr[1-1]=0;
      scr[2-1]=u;
      scr[3-1]=0;
      scr[4-1]=0;
      scr[5-1]=1;
      scr[6-1]=nep;
      scr[7-1]=nep;
      scr[8-1]=2;
      k=8;
      istart=1;
     // 595 continue
      iend=nep;
      if ((iend-istart) >= npage/2) iend=istart+npage/2-1;
      j=k-1;
      ib=istart-1;
     // 596 continue
      j=j+2;
      ib=ib+1;
      scr[j-1]=yu[1+2*ib-1];
      scr[j+1-1]=yu[2+2*ib-1]*2/sum;
      if (ib < iend) std::cout << "go to 596" << std::endl;
      nw=j+1;
      if (k == 0) std::cout << "go to 597" << std::endl;
      k=0;
      // call tab1io(0,0,nscr,scr,nb,nw)
      if (nb == 0) std::cout << "go to 598" << std::endl;
      istart=iend+1;
      std::cout << "go to 595" << std::endl;
    // 597 continue
      //call moreio(0,0,nscr,scr,nb,nw)
      if (nb == 0) std::cout << "go to 598" << std::endl;
      istart=iend+1;
      // go to 595
     // 598 continue
      ncds=ncds+1+(j*(nep+1)+5)/6;
   } // end do

   if (ie < nne) std::cout << "go to 520" << std::endl;

   // calculate incoherent inelastic scattering cross sections
   
   // 610 continue
   for (int ie = 0; ie < 20; ++ie ){ 
      ex[ie]=0;
   } // end do
   nexn=nex;
   if (nex == 2) nexn=3;
   ie=0;
   ien=1;
   while (ien > 0){
      ie=ie+1;
      // call finda(ie,ex,nex,iold,bufo,nbuf)
      if (iinc == 1) {
         // use elastic cross section  for free atom scattering.
         ex[3-1]=ex[2-1];
      }
      else {
         // use computed cross section for bound atoms.
         enow=ex[1-1];
         //xs=terp(esi,xsi,nne,enow,nlt)
         if (ie == ne) xs=0;
         ex[3-1]=xs;
      } // endif
      ien=ie;
      if (ie == ne) ien=-ien;
      // call loada(ien,ex,nexn,inew,bufn,nbuf)
   } // end do
   esi[nne+1-1]=0;
   isave=iold;
   iold=inew;
   inew=isave;
   nex=nexn;
   if (ie == ne) xs=0;

   // calcem is finished.
   // call asend(0,nscr)
   // deallocate(alpha)
   // deallocate(beta)
   // deallocate(sab)
   return;
  
}
