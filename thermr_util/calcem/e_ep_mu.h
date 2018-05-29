
#include "../extra/terp1.h"
#include "calcem_util/sigl.h"
#include "calcem_util/sigfig.h"


auto do313( int& jbeta, const int& lat, double& ep, double& enow, 
  const std::vector<double>& beta, const double& tev, const double& tevz,
  const std::vector<double>& x, int& iskip){
  bool do_313 = true;
  std::cout << 313 << std::endl;
  while ( do_313 ){
    //std::cout << std::setprecision(20) << "E'     " << ep << std::endl;
    // 313 continue
    if (jbeta == 0) jbeta=1;
    if (jbeta <= 0) {
      if (lat == 1) {
   //     std::cout << enow << "    " << beta[-jbeta-1] << "    " << tevz << "    " << jbeta << std::endl;
        ep=enow-beta[-jbeta-1]*tevz;
      }
      else {
        ep=enow-beta[-jbeta-1]*tev;
      } // endif
      if (ep == enow) {
         ep=sigfig(enow,8,-1);
      }
      else {
         ep=sigfig(ep,8,0);
      }

    } 
    else {
      if (lat == 1) {
        ep=enow+beta[jbeta-1]*tevz;
      }
      else {
        ep=enow+beta[jbeta-1]*tev;
      } // endif
      if (ep == enow) {
         ep=sigfig(enow,8,+1);
         iskip=1;
      }
      else { 
         ep=sigfig(ep,8,0);
      }

    } // endif
    if (ep > x[2-1]) { do_313 = false; } // go to 316
    jbeta=jbeta+1;
    // go to 313
  }

  //std::cout << std::setprecision(20) << "E'     " << ep << std::endl;
}





auto do310( int& ie, double& enow, std::vector<double>& egrid, const double& temp,
  const double& bk, const double& break_val, const double& therm, 
  std::vector<double>& esi, std::vector<double>& xsi, std::vector<double>& ubar,
  std::vector<double>& p2, std::vector<double>& p3, double& ep, int& jbeta,
  int& nbeta, int& iskip, int& j, std::vector<std::vector<double>>& y, 
  std::vector<double>& yt, const int nl, const int& lasym, std::vector<double>& x,
  const int& ngrid){
     std::cout << 310 << std::endl;
     double tone, elo;
     ie=ie+1;
     enow=egrid[ie-1];
     if (temp > break_val) {
       tone=therm/bk;
       elo=egrid[0];
       enow=elo*exp(log(enow/elo)*log((temp/tone)*egrid[ngrid-1]/elo)/log(egrid[ngrid-1]/elo));
     } // endif
     esi[ie-1]=enow;
     xsi[ie-1]=0;
     ubar[ie-1]=0;
     p2[ie-1]=0;
     p3[ie-1]=0;
     ep=0.0;
     x[0]=ep;
     //call sigl(enow,ep,nnl,tev,nalpha,alpha,nbeta,beta,&
     //  sab,yt,nlmax,tol)
     for ( int il = 0; il < yt.size(); ++il ){
        y[il][0] = yt[il];
     } // enddo
     jbeta=-nbeta;
     if (lasym > 0) jbeta=1;
     j=0;
     iskip=0;
}

auto do360(int& j, const int& jmax, std::vector<double>& xsi, 
    std::vector<double>& x, double& xlast, double& ylast, double& ulast,
    double& u2last, double& u3last, const double& tolmin,
    std::vector<std::vector<double>>& y, std::vector<double>& p2,
    std::vector<double>& p3,
    double& uu, double& u2, double& u3, int& nll, const int& nl, 
    std::vector<double>& p, std::vector<double>& ubar, const int ie, 
    const int i ){
  // 360 continue
  std::cout << 360 << std::endl;
  j=j+1;
  if (j >= jmax) std::cout << "call error('calcem','storage exceeded.',' ')" << std::endl;
  if (j > 1) xsi[ie-1]=xsi[ie-1]+(x[i-1]-xlast)*(y[0][i-1]+ylast)*0.5;
  if (j > 1) {
    uu=0;
    u2=0;
    u3=0;
    nll=3;
    for ( int il = 1; il < nl; ++il ){
      //call legndr(y(il,i),p,nll)
      uu=uu+p[2-1];
      u2=u2+p[3-1];
      u3=u3+p[4-1];
    } // enddo
    uu=uu/(nl-1);
    uu=uu*y[0][i-1];
    u2=u2/(nl-1);
    u3=u3/(nl-1);
    u2=u2*y[0][i-1];
    u3=u3*y[0][i-1];
    ubar[ie-1]=ubar[ie-1]+0.5*(x[i-1]-xlast)*(uu+ulast);
    p2[ie-1]=p2[ie-1]+0.5*(x[i-1]-xlast)*(u2+u2last);
    p3[ie-1]=p3[ie-1]+0.5*(x[i-1]-xlast)*(u3+u3last);
  } // endif
  // if (j != 3 or xsi[ie-1] >= tolmin) go to 380
  if (j == 3 and xsi[ie-1] < tolmin){ j = 2; }
  // j=2

}



auto e_ep_mu(int& math, int& matdp, double& teff, double& teff2, 
  std::vector<double>& scr, const int& mtref, double& za, double& awr,
  int& ncds, double& emax, double& cliq, int iinc, int lat,
  std::vector<double>& esi, std::vector<double>& xsi, const int& lasym,
  std::vector<double>& alpha, std::vector<double>& beta, 
  std::vector<std::vector<double>>& sab, const double t, const double& tol,
  const double& az, const double& az2, const double& sb, const double& sb2, 
  int& nnl ){

  int itemp,iold,inew,ne,nex;
  double temp;
  int nr,np,nwtab,nl,nlt,nlp,nlp1,jmax,nne;
  int i,ia,nl1,ltt,loc,l,jscr,ilog;
  int matd,itprnt,nb,nw,ni,nbeta=beta.size(),lt,it,nalpha=alpha.size();
  int itrunc,ib,ip,ir,idis,ie,nmu,nep,istart,iend;
  int jbeta=0,j=0,iskip,il,k,jnz,ll,jj,ii,nexn,ien,isave;
  int nll;
  double smz,t1,tmax,test,tempt,tt1,ttn,tnxt,xm;
  double uu,uum,ym,test2,xlast,ulast,xs;
  double b,diff,enow,ep,sabmin,tev,tevz,ylast;
  double u,xl,yl,sum;
  double tone,elo;
  const int ngrid=118;
  const int nlmax=65;
  const int nemax=5000;
  const int mumax=300;
  const int imax=20;
  std::vector<double> ex(imax,0.0),x(imax,0.0),yt(nlmax,0.0);
  std::vector<std::vector<double>> y(nlmax,std::vector<double> (imax,0.0));
  std::vector<double> yy(imax,0.0),yu(2*nemax,0.0),ubar(ngrid,0.0),p2(ngrid,0.0),p3(ngrid,0.0),p(4,0.0),uj(mumax,0.0),sj(mumax,0.0);
  double u2=0,u2last=0,u3=0,u3last=0;
  std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
    1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 0.000253, 0.000297, 0.000350, 0.00042, 
    0.000506, 0.000615, 0.00075, 0.00087, 0.001012, 0.00123, 0.0015, 0.0018, 
    0.00203, 0.002277, 0.0026, 0.003, 0.0035, 0.004048, 0.0045, 0.005, 0.0056, 
    0.006325, 0.0072, 0.0081, 0.009108, 0.01, 0.01063, 0.0115, 0.012397, 
    0.0133, 0.01417, 0.015, 0.016192, 0.0182, 0.0199, 0.020493, 0.0215, 0.0228, 
    0.0253, 0.028, 0.030613, 0.0338, 0.0365, 0.0395, 0.042757, 0.0465, 0.050, 
    0.056925, 0.0625, 0.069, 0.075, 0.081972, 0.09, 0.096, 0.1035, 0.111573, 
    0.120, 0.128, 0.1355, 0.145728, 0.160, 0.172, 0.184437, 0.20, 0.2277, 
    0.2510392, 0.2705304, 0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 
    0.4170351, 0.45, 0.5032575, 0.56, 0.625, 0.70, 0.78, 0.86, 0.95, 1.05, 
    1.16, 1.28, 1.42, 1.55, 1.70, 1.855, 2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 
    3.42, 3.75, 4.07, 4.46, 4.90, 5.35, 5.85, 6.40, 7.00, 7.65, 8.40, 9.15, 
    9.85, 10.00 };
  const double sabflg = -225, eps = 1.e-4, tolmin = 5.e-7, therm = 0.0253, 
    break_val = 3000, em9 = 1e-9, up = 1.1, dn = 0.9, uumin = 0.00001, yumin = 2e-7, 
    nlpmx = 10, bk =8.617385e-5;
  int mth, mfh;
  // save nwtab,sabmin,nl,nlt,nlp,nlp1,nl1,nnl,jmax,nne
  tevz = therm;


  // compute kernel and write in special mf6 energy-angle format
  // 300 continue
  std::cout << "300" << std::endl;
  // if (iform == 1) go to 510 
  // ^^ Take care of this in the upstairs function
  // This will go to the other branch ( the E-mu-E' branch )
  ltt=5;
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
  nw=6;
  // call contio(0,0,nscr,scr,nb,nw)
  ncds=ncds+1;
  scr[1-1]=1;
  scr[2-1]=1;
  scr[3-1]=-1; // LIP=-1 indicates incoherent inelastic data
  scr[4-1]=1; // LAW=1 for incoherent inelastic data
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
  scr[1-1]=temp;
  scr[2-1]=0;
  scr[3-1]=3; // lang=3 is a special code for equally probable cosines
  scr[4-1]=1;
  scr[5-1]=1;
  scr[6-1]=nne;
  scr[7-1]=nne;
  scr[8-1]=2;
  nw=8;
  // call tab2io(0,0,nscr,scr,nb,nwtab)
  ncds=ncds+2;
  cliq=0;

  if ( iinc != 1 and sab[0][0] > sab[0][1]){
    cliq=(sab[0][0]-sab[0][1])*alpha[0]/(beta[1]*beta[1]);
  }

  // loop over given incident energy grid.
  // 305 continue
  std::cout << 305 << std::endl;
  ie=0;

  // Loop over incident energy (green line in my drawing)
  bool loopE = true;
  while (loopE){

    // 310 continue
    do310( ie, enow, egrid, temp, bk, break_val, therm, esi, xsi, ubar, p2, 
      p3, ep, jbeta, nbeta, iskip, j, y, yt, nl, lasym, x, ngrid );


    bool moreBeta = true;
    while (moreBeta) {
      // set up next panel
      // 311 continue
      std::cout << 311 << std::endl;
      x[2-1]=x[1-1];

      for (int il = 0; il < y.size(); ++il ){
        y[il][2-1] = y[il][1-1];
      } // enddo


      // do 313
      do313( jbeta, lat, ep, enow, beta, tev, tevz, x, iskip);



      std::cout << 316 << std::endl;
      x[1-1]=ep;


      std::cout << "\n\n\n\n\n\n\n\n\n" << std::endl;

      /*
      std::cout << nnl << std::endl;
      std::cout << nlmax << std::endl;
      std::cout << std::setprecision(12) << enow << std::endl;
      std::cout << ep << std::endl;
      std::cout << std::endl;

      std::cout << tev << std::endl;
      std::cout << tevz << std::endl;
      std::cout << az << std::endl;
      std::cout << az2 << std::endl;
      std::cout << teff << std::endl;
      std::cout << teff2 << std::endl;
      */

      std::cout << tol << std::endl;
      std::cout << iinc << std::endl;
      std::cout << lat << std::endl;
      std::cout << lasym << std::endl;
      std::cout << cliq << std::endl;
      std::cout << sb << std::endl;
      std::cout << sb2 << std::endl;
      

      sigl( nnl, nlmax, enow, ep, tev, alpha, beta, sab, yt, tol, az, 
        tevz, iinc, lat, lasym, az2, teff2, cliq, sb, sb2, teff );

      for (int il = 0; il < yt.size(); ++il ){
        y[il][1-1]=yt[il];
      } // enddo
        std::cout << std::setw(20) << yt[0] << std::setw(20) << yt[1] << std::setw(20) << yt[2] << std::setw(20) << yt[3] << std::setw(20) << yt[4] << std::endl;
      //  std::cout << std::setw(20) << x[0] << std::setw(20) << x[1] << std::setw(20) << x[2] << std::setw(20) << x[3] << std::setw(20) << x[4] << std::endl;
      //std::cout << std::setw(20) << y[0][0] << std::setw(20) << y[0][1] << std::setw(20) << y[0][2] << std::setw(20) << y[0][3] << std::setw(20) << y[0][4] << std::endl;

      return;


      // adaptive subdivision of panel
      i=2;
      // compare linear approximation to true function
  
      bool passedTest = false;
      while (not passedTest){
        std::cout << 330 << std::endl;
        // 330 continue



        if ( i != imax ){
          if ( iskip != 1 ){
            std::cout << i << std::endl;
            std::cout << std::setw(10) << x[i-2] << std::setw(35) << x[i-1] << std::endl;
            std::cout << std::setw(10) << y[0][i-2] << std::setw(35) << y[0][i-1] << std::endl;
            if (0.5*(y[0][i-1-1]+y[0][i-1])*(x[i-1-1]-x[i-1]) >= tolmin) {
              xm=0.5*(x[i-1-1]+x[i-1]);

              if (xm > x[i-1] and xm < x[i-1-1]) {

                std::cout << ":)   " << yt[0] << "    " << yt[1] << "    " << yt[2]<< std::endl;
                sigl( nnl, nlmax, enow, xm, tev, alpha, beta, sab, yt, tol, az, 
                    tevz, iinc, lat, lasym, az2, teff2, cliq, sb, sb2, teff );

                std::cout << ":)   " << yt[0] << "    " << yt[1] << std::endl;

                   std::cout << "HERE" << std::endl; return;
                uu=0;
                uum=0;
                for ( int k = 1; k <= nl; ++k ){
                  ym = terp1(x[i-1],y[k-1][i-1],x[i-1-1],y[k-1][i-1-1],xm,2);

                  if (k > 1) uu=uu+yt[k-1];
                  if (k > 1) uum=uum+ym;
                  test=tol*abs(yt[k-1]);
                  test2=test;
                  if (k > 1) test2=tol;
                   std::cout << yt[k-1] << "    "  << ym << "    " << test2<< std::endl; return;
                  if (abs(yt[k-1]-ym) > test2) std::cout << "go to 410" << std::endl;
                } 
                std::cout << 350 << std::endl;  // 350 continue
                test=2*tol*abs(uu-1)+uumin;
                if (abs(uu-uum) > test) std::cout << "go to 410" << std::endl;
                // point passes.  save top point in stack and continue.

              } 
            } 
          } // iskip
          else {
            iskip = 0;
          }
        } 
      
        do360(j, jmax, xsi, x, xlast, ylast, ulast, u2last, u3last, tolmin, y, 
          p2, p3, uu, u2, u3, nll, nl, p, ubar, ie, i );


     return;

      } // passedTest (going back to 330)
    } // moreBeta (311)
  } //  loopE (310)

         /*
          380 continue
           jscr=7+(j-1)*(nl+1)
           scr(jscr)=x(i)
           if (y(1,i) >= em9) {
              scr(1+jscr)=sigfig(y(1,i),9,0)
           else
              scr(1+jscr)=sigfig(y(1,i),8,0)
           } // endif
           do il=2,nl
              scr(il+jscr)=sigfig(y(il,i),9,0)
              if (scr(il+jscr) > 1) {
                 // only warn for big miss, but always fix the overflow
                 //  use this same 1+0.0005 value in aceth
                 if (scr(il+jscr) > 1+0.0005) {
                    write(strng,'("1cos=",f7.4,", set to 1.&
                                  &  enow,e''=",2(1pe12.5))')&
                                    scr(il+jscr),enow,scr(jscr)
                    call mess('calcem',strng,'')
                 } // endif
                 scr(il+jscr)=1
              } // endif
              if (scr(il+jscr) < -1) {
                 // only warn for big miss, but always fix the underflow
                 if (scr(il+jscr) < -(1+0.0005)) {
                    write(strng,'("1cos=",f7.4,", set to -1.&
                                  &  enow,e''=",2(1pe12.5),i3)')&
                                    scr(il+jscr),enow,scr(jscr)
                    call mess('calcem',strng,'')
                 } // endif
                 scr(il+jscr)=-1
              } // endif
           } // enddo
           xlast=x(i)
           ylast=y(1,i)
           if (ylast != 0) jnz=j
           ulast=0
           u2last=0
           u3last=0
           nll=3
           do il=2,nl
              call legndr(y(il,i),p,nll)
              ulast=ulast+p(2)
              u2last=u2last+p(3)
              u3last=u3last+p(4)
           } // enddo
           ulast=ulast*y(1,i)/(nl-1)
           u2last=u2last*y(1,i)/(nl-1)
           u3last=u3last*y(1,i)/(nl-1)
           i=i-1
           // if (i >= 2) go to 330
           if (i >= 2) continue;
           jbeta=jbeta+1
           if (jbeta <= nbeta) go to 311
           do il=1,nl
              y(il,i)=0
           } // enddo
           go to 430
           // test fails.  add point to stack and continue.
      */
 //     } // Loops through all times that 
      /*
          410 continue
           i=i+1
           x(i)=x(i-1)
           x(i-1)=xm
           do il=1,nl
              y(il,i)=y(il,i-1)
              y(il,i-1)=yt(il)
           } // enddo
           go to 330
           // linearization complete.  write out result.
          430 continue
           j=j+1
           xsi(ie)=xsi(ie)+(x(i)-xlast)*(y(1,i)+ylast)*0.5
           uu=0
           u2=0
           u3=0
           ubar(ie)=ubar(ie)+0.5*(x(i)-xlast)*(uu+ulast)
           p2(ie)=p2(ie)+0.5*(x(i)-xlast)*(u2+u2last)
           p3(ie)=p3(ie)+0.5*(x(i)-xlast)*(u3+u3last)
           xsi(ie)=sigfig(xsi(ie),9,0)
           scr(7+(nl+1)*(j-1))=x(i)
           jscr=7+(j-1)*(nl+1)
           if (y(1,i) >= em9) {
              scr(1+jscr)=sigfig(y(1,i),9,0)
           else
              scr(1+jscr)=sigfig(y(1,i),8,0)
           } // endif
           do il=2,nl
              scr(il+jscr)=sigfig(y(il,i),9,0)
              if (scr(il+jscr) > 1) {
                 // only warn for big miss, but always fix the overflow
                if (scr(il+jscr) > 1+0.0005) {
                  write(strng,'("2cos",f7.4,", set to 1.&
                          &  enow,e''=",2(1pe12.5))')&
                            scr(il+jscr),enow,scr(jscr)
            call mess('calcem',strng,'')
         } // endif
         scr(il+jscr)=1
      } // endif
      if (scr(il+jscr) < -1) {
         // only warn for big miss, but always fix the underflow
         if (scr(il+jscr) < -(1+0.0005)) {
            write(strng,'("2cos",f7.4,", set to -1.&
                          &  enow,e''=",2(1pe12.5),i3)')&
                            scr(il+jscr),enow,scr(jscr)
            call mess('calcem',strng,'')
         } // endif
         scr(il+jscr)=-1
      } // endif
   } // enddo
   if (y(1,1) != 0) jnz=j
   if (jnz < j) j=jnz+1
   if (iprint == 2) {
      ubar(ie)=ubar(ie)/xsi(ie)
      p2(ie)=p2(ie)/xsi(ie)
      p3(ie)=p3(ie)/xsi(ie)
      write(nsyso,'(/,1x,"incident energy =",1pe13.6,&
                   &  5x,"cross section =",1pe13.6,&
                   &  5x,"mubar,p2,p3 =",3(1pe12.4))')&
                          enow,xsi(ie),ubar(ie),p2(ie),p3(ie)
      write(nsyso,'(/,5x,"exit energy",11x,"pdf",7x,"cosines")')
      write(nsyso,'(  3x,"---------------",5x,"-----------",2x,88("-"))')
      ll=6
      do jj=1,j
         write(nsyso,'(2x,1pe15.8,5x,1pe12.5,0p,8f11.6)')&
           (scr(ll+ii),ii=1,nlp)
         if (nl1 > nlp) write(nsyso,'(34x,8f11.6)')&
           (scr(ll+ii),ii=nlp1,nl1)
         ll=ll+nl1
      } // enddo
   } // endif
   scr(1)=0
   scr(2)=enow
   scr(3)=0
   scr(4)=0
   scr(5)=(nl+1)*j
   scr(6)=nl+1
   ncds=ncds+1+(j*(nl+1)+5)/6
   call listio(0,0,nscr,scr,nb,nw)
   loc=1
   do while (nb != 0)
      loc=loc+nw
      call moreio(0,0,nscr,scr(loc),nb,nw)
   } // enddo
   if (ie < nne) go to 310
   go to 610

   */
//   } // While looking over incident energies E

}
