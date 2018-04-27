#include "sig.h"
#include "../../coh_util/sigcoh_util/legndr.h"


auto do_110_120_130(int& i, int& imax, double& xm, double& ym, std::vector<double>& x,
  std::vector<double>& y, double& e, double& ep, double& tev, double& tevz, 
  int& nalpha, int& nbeta, std::vector<double>& alpha,std::vector<double>& beta,
  std::vector<std::vector<double>>& sab, double& bbm, double& az, double& az2,
  int& lasym, double& teff, double& teff2, int& lat, double& cliq, double& sb,
  double& sb2, int& iinc, double& yt, double& sum, int& nl, double& sigmin,
  std::vector<double>& s, int& nbin, double& fract, double& rnbin, double& xl,
  int& j, double& ymax, double& eps, double& rfract, double& seep, double& gral, 
  double& yl, double& test2, double& s1bb, double& test, double& tol, double& xtol){

  bool do110 = true;
  bool do110_inner = true;
  bool do120 = true;
  while ( do110 ){
    while ( do110_inner ){ 
      std::cout << "110" << std::endl;
      // 110 continue
      
      // if (i == imax) go to 120
      if (i == imax) break;

      xm=0.5*(x[i-2]+x[i-1]);
      ym=0.5*(y[i-2]+y[i-1]);

      yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
        az2,teff2,lat,cliq,sb,sb2,teff,iinc);
      
      test=tol*abs(yt)+tol*ymax/50;
      test2=ym+ymax/100;
      if (abs(yt-ym) <= test and abs(y[i-2]-y[i-1]) <= test2 and
        (x[i-2]-x[i-1]) < 0.5 ) { break; } // go to 120
      if (x[i-2]-x[i-1] < xtol) { break; } // go to 120
      i=i+1;
      x[i-1]=x[i-2];
      y[i-1]=y[i-2];
      x[i-2]=xm;
      y[i-2]=yt;
    }  // do 110 inner. This corresponds to    // go to 110 

    while ( do120 ){
      std::cout << "120" << std::endl;

      // 120 continue
      sum=sum+0.5*(y[i-1]+yl)*(x[i-1]-xl);
      xl=x[i-1];
      yl=y[i-1];
      i=i-1;

      // if (i > 1) go to 110
      if (i > 1) { break; }

      // if (i == 1) go to 120
      if (i != 1) { // don't go to 120 - either go to 130 or return 
        s[0]=sum; 

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
  rnbin=1.0/nbin;
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
  x[2]=-1;
  xl=x[2];
  y[2]=sig(e,ep,x[2],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  if (ep == 0.0) x[1]=0;
  if (ep != 0.0) x[1]=0.5*(e+ep-(s1bb-1)*az*tev)*seep;
  if (abs(x[1]) > 1-eps) x[1]=0.99e0;
  y[1]=sig(e,ep,x[1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  x[0]=+1;
  y[0]=sig(e,ep,x[0],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  ymax=y[0];
  if (y[1] > ymax) ymax=y[1];
  if (y[2] > ymax) ymax=y[2];
  if (ymax < eps) ymax=eps;

}




void do150( double& xm, double& ym, std::vector<double>& x, std::vector<double>& y,
  double& yt, double& e, double& ep, double& tev, int nalpha, int nbeta,
  std::vector<double>& alpha, std::vector<double>& beta, std::vector<std::vector<double>>& sab,
  double& test, int& i, double& test2, double& ymax, int iinc, double teff,
  double tol, double teff2, double az2, double xtol, int lasym, double tevz, 
  double az, int lat, int bbm, double cliq, double sb, double sb2, int imax ){

  while ( i < imax ){
    // 150 continue
    std::cout << "150" << std::endl;
    //if (i == imax) go to 160
    xm=0.5*(x[i-2]+x[i-1]);
    ym=0.5*(y[i-2]+y[i-1]);
    yt=sig(e,ep,xm,tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
      az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  
    test=tol*abs(yt)+tol*ymax/50;
    test2=ym+ymax/100;
    if (abs(yt-ym) <= test and abs(y[i-2]-y[i-1]) <= test2 and
      (x[i-2]-x[i-1]) < 0.5) return;
    if (x[i-2]-x[i-1] < xtol) return;
    i=i+1;
    x[i-1]=x[i-2];
    y[i-1]=y[i-2];
    x[i-2]=xm;
    y[i-2]=yt;
    // go to 150
  }
}


int do250( std::vector<double>& x, std::vector<double>& y, double& xl, double& yl,
  int& i ){
  std::cout << "250" << std::endl;
  // go to 250
  //  250 continue
  xl=x[i-1];
  yl=y[i-1];
  i=i-1;
  if (i > 1) {
    // go to 150
    std::cout << "go to 150 from 250" << std::endl;
    return 150;
    //continue;
  }
  if (i == 1) {
    // go to 160
    std::cout << "go to 160 from 250" << std::endl;
    return 160;
  }

  std::cout << "go to 260 from 250" << std::endl;
  // 260 continue
  return 260;

}


int do160(double add, std::vector<double>& x, std::vector<double>& y, double& xl,
  double& yl, int& i, double& xil, int& j, double& fract, int& nbin, double& sum,
  double& gral, double& xn, double& shade ){

  while ( true ){
    std::cout << "160" << std::endl;
    add=0.5*(y[i-1]+yl)*(x[i-1]-xl);

    if (x[i-1] == xl) {
      // returns either 150, 160, or 260
      int what_next = do250( x, y, xl, yl, i ); 
      if ( what_next == 160 ){ continue; }
      return what_next;
    } 

    xil=1/(x[i-1]-xl);
   
    if (i == 1 and j == nbin-1) {
      // go to 165
      std::cout << "165" << std::endl;
      // 165 continue
      xn=x[i-1];
      j=j+1;
      // go to 190
      return 190;

    }

    if (sum+add >= fract*shade and j < nbin-1){ 
      // go to 170
      return 170;

      //std::cout << "go to 190" << std::endl;
       
    }
    sum=sum+add;
    gral=gral+0.5*(yl*x[i-1]-y[i-1]*xl)*(x[i-1]+xl)
      +(y[i-1]-yl)*(x[i-1]*x[i-1]+x[i-1]*xl+xl*xl)/3.0;
     
    // go to 250
    // returns either 150, 160, or 260
    int what_next = do250( x, y, xl, yl, i ); 
    if ( what_next != 160 ){ 
      // returns either 150 or 260
      return what_next;
    }
  }
}



auto do170(int& j, double& test, double& fract, double& sum, std::vector<double>& y,
  std::vector<double>& x, double& yl, double& xn, double& xl, double& f, double& disc,
  double& ytol, double& rf, int& i, double& xil, const double& sigmin){
  // 170 continue
  std::cout << 170 << std::endl;
  j=j+1;

  if (abs(test) > ytol){ 
    test=(fract-sum)*(y[i-1]-yl)/((x[i-1]-xl)*yl*yl);
  }
  if (abs(test) > ytol or yl < sigmin ){
    std::cout << 175 << std::endl;
    // 175 continue
    f=(y[i-1]-yl)*xil;
    rf=1/f;
    disc=(yl*rf)*(yl*rf)+2*(fract-sum)*rf;
    if (disc < 0.0) {
       // write(strng,'(''disc='',1p,e12.4)') disc
       // call mess('sigl',strng,'set to abs value and continue')
       disc=abs(disc);
    } // end if
    if (f > 0.0) xn=xl-(yl*rf)+sqrt(disc);
    if (f < 0.0) xn=xl-(yl*rf)-sqrt(disc);
    if (xn > xl and xn <= x[i-1]) { 
      // go to 190
      return;
    }
    if (xn > xl and xn < (x[i-1]+ytol*(x[i-1]-xl))) {
      // go to 180
      std::cout << 180 << std::endl;
      // 180 continue
      xn=x[i-1];
      // go to 190
      return;
    } 
    std::cout << "call error('sigl','no legal solution (quadratic path).',' ')" << std::endl;

  }

  xn=xl+(fract-sum)/yl;
  if (xn > x[i-1]) {
    // go to 180
    std::cout << 180 << std::endl;
    // 180 continue
    xn=x[i-1];
    // go to 190
    return;
  }

  // if (xn >= xl and xn <= x(i)) go to 190
  if (xn < xl or xn > x[i-1]) {
    std::cout << "call error('sigl','no legal solution.',' ')" << std::endl;
  }
  // go to 190
  return;
}


auto sigl( int nlin, int nalpha, int nbeta, int nlmax, double& e, double& ep,
  double& tev, std::vector<double>& alpha, std::vector<double>& beta,
  std::vector<std::vector<double>>& sab, std::vector<double>& s, double& tolin,
  double& az, double& tevz, int iinc, int lat, double& bbm, 
  int lasym, double& az2, double& teff2, double& cliq, double& sb,
  double& sb2, double& teff ){

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
  double xtol = 0.00001;
  double ytol = 0.001;
  double sigmin = 1.0e-32, eps = 1.0e-3;
  double shade = 0.99999999; 

  // constant factors
  b=(ep-e)/tev;
  tol=0.5*tolin;
  nl=nlin;
  if (nl < 0) nl=-nl;
  b=abs(b);
  if (lat == 1 and iinc == 2) b=b*tev/tevz;
  s1bb=sqrt(1+b*b);
  if (ep != 0.0) seep=1/sqrt(e*ep);

  // adaptive calculation of cross section
  i=3;
  sum=0;
  x[2]=-1;
  xl=x[2];
  y[2]=sig(e,ep,x[2],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
      az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  yl=y[2];
  if (ep == 0.0) x[1]=0;
  if (ep != 0.0) x[1]=0.5*(e+ep-(s1bb-1)*az*tev)*seep;
  if (abs(x[1]) > 1-eps) x[1]=0.99;
  y[1]=sig(e,ep,x[1],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
      az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  x[0]=+1;
  y[0]=sig(e,ep,x[0],tev,nalpha,alpha,nbeta,beta,sab,bbm,az,tevz,lasym,
      az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  ymax=y[1];
  if (y[0] > ymax) ymax=y[0];
  if (y[2] > ymax) ymax=y[2];
  if (ymax < eps) ymax=eps;


  do_110_120_130(i, imax, xm, ym, x, y, e, ep, tev, tevz, nalpha, nbeta, alpha,beta, sab, bbm, az, az2, lasym, teff, teff2, lat, cliq, sb, sb2, iinc, yt, sum, nl, sigmin, s, nbin, fract, rnbin, xl, j, ymax, eps, rfract, seep, gral, yl, test2, s1bb, test, tol, xtol);


  bool go_straight_to_150_from_190 = true;
  while ( true ){ 
    if (go_straight_to_150_from_190){
      do150( xm, ym, x, y, yt, e, ep, tev, nalpha, nbeta, alpha, beta, sab, test, i, test2, ymax, iinc, teff, tol, teff2, az2, xtol, lasym, tevz, az, lat, bbm, cliq, sb, sb2, imax );
    }
    go_straight_to_150_from_190 = true; 

    // check bins for this panel

    int what_next = do160(add, x, y, xl, yl, i, xil, j, fract, nbin, sum, gral, xn, shade );
    if (what_next == 150){ continue; }
    if (what_next == 260){ return; }
    if (what_next == 170 ){ do170(j, test, fract, sum, y, x, yl, xn, xl, f, disc, ytol, rf, i, xil, sigmin); }

    //190 continue
    std::cout << "190" << std::endl;
    yn=yl+(y[i-1]-yl)*(xn-xl)*xil;
    gral=gral+(xn-xl)*(yl*0.5*(xn+xl)
      +(y[i-1]-yl)*xil*(-xl*0.5*(xn+xl)
      +(xn*xn+xn*xl+xl*xl)/3.0));
    xbar=gral*rfract;

    // compute legendre components
    if (nlin >= 0) {
       legndr(xbar,p,nl);
       for ( int il = 1; il < nl; ++il ){
          s[il]=s[il]+p[il]*rnbin;
       } // end do
    }
    // output equally probable angles
    else {
       s[j+1]=xbar;
    } // end if

    // continue bin loop and linearization loop
    xl=xn;
    yl=yn;
    sum=0;
    gral=0;
    if (j == nbin) return;
    if (xl < x[i-1]) {
      // go to 160
      go_straight_to_150_from_190 = false; 
      continue;
    
    }
   
    // returns either 150, 160, or 260
    what_next = do250( x, y, xl, yl, i ); 
    if ( what_next == 160 ){  
      go_straight_to_150_from_190 = false; 
      continue;
    }
    if ( what_next == 150 ){
      continue;
    }
    if (what_next == 260 ){
      return;
    }

  }

}
