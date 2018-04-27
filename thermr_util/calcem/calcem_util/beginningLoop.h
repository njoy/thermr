#include "sig.h"

void do_110_120_130(int& i, int& imax, std::vector<double>& x,
  std::vector<double>& y, double& e, double& ep, double& tev, double& tevz, 
  std::vector<double>& alpha,std::vector<double>& beta,
  std::vector<std::vector<double>>& sab, double& bbm, double& az, double& az2,
  int& lasym, double& teff, double& teff2, int& lat, double& cliq, double& sb,
  double& sb2, int& iinc, double& sum, int& nl, double& sigmin,
  std::vector<double>& s, int& nbin, double& fract, double& xl,
  int& j, double& ymax, double& eps, double& seep, double& gral, 
  double& yl, double& s1bb, double& tol, double& xtol){


  double xm, ym, yt, test, test2;

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

      yt=sig(e,ep,xm,tev,alpha,beta,sab,bbm,az,tevz,lasym,
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
  fract=sum/nbin;
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
  y[2]=sig(e,ep,x[2],tev,alpha,beta,sab,bbm,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);
  if (ep == 0.0) x[1]=0;
  if (ep != 0.0) x[1]=0.5*(e+ep-(s1bb-1)*az*tev)*seep;
  if (abs(x[1]) > 1-eps) x[1]=0.99e0;
  y[1]=sig(e,ep,x[1],tev,alpha,beta,sab,bbm,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  x[0]=+1;
  y[0]=sig(e,ep,x[0],tev,alpha,beta,sab,bbm,az,tevz,lasym,
    az2,teff2,lat,cliq,sb,sb2,teff,iinc);

  ymax=y[0];
  if (y[1] > ymax) ymax=y[1];
  if (y[2] > ymax) ymax=y[2];
  if (ymax < eps) ymax=eps;

}

