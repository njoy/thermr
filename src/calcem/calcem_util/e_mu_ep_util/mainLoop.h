#include <iostream>
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"

template <typename Range, typename Float>
bool addPoint( Range& x, Range& yy, const Range& yu, const Float& xm, const Float& ym, const Float& tol, int& i ){
  if (x[i-2]-x[i-1] > 0.25 or (abs(yu[0]-ym ) > 2.0*tol*ym + 5e-7)) { 
    ++i;
    x[i-1] = x[i-2];
    x[i-2] = xm;
    yy[i-1] = yy[i-2];
    yy[i-2] = yu[0];
    return true;
  }
  return false;
}



template <typename Range, typename Float>
auto do_530_etc( Float enow, const Float& tev, const Float& tol, 
  const int lat, const int iinc, const int lasym, const Range& alphas, 
  const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, 
  const Float& sigma_b2, const Float& teff ){
  std::cout.precision(15);

  int imax = 20, mumax = 300, nemax = 50;
  Range x(imax,0.0); x[0] = 1.0; x[1] = -1.0;
  Range yy(imax,0.0), yu(2*nemax);
  Range uj(mumax,0.0), sj(mumax,0.);
  Float xm, ym, yl=0.0, xl=0.0;
  int i = 2;
  int j = 0;
  Float u = -1.0, sum = 0.0;

  u = -1.0;
  sigu( enow, u, tev, alphas, betas, sab, tol, az, iinc, lat, lasym, sigma_b, sigma_b2, teff, yu);
  yy[1] = yu[0];
  xl = x[1];
  yl = yy[1];

  u = 1.0;
  sigu( enow, u, tev, alphas, betas, sab, tol, az, iinc, lat, lasym, sigma_b, sigma_b2, teff, yu );
  yy[0]=yu[0];
  i = 2;

  do { // 530 
    if ( i != imax ){
      xm = 0.5*(x[i-2]+x[i-1]);
      xm = sigfig(xm,7,0);
      if ( xm > x[i-1] and xm < x[i-2] ){
        sigu( enow, xm, tev, alphas, betas, sab, tol, az, iinc, lat, lasym, 
              sigma_b, sigma_b2, teff, yu );
        ym = yy[i-1]+(xm-x[i-1])*(yy[i-2]-yy[i-1])/(x[i-2]-x[i-1]);

        if (addPoint(x,yy,yu,xm,ym,tol,i)){ continue; }

      }
    }
    
    ++j;
    uj[j-1] =  x[i-1];
    sj[j-1] = yy[i-1];
    if ( j > 1 ){ 
        sum += 0.5*(yy[i-1]+yl)*(x[i-1]-xl);
        xl =  x[i-1];
        yl = yy[i-1];
    }
    --i;

  } while ( i >= 2 );

  ++j; // 580
  uj[j-1] = x[0];
  sj[j-1] = yy[0];

  int nmu = j;

  sum += 0.5*(yy[0]+yl)*(x[0]-xl);

  Float xsi = sum*0.5;
  Float ubar = 0.0;
  for ( int i = 2; i <= nmu; ++i ){
    ubar += 0.5*(uj[i-1]-uj[i-2])*(sj[i-1]+sj[i-2])*(uj[i-1]+uj[i-2]);
  }
  ubar *= 0.5/sum;
  uj.resize(j);
  return std::make_tuple(xsi,ubar,uj);
}







template <typename Range, typename Float>
auto mu_ep( Float& enow, const Float& tev, const Float& tol, 
  const int lat, const int iinc, const int lasym, const Range& alphas, 
  const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, 
  const Float& sigma_b2, const Float& teff ){
  std::cout.precision(15);

  auto out = do_530_etc( enow, tev, tol, lat, iinc, lasym, alphas, betas, sab, az, sigma_b, 
              sigma_b2, teff );
  auto uj = std::get<2>(out);
  int imax = 20, mumax = 300, nemax = 5000;
  int i = 2, j = 0;
  Range x(imax,0.0); x[0] = 1.0; x[1] = -1.0;
  Range yy(imax,0.0), yu(2*nemax);
  Float xm, ym, yl=0.0, xl=0.0;
  Float u = -1.0, sum = 2.0*std::get<0>(out);
  std::vector<Range> totalSCR(uj.size());
  
  for (size_t il = 0; il < uj.size(); ++il ){

    Range scr(500,0.0);

    sigu( enow, uj[il], tev, alphas, betas, sab, tol, az, iinc, lat, lasym, 
               sigma_b, sigma_b2, teff, yu );

    int nep = int(yu[1]);
    j = 0;
    for (int i = 1; i <= nep; ++i ){
      j = nep-i;
      if (yu[2*(nep-i)+4-1]/sum > 2e-7){ break; }
    }
    nep = j;
    int istart = 1;
    int nw;
    int iend = nep;
    int ib = istart - 1;
    j = 0;
    do { // std::cout << " --- 596 --- " << std::endl;
      j  += 2;
      ib += 1;
      if ( int(scr.size()) < j ){ scr.resize(2*scr.size()); }
      scr[j-1-1] = yu[1+2*ib-1];
      scr[j+0-1] = yu[2+2*ib-1]*2/sum;
    } while ( ib < iend );
    scr.resize(nep*2);
    totalSCR[il] = scr;
  }
  return std::make_tuple(uj,totalSCR);
}












  



