#include <iostream>
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_ep_mu_util/sigl.h"
#include "coh/coh_util/sigcoh_util/legndr.h"



template <typename Range, typename Float>
auto do_313( const int& lat, int& jbeta, const Float& enow, //Float& ep, 
  const Range& betas, const Range& x, const Float& tev ){//, int& iskip ) {

  Float ep = 0.0;
  while ( true ){ 
    std::cout << " --- 313 --- " << std::endl;
    if ( jbeta == 0 ){ jbeta = 1; }
    if ( jbeta <= 0 ){ 
      if ( lat == 1 )  { ep = enow - betas[-jbeta-1]*0.0253; }
      else             { ep = enow - betas[-jbeta-1]*tev ; }
      if ( ep == enow ){ ep = sigfig(enow,8,-1); }
      else             { ep = sigfig(ep,  8, 0); }
    }
    else { 
      if ( lat == 1 )  { ep = enow + betas[jbeta-1]*0.0253; }
      else             { ep = enow + betas[jbeta-1]*tev ; }
      if ( ep == enow ){ ep = sigfig(enow,8, 1); }
      else             { ep = sigfig(ep,  8, 0); }
    }
    //if ( ep > x[1] ){ break; }
    if ( ep > x[1] ){
        return ep; 
    }
    jbeta += 1;
  }
  return ep;
} 




template <typename Range, typename Float>
auto do_410(int& i, Range& x, Range& y, Range& s, Float xm, int nl, int imax) {
  std::cout << " --- 410 ---" << std::endl;
  i += 1;
  x[i-1] = x[i-2];
  x[i-2] = xm;
  for ( int il = 0; il < nl; ++il ){
    y[il*imax+i-2] = y[il*imax+i-2];
    y[il*imax+i-2] = s[il];
  }

}

template <typename Range, typename Float> 
auto do_380( int& i, const int& imax, const int& j, int& jnz, int nl, Range& scr, 
  const Range& x, const Range& y, Float& ulast, Float& u2last, Float& u3last,
  Float& xlast, Float& ylast) {
  std::cout << " --- 380 --- " << std::endl;
  int jscr = 7 + (j-1)*(nl+1);
  scr[jscr-1] = x[i-1];
  if (y[i-1] >= 1e-9 ){ scr[jscr] = sigfig(y[i-1],9,0); }
  else                { scr[jscr] = sigfig(y[i-1],8,0); }

  for ( int il = 2; il <= nl; ++il ){
    scr[il+jscr-1] = sigfig(y[(il-1)*imax+i-1],9,0);
    if (scr[il+jscr-1] > 1.0 ){ 
      if (scr[il+jscr-1] > 1.0+0.0005 ){
        std::cout << "call mess????" << std::endl; throw std::exception(); }
      scr[il+jscr-1] = 1.0; }
    if ( scr[il+jscr-1] < -1.0 ){
      if (scr[il+jscr-1] < -(1.0+0.0005) ){
        std::cout << "call mess????" << std::endl; throw std::exception(); }
      scr[il+jscr-1] = -1.0; }
  }

  xlast = x[i-1];
  ylast = y[0*imax +  i-1];
  if (ylast != 0.0){ jnz = j; }
  ulast = 0.0;
  u2last = 0.0;
  u3last = 0.0;
  //nll = 3;
  Range p (4,0.0);
  for (int il = 2; il <= nl; ++il){
    legndr(y[(il-1)*imax+(i-1)],p,3);
    ulast  += p[1];
    u2last += p[2];
    u3last += p[3];
  }
  ulast  *= y[i-1]/(nl-1);
  u2last *= y[i-1]/(nl-1);
  u3last *= y[i-1]/(nl-1);
  i -= 1;
}



template <typename Range, typename Float>
auto do_360( Range& xsi, Range& x, Range& y, Float& xlast, Float& ylast, int& i, 
  int& j, int nl, Float& ulast, Float& u2last, Float& u3last,
  Range& ubar, Range& p2, Range& p3, int imax, int jmax, int ie, int& jnz, int& jbeta, 
  int nbeta, Range& scr){
  std::cout << " --- 360 --- " << std::endl;
  j += 1;
  Float uu,u2,u3;
  int nll;
  if ( j >= jmax ){ std::cout << "j is too big" << std::endl; throw std::exception(); }
  if ( j > 1 ){ 
      xsi[ie-1] += (x[i-1]-xlast)*(y[i-1]+ylast)*0.5;
      uu = 0;
      u2 = 0;
      u3 = 0;
      nll = 3;
      Range p (4,0.0);
      for ( int il = 1; il < nl; ++il ){
        legndr( y[il*imax+i-1], p, nll );
        uu += p[1];
        u2 += p[2];
        u3 += p[3];
      }
      uu /= (nl-1);
      uu *= y[i-1];
      u2 /= (nl-1);
      u2 *= y[i-1];
      u3 /= (nl-1);
      u3 *= y[i-1];
      ubar[ie-1] += 0.5*(x[i-1]-xlast)*(uu+ulast);
      p2[ie-1] += 0.5*(x[i-1]-xlast)*(u2+u2last);
      p3[ie-1] += 0.5*(x[i-1]-xlast)*(u3+u3last);


  }
  if ( j == 3 and xsi[ie-1] >= 5e-7 ){
    j = 2;
  }



  do_380( i, imax, j, jnz, nl, scr, x, y, ulast, u2last, u3last, xlast, ylast);

  if ( i >= 2 ){ 
      return 330;
  }
  jbeta += 1;
  if (jbeta <= nbeta){
      //std::cout << " got to 311 " << std::endl;
      return 311;
  }
  for ( int il = 0; il < nl; ++il ){
    y[il*imax+i-1] = 0.0;
  }
  //std::cout << " go to 430 " << std::endl;
  return 430;





}



template <typename Range, typename Float>
//template <typename Float>
auto e_ep_mu( Float T, Float& teff, Float& teff2, int jmax, int nne, int nnl, int nl, Float tol, Float& sigma_b, Float& sigma_b2, Float az, int lasym, int lat, int iinc, const Range& alphas, const Range& betas, const Range& sab ){
  // should only be here if iform = 0
  std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
    1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 2.53e-4, 2.97e-4, 3.5e-4, 4.2e-4, 5.06e-4, 
    6.15e-4, 7.5e-4, 8.7e-4, 1.012e-3, 1.23e-3, 1.5e-3, 1.8e-3, 2.03e-3, 
    2.277e-3, 2.6e-3, 3e-3, 3.5e-3, 4.048e-3, 4.5e-3, 5e-3, 5.6e-3, 6.325e-3, 
    7.2e-3, 8.1e-3, 9.108e-3, 0.01, 0.01063, 0.0115, 0.012397, 0.0133, 0.01417, 
    0.015, 0.016192, 0.0182, 0.0199, 0.020493, 0.0215, 0.0228, 0.0253, 0.028, 
    0.030613, 0.0338, 0.0365, 0.0395, 0.042757, 0.0465, 0.05, 0.056925, 0.0625, 
    0.069, 0.075, 0.081972, 0.09, 0.096, 0.1035, 0.111573, 0.12, 0.128, 0.1355, 
    0.145728, 0.16, 0.172, 0.184437, 0.2, 0.2277, 0.2510392, 0.2705304, 
    0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 
    0.5032575, 0.56, 0.625, 0.7, 0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 
    1.7, 1.855, 2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 
    4.9, 5.35, 5.85, 6.4, 7.0, 7.65, 8.4, 9.15, 9.85, 10.0 };

  Range scr(500000,0.0);
  Float xm, ym;
  Float ylast = 0.0, xlast = 0.0, ulast = 0.0, u2last = 0.0, u3last = 0.0;
  Float tolmin = 5e-7;
  int i;
  
  Float kb = 8.6173303E-5;
  Float tev = kb*T;
  teff  *= kb;
  teff2 *= kb;
  // write ( in the 6222 section of thermr output )
  // CONTIO
  // za  awr   0  5  1  0    (this is funky bc the 5 shows up as 1 in the output)
  // TAB1IO
  // 1    1   -1  1  1  2  
  // 2    2  
  // 1e-5 1 emax 1
  // TAB2IO
  // T    0    3  1  1  nne
  // nne  2
  
  int ie;
  int jnz = 0;
  Float enow, ep;

  int imax = 20;
  std::vector<double> esi(nne+1), xsi(nne+1), 
  ubar(egrid.size()), p2(egrid.size()), p3(egrid.size()), x(imax), y(65*imax,0.0); // This here is nlmax = 65

  int j, jbeta, iskip;
  // loop over given incident energy grid
  std::cout << " --- 305 --- " << std::endl;
  ie = 0;

  while (true){
    std::cout << " --- 310 --- " << std::endl;
    ie = ie+1;
    enow = egrid[ie-1];
    if ( T > 3000.0 ){
      std::cout << "do temp approx" << std::endl;
      Float tone = 0.0253/kb;
      Float elo = egrid[0];
      enow = elo*exp(log(enow/elo)*log((T/tone)*egrid[egrid.size()-1]/elo)/log(egrid[egrid.size()-1]/elo));
    }
    enow = sigfig(enow,8,0);
    esi[ie-1] = enow;
    xsi[ie-1] = 0.0;
    ubar[ie-1] = 0.0;
    p2[ie-1] = 0.0;
    p3[ie-1] = 0.0;
    ep = 0.0;
    x[0] = ep;

    auto s  = sigl(ep,enow,tev,tol,nnl,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);

    for ( int il = 0; il < nl; ++il ){
      y[il*imax+0] = s[il];
      //y[i*imax+j] = y(i,j) where this goes up to y(nlmax,imax)
    }
    jbeta = -betas.size();
    if (lasym > 0) jbeta = 1;
    j = 0;
    iskip = 0;


    int output_380 = 0;
    while (true){ 
      std::cout << " --- 311 --- " << std::endl;
      x[1] = x[0];
      for ( int il = 0; il < nl; ++il ){
        y[il*imax+1] = y[il*imax+0];
      }


      ep = do_313( lat, jbeta, enow, betas, x, tev );//, iskip );



      std::cout << " --- 316 ---" << std::endl; 
      ep = sigfig(ep,8,0);
      x[0] = ep;
      s = sigl(ep,enow,tev,tol,nnl,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
      for ( int il = 0; il < nl; ++il ){ y[il*imax+0] = s[il]; }

      // adaptive subdivision of panel
      i = 2;
      Float uu  = 0;
      Float uum = 0;
      
      while (true){ 
        std::cout << " --- 330 --- " << "   " << i << "   " << y[18] << std::endl;
        Float quickTest = 0.5 * ( y[0*imax+(i-1)-1] + y[0*imax+(i)-1] ) * ( x[(i-1)-1] - x[(i)-1] ); 
        xm = 0.5*(x[(i-1)-1]+x[(i)-1]); 
        xm = sigfig(xm,8,0);

        bool go_to_330 = false;
        if ( not ( i == imax or iskip == 1 or quickTest < tolmin or xm <= x[i-1] or xm >= x[i-2] )){ 
          // Don't immediately go to 360
          for ( int k = 0; k < nl; ++k ){
            ym = ( x[i-2] == x[i-1] ) ? 
              y[k*imax+i-1] :
              y[k*imax+i-1] + (xm-x[i-1])*(y[k*imax+i-2]-y[k*imax+i-1])/(x[i-2]-x[i-1]);
            
            if ( k > 0 ){ uu  += s[k]; uum += ym;  }

            Float test2 = ( k > 0 ) ? tol : tol*abs(s[k]);
            if ( abs(s[k]-ym) > test2 ){
              do_410(i, x, y, s, xm, nl, imax);
              go_to_330 = true;
            }
          }

          std::cout << " --- 350 --- " << std::endl;
          if ( abs(uu-uum) > 2.0*tol*abs(uu)+0.00001 ){ 
              do_410(i, x, y, s, xm, nl, imax);
              go_to_330 = true;
          }

        }
        if ( go_to_330 == true ){ continue; }

        if ( i != imax and iskip == 1 ){ iskip = 0; }
        output_380 = do_360( xsi, x, y, xlast, ylast, i, j,  nl, ulast, u2last, u3last, 
                             ubar, p2, p3, imax, jmax, ie, jnz, jbeta, betas.size(), scr );

        if ( output_380 == 330 ){ continue; }
        break;

      } // 330 LOOP
      if ( output_380 == 311 ){ continue; } 
      break;

    } // 311 LOOP

    std::cout << " --- 430 --- " << std::endl;
    return;

  } // 310 LOOP


}
