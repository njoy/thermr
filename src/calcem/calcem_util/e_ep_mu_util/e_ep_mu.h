#include <iostream>
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_ep_mu_util/sigl.h"
#include "coh/coh_util/sigcoh_util/legndr.h"

template <typename Float>
Float highTempApprox( const Float& T, const Float& enow, const Float& egrid_first, const Float& egrid_last ){
  std::cout << "do temp approx" << std::endl;
  Float kb = 8.6173303E-5;
  Float tone = 0.0253/kb, elo = egrid_first;
  return elo*exp(log(enow/elo)*log((T/tone)*egrid_last/elo)/log(egrid_last/elo));
}

template <typename Range, typename Float>
bool needMidpoint(const Range& x, const Range& y, const Float& xm, //const int& imax, 
  const int& i, const int& nl, const Range& s, const Float& tol ){

  int imax = x.size();
  Float quickTest = 0.5*(y[0*imax+i-2] + y[0*imax+i-1])*(x[i-2] - x[i-1]), ym;

  if ( not ( i == imax or quickTest < 5e-7 or xm <= x[i-1] or xm >= x[i-2] )){ 
    Float uu  = 0, uum = 0;
    for ( int k = 0; k < nl; ++k ){
      ym = ( x[i-2] == x[i-1] ) ? 
        y[k*imax+i-1] :
        y[k*imax+i-1] + (xm-x[i-1])*(y[k*imax+i-2]-y[k*imax+i-1])/(x[i-2]-x[i-1]);

      if ( k > 0 ){ uu  += s[k]; uum += ym;  }

      Float test2 = ( k > 0 ) ? tol : tol*abs(s[k]);
      if ( abs(s[k]-ym) > test2 ){ return true; } // need midpoint
    }

    // 350
    if (abs(uu-uum) > 2*tol*abs(uu)+0.00001){ return true; } // need midpoint
  }
  return false; // do not need to go to 410
}





template <typename Range, typename Float>
auto findFirstEprime( const int& lat, int& jbeta, const Float& E, //Float& Ep, 
  const Range& betas, const Range& x, const Float& tev ){//, int& iskip ) {
  // Given an incoming neutron energy E and some jbeta value, we're going to
  // look at all beta values (positive and negative) to see which is the first
  // one that could give us the desired outgoing energy (x[1]).
  //
  // Basically, we're trying to figure out the Ep the index of the lowest beta
  // value that can be invoked, and the corresponding E' value.

  using std::abs;
  int sign; Float Ep;

  do {
    //std::cout << " --- 313 --- " << std::endl;
    if ( jbeta == 0 ){ jbeta = 1; }
    sign = abs(jbeta) / jbeta;
    Ep = ( lat == 1 ) ? E + sign*betas[abs(jbeta)-1]*0.0253 
                      : E + sign*betas[abs(jbeta)-1]*tev ;
    Ep = ( Ep == E )  ? sigfig(E, 8, sign) : sigfig(Ep, 8, 0);
    // This is the E' that we get with a beta value of betas[|jbeta|-1]
    ++jbeta;
  } while ( Ep <= x[1] );

  jbeta -= 1;
  return Ep;
} 




template <typename Range, typename Float>
auto insertPoint(int& i, Range& x, Range& y, const Range& s, const Float& xm, int nl ) {
  std::cout << " --- 410 ---" << std::endl;
  i += 1;
  x[i-1] = x[i-2];
  x[i-2] = xm;
  int imax = x.size();
  for ( int il = 0; il < nl; ++il ){
    y[il*imax+i-1] = y[il*imax+i-2];
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
    if ( abs(scr[il+jscr-1]) > 1.0 ){ 
      if (abs(scr[il+jscr-1]) > 1.0+0.0005 ){std::cout<<"call mess"<<std::endl; throw std::exception(); }
      scr[il+jscr-1] = scr[il+jscr-1]/abs(il+jscr-1); }
  }

  xlast = x[i-1];
  ylast = y[0*imax +  i-1];
  if (ylast != 0.0){ jnz = j; }
  ulast = 0.0;
  u2last = 0.0;
  u3last = 0.0;

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
  --i;
}



template <typename Range, typename Float>
auto do_360( Range& xsi, Range& x, Range& y, Float& xlast, Float& ylast, int& i, 
  int& j, int nl, Float& ulast, Float& u2last, Float& u3last, Range& ubar,
  Range& p2, Range& p3, int imax, int jmax, int ie, int& jnz, int& jbeta, 
  int nbeta, Range& scr){

  std::cout << " --- 360 --- " << std::endl;
  ++j;

  if ( j >= jmax ){ std::cout << "j is too big" << std::endl; throw std::exception(); }

  if ( j > 1 ){ 
      xsi[ie] += (x[i-1]-xlast)*(y[i-1]+ylast)*0.5;
      Float uu = 0, u2 = 0, u3 = 0;
      Range p (4,0.0);
      for ( int il = 1; il < nl; ++il ){
        legndr( y[il*imax+i-1], p, 3 );
        uu += p[1]; u2 += p[2]; u3 += p[3];
      }
      uu *= y[i-1]/(nl-1); ubar[ie] += 0.5*(x[i-1]-xlast)*(uu+ulast);
      u2 *= y[i-1]/(nl-1); p2[ie]   += 0.5*(x[i-1]-xlast)*(u2+u2last);
      u3 *= y[i-1]/(nl-1); p3[ie]   += 0.5*(x[i-1]-xlast)*(u3+u3last);
  }

  if ( j == 3 and xsi[ie] >= 5e-7 ){ j = 2; }

  do_380( i, imax, j, jnz, nl, scr, x, y, ulast, u2last, u3last, xlast, ylast);

  if ( i >= 2 )      { return 330; }
  ++jbeta;
  if  ( jbeta <= nbeta){ return 311; }
  for ( int il = 0; il < nl; ++il ){ y[il*imax+i-1] = 0.0; }
  return 430;
}







template <typename Range, typename Float>
auto do_330( const Float& enow, Range& x, Range& y, int& i, const Float& tev, const Float& tol, const int lat, const int iinc, const int lasym, const Range& alphas, const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, const Float& sigma_b2, const Float& teff, const int nnl, const int nl ){
  int imax = x.size();
  while (true){ 
    std::cout << " --- 330 --- " << std::endl;

    Float xm = 0.5*(x[i-2]+x[i-1]); xm = sigfig(xm,8,0);
    Range s = sigl(xm,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,abs(nnl)-1,true);

    //std::cout << (s|ranges::view::all) << std::endl;
    if ( needMidpoint(x, y, xm, i, nl, s, tol) == true ){ 
      std::cout << "need point" << std::endl;

      /*
      */
      std::cout << y[0*imax+0] << std::endl;
      std::cout << y[0*imax+1] << std::endl;
      std::cout << y[0*imax+2] << std::endl;
      std::cout << std::endl;
      std::cout << y[1*imax+0] << std::endl;
      std::cout << y[1*imax+1] << std::endl;
      std::cout << y[1*imax+2] << std::endl;
      //std::cout << x[0] << "   " << x[1] << "   " << x[2]  << std::endl;

      std::cout << std::endl;

      // make it so that sigl can output the sum (PDF) so taht you can put thatin the 0*imax index of y

      std::cout << s[0] << "  " << s[1] << "  " << s[2] << std::endl;

      insertPoint(i, x, y, s, xm, nl );

      //std::cout << x[0] << "   " << x[1] << "   " << x[2]  << std::endl;
      std::cout << y[0*imax+0] << std::endl;
      std::cout << y[0*imax+1] << std::endl;
      std::cout << y[0*imax+2] << std::endl;
      std::cout << std::endl;
      std::cout << y[1*imax+0] << std::endl;
      std::cout << y[1*imax+1] << std::endl;
      std::cout << y[1*imax+2] << std::endl;
      /*
      */


      return;



      continue; 
    }
    return;

    //output_380 = do_360( xsi, x, y, xlast, ylast, i, j,  nl, ulast, u2last,
    //  u3last, ubar, p2, p3, imax, jmax, ie, jnz, jbeta, betas.size(), scr );

    if ( i < 2 ){ break; }

  } // 330 LOOP
}






















template <typename Range, typename Float>
auto e_ep_mu( Float T, Float& teff, Float& teff2, int jmax, int nne, int nnl, 
  int nl, Float tol, Float& sigma_b, Float& sigma_b2, Float az, int lasym, 
  int lat, int iinc, const Range& alphas, const Range& betas, const Range& sab ){
  // nl = nbin + 1 where nbin is the number of (user-defined) equiprobable 
  // angles that they want. 
  // nnl = -nl to the best of my knowledge. The negative is what flags us that 
  // we need to do the equiprobable angles and not the legendre components

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
  
  int jnz = 0;
  Float enow, ep;

  int imax = 20;
  std::vector<double> esi(nne+1), xsi(nne+1), 
  ubar(egrid.size()), p2(egrid.size()), p3(egrid.size()), x(imax), y(65*imax,0.0); // This here is nlmax = 65

  int j, jbeta;

  // loop over given incident energy grid
  std::cout << " --- 305 --- " << std::endl;
  for ( size_t ie = 0; ie < egrid.size(); ++ie ){
    std::cout << " --- 310 --- " << std::endl;

    enow = (T <= 3000) ? egrid[ie] 
                       : highTempApprox(T,egrid[ie],egrid[0],egrid[egrid.size()-1]);
    enow = sigfig(enow,8,0);

    esi[ie] = enow; xsi[ie] = 0.0; ubar[ie] = 0.0;
    p2 [ie] = 0.0;  p3[ie] = 0.0;  ep       = 0.0; x[0] = ep;

    auto s  = sigl(ep,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                   sigma_b2,teff,abs(nnl)-1,true);

    // Vector of xs for E->E' for nnl equiprobable angles (if nnl < 0) else the
    // legendre components (?)

    for ( int il = 0; il < nl; ++il ){
      y[il*imax+0] = s[il]; //y[i*imax+j] = y(i,j). max(i) = nlmax, max(j) = imax
    }
    jbeta = (lasym > 0) ? 1 : -betas.size();
    j     = 0;


    int output_380 = 0;
    while (true){ 
      std::cout << " --- 311 --- " << std::endl;
      x[1] = x[0];
      for ( int il = 0; il < nl; ++il ){
        y[il*imax+1] = y[il*imax+0];
      }

      ep = findFirstEprime( lat, jbeta, enow, betas, x, tev ); ep = sigfig(ep,8,0);

      std::cout << " --- 316 ---" << std::endl; 
      x[0] = ep;
      s = sigl(ep,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,abs(nnl)-1,true);
      for ( int il = 0; il < nl; ++il ){ y[il*imax+0] = s[il]; }

      // adaptive subdivision of panel
      i = 2;
     

      while (true){ 
        std::cout << " --- 330 --- " << "   " << i << "   " << y[18] << std::endl;

        Float xm = 0.5*(x[i-2]+x[i-1]); xm = sigfig(xm,8,0);
        s = sigl(xm,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,abs(nnl)-1,true);

        if ( needMidpoint(x, y, xm, imax, i, nl, s, tol) == true ){ 
          insertPoint(i, x, y, s, xm, nl, imax);
          continue; 
        }

        output_380 = do_360( xsi, x, y, xlast, ylast, i, j,  nl, ulast, u2last,
          u3last, ubar, p2, p3, imax, jmax, ie, jnz, jbeta, betas.size(), scr );

        if ( i < 2 ){ break; }
        //if ( i >= 2 ){ continue; }   //if ( output_380 == 330 ){ continue; }
        //break;

      } // 330 LOOP

      if ( output_380 == 311 ){ continue; } 
      break;

    } // 311 LOOP

    std::cout << " --- 430 --- " << std::endl;
    return;

  } // 310 LOOP


}
