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
  const int& i, const int& nl, const Range& s, const Float& tol, const Float& pdf ){

  int imax = x.size();
  Float quickTest = 0.5*(y[0*imax+i-2] + y[0*imax+i-1])*(x[i-2] - x[i-1]), ym;

  if ( not ( i == imax or quickTest < 5e-7 or xm <= x[i-1] or xm >= x[i-2] )){ 
    Float uu  = 0, uum = 0;

    for ( int k = 0; k < nl; ++k ){
      ym = ( x[i-2] == x[i-1] ) ? 
        y[k*imax+i-1] :
        y[k*imax+i-1] + (xm-x[i-1])*(y[k*imax+i-2]-y[k*imax+i-1])/(x[i-2]-x[i-1]);

      if ( k > 0 ){ uu  += s[k-1]; uum += ym;  }

      Float test2 = ( k > 0 ) ? tol : tol*abs(pdf);

      if ( k == 0 ){
       // std::cout << " -------- " << pdf << "    " << ym << "   " << abs(pdf-ym) << std::endl;
        if ( abs(pdf-ym) > test2 ){ return true; } // need midpoint
      }
      else {
        //std::cout << " -------- " << s[k-1] << "    " << ym << std::endl;
        if ( abs(s[k-1]-ym) > test2 ){ return true; } // need midpoint
      }
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
  do { //std::cout << " --- 313 --- " << std::endl;
    if ( jbeta == 0 ){ jbeta = 1; }
    sign = abs(jbeta) / jbeta;
    //std::cout << E << "   " << betas[abs(jbeta)-1]*tev << "    " << E+ betas[abs(jbeta)-1]*tev << std::endl;
    Ep = ( lat == 1 ) ? E + sign*betas[abs(jbeta)-1]*0.0253 
                      : E + sign*betas[abs(jbeta)-1]*tev ;
    //std::cout << " Ep (inside)" << Ep << std::endl;
    //Ep = ( Ep == E )  ? sigfig(E, 8, sign) : sigfig(Ep, 8, 0);
    Ep = sigfig(Ep,8,0);
    //std::cout << " Ep (inside)" << Ep << std::endl;
    // This is the E' that we get with a beta value of betas[|jbeta|-1]
    ++jbeta;
  } while ( Ep <= x[1] );

  jbeta -= 1;
  return Ep;
} 




template <typename Range, typename Float>
auto insertPoint(int& i, Range& x, Range& y, const Range& s, const Float& xm, int nl, const Float& pdf ) {
  int imax = x.size();
  //std::cout << " --- 410 --- " << std::endl;

  i += 1;
  x[i-1] = x[i-2];
  x[i-2] = xm;

  y[0*imax+i-1] = y[0*imax+i-2];
  y[0*imax+i-2] = pdf;

  for ( int il = 1; il < nl; ++il ){
    y[il*imax+i-1] = y[il*imax+i-2];
    y[il*imax+i-2] = s[il-1];
  }
}



template <typename Range, typename Float> 
auto getMoments(Float& ulast, Float& u2last, Float& u3last, const Range& y, const int& i, const int& imax, const int& nl ){
  ulast  = 0; 
  u2last = 0; 
  u3last = 0;

  int nll = 3;
  Range p(4,0.0);
  for ( int il = 1; il < nl; ++il ){
    legndr(y[il*imax+i-1],p,nll);
    ulast  += p[1];
    u2last += p[2];
    u3last += p[3];
  }
  ulast  *= (y[0*imax+i-1]/(nl-1));
  u2last *= (y[0*imax+i-1]/(nl-1));
  u3last *= (y[0*imax+i-1]/(nl-1));
  
}






template <typename Range, typename Float>
auto do_330( const Float& enow, Range& x, Range& y, int& j, const Float& tev, const Float& tol, const int lat, const int iinc, const int lasym, const Range& alphas, const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, const Float& sigma_b2, const Float& teff, const int nnl, const int nl, int& jbeta, Range& scr, Range& xsi, int& ie, Float& xlast, Float& ylast ){//, int& counter ){
  int imax = x.size();
  int i = 2;
  while (true){  // 330

    if ( i < imax ){
      Float xm = 0.5*(x[i-2]+x[i-1]); 
      xm = sigfig(xm,8,0);
      Range s(abs(nnl)-1,0.0);
      Float pdf = sigl(xm,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,s,true);
  
      if ( needMidpoint(x, y, xm, i, nl, s, tol, pdf) == true ){ 
        insertPoint(i, x, y, s, xm, nl, pdf );
        continue; 
      }
    }

    //std::cout << " --- 360 --- " << std::endl;
    ++j;
    if (j > 1) { xsi[ie] += (x[i-1]-xlast)*(y[0*imax+i-1]+ylast)*0.5; }
    //if (j > 1) { std::cout << " ---  " << (x[i-1]-xlast)*(y[0*imax+i-1]+ylast)*0.5 << std::endl; }
    //if (j > 1) { std::cout << " ---  " << xsi[ie] << std::endl; }
    if ( j == 3 and xsi[ie] < 5e-7 ){ j = 2; }


    int jscr = (j-1)*(nl+1)+1;

    if ( (unsigned) jscr+nl > scr.size() ){ scr.resize((jscr+nl)*2); }

    scr[jscr-1]=x[i-1];
    scr[jscr] = (y[0*imax+i-1] < 1e-9) ? 
                 sigfig(y[0*imax+i-1],8,0) 
               : sigfig(y[0*imax+i-1],8,0) ;

    xlast = x[i-1];
    ylast = y[0*imax+i-1];


    for ( int il = 1; il < nl; ++il ){
      scr[il+jscr] = sigfig(y[il*imax+i-1],9,0);
      if (scr[il+jscr] > 1.0){ y[il*imax+i-1] = 1.0; }
      if (scr[il+jscr] <-1.0){ y[il*imax+i-1] =-1.0; }
    }

    Float ulast, u2last, u3last;
    getMoments(ulast, u2last, u3last, y, i, x.size(), nl);


    --i;
    if ( i >= 2 ){ continue; } // go to 330

    ++jbeta;
    //if (jbeta > int(betas.size())){ // we will be going to 430 soon
    //  for ( auto& yVal : y ){ yVal = 0.0; }
    //  xsi[ie] += (x[i-1]-xlast)*(y[0*imax+i-1]+ylast)*0.5;
    //}

    return std::make_tuple(ulast,u2last,u3last);  // go to 311 or 430
    
  } // 330 LOOP
}









template <typename Range, typename Float>
auto do_330_extra( const Float& enow, int& j, const Float& tev, const Float& tol, const int lat, const int iinc, const int lasym, const Range& alphas, const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, const Float& sigma_b2, const Float& teff, const int nbin, int& jbeta, Range& scr, Range& xsi, int ie, Float& xlast, Float& ylast ){
    std::cout.precision(15);

  int nl = nbin + 1;
  int nnl = -nl;


  Range x(20,0.0), y(20*65,0.0);
  int imax = x.size();
  Float ep;

  do {

    x[1] = x[0];
  
    for ( int il = 0; il < nl; ++il ){ y[il*imax+1] = y[il*imax+0]; }

    ep = findFirstEprime( lat, jbeta, enow, betas, x, tev ); // 313
    ep = sigfig(ep,8,0);
    x[0] = ep;

    Range s(abs(nnl)-1,0.0);
    Float pdf = sigl(ep,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,s,true);
  
    y[0*imax+0] = pdf;
    for ( int il = 1; il < nl; ++il ){ y[il*imax+0] = s[il-1]; }

    do_330(enow,x,y,j,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff,nnl,nl,jbeta,scr,xsi,ie,xlast,ylast);


  } while( jbeta <= int(betas.size()));


  for ( auto& yVal : y ){ yVal = 0.0; }

  int jscr = (j)*(nl+1)+1;
  scr[jscr-1] = ep;
  scr.resize((j+1)*(nl+1)+1);

  double sum = 0.0;
  int lengthRow = nl+1;

  for ( size_t k = 0; k < scr.size(); ++k ){
    if ( (k-1)%lengthRow == 0 and (k-1) > 0 ){
      sum += (scr[k]+scr[k-lengthRow])*0.5*(scr[k-1]-scr[k-(lengthRow+1)]);
    }
  }
  for ( size_t k = 1; k < scr.size(); ++k ){
    if ( (k-1)%lengthRow == 0 ){ scr[k] = scr[k]/sum; }
  }


  // 430
  ++j;

  xsi[ie] += (x[0]-xlast)*(y[0*imax+0-1]+ylast)*0.5;




  return xsi[ie];

}









