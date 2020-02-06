#include <iostream>
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_ep_mu_util/sigl.h"
#include "coh/coh_util/sigcoh_util/legndr.h"

template <typename Float>
Float highTempApprox( const Float& T, const Float& enow, const Float& egrid_first, const Float& egrid_last ){
  Float kb = 8.6173303E-5;
  Float tone = 0.0253/kb;
  return egrid_first*exp(log(enow/egrid_first) 
       * log((T/tone)*egrid_last/egrid_first)/log(egrid_last/egrid_first));
}

template <typename Range, typename Float>
bool needMidpoint(const Range& x, const Range& y, const Float& xm, const int i, 
  const int nbin, const Range& s, const Float& tol, const Float& pdf ){

  int imax = x.size();
  Float quickTest = 0.5*(y[0*imax+i-2] + y[0*imax+i-1])*(x[i-2] - x[i-1]), ym;

  if ( not ( i == imax or quickTest < 5e-7 or xm <= x[i-1] or xm >= x[i-2] )){ 
    Float uu  = 0, uum = 0;

    for ( int k = 0; k < nbin+1; ++k ){
      ym = ( x[i-2] == x[i-1] ) ? 
        y[k*imax+i-1] :
        y[k*imax+i-1] + (xm-x[i-1])*(y[k*imax+i-2]-y[k*imax+i-1])/(x[i-2]-x[i-1]);

      if ( k > 0 ){ uu  += s[k-1]; uum += ym;  }

      Float test2       = ( k  > 0 ) ? tol : tol*abs(pdf),
            valueToTest = ( k == 0 ) ? pdf : s[k-1];
      if ( abs(valueToTest-ym) > test2 ){ return true; }

    }

    if (abs(uu-uum) > 2*tol*abs(uu)+0.00001){ return true; } // 350
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
  int sign; Float Ep; // 313
  do { 
    if ( jbeta == 0 ){ jbeta = 1; }
    sign = abs(jbeta) / jbeta;
    Ep = ( lat == 1 ) ? E + sign*betas[abs(jbeta)-1]*0.0253 
                      : E + sign*betas[abs(jbeta)-1]*tev ;
    Ep = sigfig(Ep,8,0);
    ++jbeta;
  } while ( Ep <= x[1] );

  jbeta -= 1;
  return Ep;
} 




template <typename Range, typename Float>
auto insertPoint(int& i, Range& x, Range& y, const Range& s, const Float& xm, int nbin, const Float& pdf ) {
  int imax = x.size(); // 410 
  ++i;
  x[i-1] = x[i-2];
  x[i-2] = xm;

  y[0*imax+i-1] = y[0*imax+i-2];
  y[0*imax+i-2] = pdf;

  for ( int il = 1; il < nbin+1; ++il ){
    y[il*imax+i-1] = y[il*imax+i-2];
    y[il*imax+i-2] = s[il-1];
  }
}



template <typename Range>
auto addToMoments( const Range& y, Range& vecOfVals, const int nbin, const int imax, const int i ){
  Range p(4,0.0);
  for (int il = 1; il < nbin+1; ++il){
    legndr(y[il*imax+i-1],p,3);
    vecOfVals[2] += p[1]; // uu
    vecOfVals[3] += p[2]; // u2
    vecOfVals[4] += p[3]; // u3
  }
  vecOfVals[2] *= y[0*imax+i-1]/(nbin);
  vecOfVals[3] *= y[0*imax+i-1]/(nbin);
  vecOfVals[4] *= y[0*imax+i-1]/(nbin);
}







template <typename Range, typename Float>
auto do_330( const Float& enow, Range& x, Range& y, int& j, const Float& tev, 
  const Float& tol, const int lat, const int iinc, const int lasym, 
  const Range& alphas, const Range& betas, const Range& sab, const Float& az, 
  const Float& sigma_b, const Float& sigma_b2, const Float& teff, 
  const int nbin, int& jbeta, Range& scr, Range& out, Range& lastVals ){
  int imax = x.size();

  int i = 2; // 330 
  while ( i >= 2 ){

    if ( i < imax ){
      Float xm = 0.5*(x[i-2]+x[i-1]); 
      xm = sigfig(xm,8,0);
      Range s(nbin,0.0);
      Float pdf = sigl(xm,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,
                       sigma_b,sigma_b2,teff,s,true);
  
      if ( needMidpoint(x, y, xm, i, nbin, s, tol, pdf) == true ){ 
        insertPoint(i, x, y, s, xm, nbin, pdf );
        continue; 
      }
    }

    ++j; // 360 
    //std::cout << "-   " << y[1*imax+0] << "  " << y[1*imax+1] <<  "   " << y[1*imax+2] << std::endl;
    //std::cout << std::endl;


    if (j > 1) { 
      Range currentMoments (lastVals.size(),0.0);
      addToMoments( y, currentMoments, nbin, imax, i );
      out[0] += 0.5*(x[i-1]-lastVals[0])*(y[0*imax+i-1]    +lastVals[1]); // xs
      out[1] += 0.5*(x[i-1]-lastVals[0])*(currentMoments[2]+lastVals[2]); // ubar
      out[2] += 0.5*(x[i-1]-lastVals[0])*(currentMoments[3]+lastVals[3]); // p2
      out[3] += 0.5*(x[i-1]-lastVals[0])*(currentMoments[4]+lastVals[4]); // p3
    }
 
    if ( j == 3 and out[0] < 5e-7 ){ j = 2; }
    int jscr = (j-1)*(nbin+2)+1;
    if ( (unsigned) jscr+nbin+1 > scr.size() ){ scr.resize((jscr+nbin+1)*2); }

    scr[jscr-1]=x[i-1];
    scr[jscr] = (y[0*imax+i-1] < 1e-9) ? 
                 sigfig(y[0*imax+i-1],8,0) 
               : sigfig(y[0*imax+i-1],8,0) ;

    for ( int il = 1; il < nbin+1; ++il ){
      scr[il+jscr] = sigfig(y[il*imax+i-1],9,0);
      if (scr[il+jscr] > 1.0){ y[il*imax+i-1] = 1.0; }
      if (scr[il+jscr] <-1.0){ y[il*imax+i-1] =-1.0; }
    }

   
    lastVals[0] = x[i-1]; lastVals[1] = y[0*imax+i-1];
    lastVals[2] = 0.0;    lastVals[3] = 0.0;    lastVals[4] = 0.0;

    addToMoments( y, lastVals, nbin, imax, i );


    --i;
  }

  ++jbeta;

}









template <typename Range, typename Float>
auto do_330_extra( const Float& enow, int& j, const Float& tev, const Float& tol, 
  const int lat, const int iinc, const int lasym, const Range& alphas, 
  const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, 
  const Float& sigma_b2, const Float& teff, const int nbin, int& jbeta, 
  Range& scr, Range& lastVals, Range& y ){

  Range x(20,0.0), out(4,0.0);
  int imax = x.size();
  Float ep;

  do {
    x[1] = x[0];
    for ( int il = 0; il < nbin+1; ++il ){ y[il*imax+1] = y[il*imax+0]; }

    ep = findFirstEprime( lat, jbeta, enow, betas, x, tev ); // 313
    ep = sigfig(ep,8,0);
    x[0] = ep;
    Range s(nbin,0.0);
    Float pdf = sigl(ep,enow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                     sigma_b2,teff,s,true);
    y[0*imax+0] = pdf;
    for ( int il = 1; il < nbin+1; ++il ){ y[il*imax+0] = s[il-1]; }
    do_330(enow,x,y,j,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,
           teff,nbin,jbeta,scr,out,lastVals);

  } while( jbeta <= int(betas.size()));

  for ( int il = 0; il < nbin+1; ++il ){
    y[il*imax] = 0.0;
  }

  scr[(j)*(nbin+2)] = ep;
  scr.resize((j+1)*(nbin+2));

  Float sum = 0.0;
  int lengthRow = nbin+2;

  for ( size_t k = 0; k < scr.size(); ++k ){
    if ( (k-1)%lengthRow == 0 and (k-1) > 0 ){
      sum += (scr[k]+scr[k-lengthRow])*0.5*(scr[k-1]-scr[k-(lengthRow+1)]);
    }
  }
  for ( size_t k = 1; k < scr.size(); ++k ){
    if ( (k-1)%lengthRow == 0 ){ scr[k] = scr[k]/sum; }
  }

  ++j; // 430

  out[0] += (x[0]-lastVals[0])*(y[0*imax+1-1]+lastVals[1])*0.5;
  out[1] += 0.5*(x[0]-lastVals[0])*(0.0+lastVals[2]); 
  out[2] += 0.5*(x[0]-lastVals[0])*(0.0+lastVals[3]); 
  out[3] += 0.5*(x[0]-lastVals[0])*(0.0+lastVals[4]); 
  Float invXS = 1.0/out[0];
  out[0] = sigfig(out[0],9,0);
  out[1] = sigfig(out[1]*invXS,5,0);
  out[2] = sigfig(out[2]*invXS,5,0);
  out[3] = sigfig(out[3]*invXS,5,0);

  return out;
}



template <typename Range, typename Float>
auto e_ep_mu_MAIN( Range eVec, const Float& tev, const Float& tol, 
  const int lat, const int iinc, const int lasym, const Range& alphas, 
  const Range& betas, const Range& sab, const Float& az, const Float& sigma_b, 
  const Float& sigma_b2, const Float& teff, const int nbin, const Float& temp ){
  std::cout.precision(15);

  Range lastVals(5,0.0);
  Float eNow, ePrime;
  int imax = 20;
  Range y(20*65,0.0);


  std::vector<Range> total_SCR(eVec.size());
  std::vector<Range> total_OutputData(eVec.size());


  for ( size_t iEnergy = 0; iEnergy < eVec.size(); ++iEnergy ){
    Range scr(y.size()*5,0.0);
    eNow = eVec[iEnergy];
    if (temp > 3000.0){ eNow = highTempApprox(temp,eNow,eVec[0],eVec[eVec.size()-1]); }
    eNow = sigfig(eNow,8,0);
    eVec[iEnergy] = eNow;

    ePrime = 0.0;
    Range s(nbin,0.0);
    Float pdf = sigl(0.0,eNow,tev,tol,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                     sigma_b2,teff,s,true);
    y[0*imax+0] = pdf;
    for (int il = 1; il < nbin+1; ++il){y[il*imax+0] = s[il-1];}

    int j = 0, jbeta = -int(betas.size());

    auto out = do_330_extra( eNow, j, tev, tol, lat, iinc, lasym, 
    alphas, betas, sab, az, sigma_b, sigma_b2, teff, nbin, jbeta, scr, lastVals, y);

    total_SCR[iEnergy]        = scr;
    total_OutputData[iEnergy] = out;

  }
  return std::make_tuple(eVec,total_SCR,total_OutputData);



}








