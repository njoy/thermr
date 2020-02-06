#include "calcem/calcem_util/e_ep_mu_util/ep_mu.h"



template <typename Range, typename Float>
auto e_ep_mu( Range eVec, const Float& tev, const Float& tol, 
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

    auto out = prepareEpMu( eNow, j, tev, tol, lat, iinc, lasym, 
    alphas, betas, sab, az, sigma_b, sigma_b2, teff, nbin, jbeta, scr, lastVals, y);

    total_SCR[iEnergy]        = scr;
    total_OutputData[iEnergy] = out;

  }
  return std::make_tuple(eVec,total_SCR,total_OutputData);

}




