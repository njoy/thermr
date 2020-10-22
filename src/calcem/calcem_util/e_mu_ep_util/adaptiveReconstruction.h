#include "general_util/sigfig.h"
#include "generalTools/constants.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"
#include "ENDFtk.hpp"

using namespace njoy::ENDFtk;
using LaboratoryAngleEnergy = section::Type<6>::LaboratoryAngleEnergy;
using AngularDistribution   = section::Type<6>::LaboratoryAngleEnergy::AngularDistribution;
using EnergyDistribution    = section::Type<6>::LaboratoryAngleEnergy::EnergyDistribution;

template <typename A, typename F>
auto addToStack(int& i, A& x, A& yy, const A& yu, const F& xm ){
  // test fails. add to stack and continue (575)
  ++i;
  x[i-1]  = x[i-2];
  x[i-2]  = xm;
  yy[i-1] = yy[i-2];
  yy[i-2] = yu[0];
}



template <typename Float, typename Range>
auto adaptiveReconstruction( const Range& alphas, const Range& betas, const Range& sab,
  int iinc, const Range& egrid, const Float& temp, int nbin, const double& emax,
  const Float& tol, int lat, int lasym, const Float& az, const Range& boundXsVec,
  const Float& teff ){

  Range x(20,0), y(20,0), yu(100,0), yy(20,0), uj(300,0), sj(300,0), 
        ubar(egrid.size(),0);

  Float xm, ym, tolmin = 5e-7, sum = 0, u = -1, tev = temp*kb; 
  
  int i = 2, j = 0, nne = egrid.size();

  for (size_t i = 0; i < egrid.size(); ++i){
    if (egrid[i] > emax){
      nne = i + 1;
      break; 
    }
  }

  Range esi(nne,0.0), xsi(nne,0.0);

  std::vector<AngularDistribution> angularDistVec;

  for (size_t ie = 0; ie < esi.size(); ++ ie){
    Float enow = egrid[ie];
    if (ie > 0 and temp > 3000.0){ enow = enow*temp/3000.0; }

    enow = sigfig(enow,8,0);
    esi[ie] = enow; 

    j = 0;
    sum = 0; 

    // adaptive reconstruction of angular cross section
    u = -1;
    
    x[1] = u;
    sigu(enow,u,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,yu);
    yy[1] = yu[0];

    Float xl =  x[1],
          yl = yy[1];

    u = 1.0;
    x[0] = u;
    sigu(enow,u,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,yu);
    yy[0] = yu[0];

    i = 2;

    while (i > 1){ // 530 

      xm = 0.5*(x[i-2]+x[i-1]);
      xm = sigfig(xm,7,0);

      if (i < 20 and xm > x[i-1] and xm < x[i-2]){ 
        sigu(enow,xm,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,yu);
  
        ym = yy[i-1]+(xm-x[i-1])*(yy[i-2]-yy[i-1])/(x[i-2]-x[i-1]);
        if (x[i-2]-x[i-1] > 0.25 or abs(yu[0]-ym) > 2.0*tol*ym + tolmin){
          addToStack(i, x, yy, yu, xm );
          continue;
        }
      } 
      // 560
      j++;
      uj[j-1] = x[i-1];
      sj[j-1] = yy[i-1];
      if (j > 1){
        sum += 0.5*(yy[i-1]+yl)*(x[i-1]-xl);
        xl = x[i-1];
        yl = yy[i-1];
      }
      --i;
    } // while i > 1

    // 580
    ++j;
    uj[j-1] = x[0];
    sj[j-1] = yy[0];
    int numMus= j;
    sum += 0.5*(yy[0]+yl)*(x[0]-xl);
    xsi[ie-1] = sum*0.5;

    for (int i = 2; i <= numMus; ++i){
      ubar[ie] += 0.5*(uj[i-1]-uj[i-2])*(sj[i-1]+sj[i-2])*(uj[i-1]+uj[i-2]);
    }
    ubar[ie] = 0.5*ubar[ie]/sum;

    // Add in a dummy energy distribution to initialize the vector
    std::vector<EnergyDistribution> energyDistributionsVec {
        EnergyDistribution( 1.0, { 4 }, { 2 }, { 1e-5, 1.1e+7, 1.147e+7, 3e+7 },
                                                { 0., 2., 4., 6. } ) };
    for (int il = 0; il < numMus; ++il ){
      u = uj[il];
      sigu(enow,u,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,yu);
      int nep = int(yu[1]);
      j = 0;
      for (int i = 1; i <= nep; ++i){
        j = nep-i;
        if (yu[2*(nep-i)+4-1]/sum > 2e-7){
          break;
        }
      }
      nep = j;

      // 595
      
      if (nep > 153){ nep = 153; }
      Range EpVec(nep), ProbVec(nep);

      // 596
      
      for (size_t ib = 1; ib < nep+1; ++ib){
        EpVec[ib-1]   = yu[  2*ib];
        ProbVec[ib-1] = yu[1+2*ib]*2/sum;
      }

      std::vector<long> interpolants {(long) EpVec.size()}, boundaries {2};

      if ( il == 0 ){
        energyDistributionsVec[0] = EnergyDistribution ( u, 
          std::move(interpolants), std::move(boundaries), std::move(EpVec), 
          std::move(ProbVec) );
      }
      else {
        energyDistributionsVec.push_back(EnergyDistribution ( u, 
          std::move(interpolants), std::move(boundaries), std::move(EpVec), 
          std::move(ProbVec) ));
      }

    }

    std::vector<long> interpolants {2};
    std::vector<long> boundaries   {(long) energyDistributionsVec.size()};
    angularDistVec.push_back(AngularDistribution( enow, std::move(boundaries),
                std::move(interpolants), std::move(energyDistributionsVec) ) );

  }

  std::vector< long > boundaries = { (long) esi.size() };
  std::vector< long > interpolants = { 1 };

  LaboratoryAngleEnergy laboratoryAngleEnergy( 
    std::move( boundaries ),
    std::move( interpolants ),
    std::move( angularDistVec) );

  return laboratoryAngleEnergy;

} 


