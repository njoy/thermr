#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"
#include "ENDFtk.hpp"



using namespace njoy::ENDFtk;
//using Variant               = njoy::ENDFtk::section::Type<6>::ContinuumEnergyAngle::Variant;
//using ContinuumEnergyAngle  = njoy::ENDFtk::section::Type<6>::ContinuumEnergyAngle;
//using ThermalScatteringData = njoy::ENDFtk::section::Type<6>::ContinuumEnergyAngle::ThermalScatteringData;
using LaboratoryAngleEnergy =
section::Type< 6 >::LaboratoryAngleEnergy;
using AngularDistribution =
section::Type< 6 >::LaboratoryAngleEnergy::AngularDistribution;
using EnergyDistribution =
section::Type< 6 >::LaboratoryAngleEnergy::EnergyDistribution;



template <typename A, typename F>
auto do575(int& i, A& x, A& yy, const A& yu, const F& xm ){
  // test fails. add to stack and continue
  // 575 continue
  //std::cout << 575 << std::endl;
  ++i;
  //std::cout << "i      " << i << std::endl;
  x[i-1] = x[i-2];
  x[i-2] = xm;
  yy[i-1]= yy[i-2];
  yy[i-2]= yu[0];

}



template <typename Float, typename Range>
auto adaptiveReconstruction( const Range& alphas, const Range& betas, const Range& sab,
  int iinc, const Range& egrid, const Float& temp, int nbin, const double& emax,
  const Float& tol, int lat, int lasym, const Float& az, const Range& boundXsVec,
  const Float& teff ){

  std::cout.precision(15);
  Range x(20,0.0), y(20,0.0), yu(100,0.0), yy(20,0.0),
        uj(300,0.0), sj(300,0.0), ubar(egrid.size(),0.0);// scr(500000);
  Float xm, ym, tolmin = 5e-7;
  
  int i = 2, j = 0;
  Float sum = 0, u = -1;

  int nne = 0;
  for (size_t i = 0; i < egrid.size(); ++i){
    if (egrid[i] > emax){
      nne = i + 1;
      break; 
    }
  }

  Range esi(nne,0.0);
  Range xsi(nne,0.0);
  Float bk = 8.6173303E-5;
  Float tev = temp * bk; 


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

    Float xl =  x[1];
    Float yl = yy[1];

    u = 1.0;
    x[0] = u;
    sigu(enow,u,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,yu);
    yy[0] = yu[0];


    i = 2;


    bool do_530 = true;
    while (do_530){
      // 530 
      //if ( i == 20 ){ std::cout << "go to 560" << std::endl; }
      xm = 0.5*(x[i-2]+x[i-1]);
      xm = sigfig(xm,7,0);
      //if (xm <= x[i-1] or xm >= x[i-2]){ std::cout << " go to 560 " << std::endl; }
      if (i < 20 and xm > x[i-1] and xm < x[i-2]){ 
        sigu(enow,xm,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,yu);
  
        if (x[i-2]-x[i-1] > 0.25){ 
          //std::cout << "go to 575" << std::endl; 
          do575(i, x, yy, yu, xm );
          continue;
        }

  
        ym = yy[i-1]+(xm-x[i-1])*(yy[i-2]-yy[i-1])/(x[i-2]-x[i-1]);
        if (abs(yu[0]-ym) > 2.0*tol*ym + tolmin){ 
          //std::cout << "go to 575" << std::endl; 
          do575(i, x, yy, yu, xm );
          continue;
        }
      } 


      // 560
  
      //std::cout << "560" << std::endl;

      j++;
      uj[j-1] = x[i-1];
      sj[j-1] = yy[i-1];
      if (j > 1){
        sum += 0.5*(yy[i-1]+yl)*(x[i-1]-xl);
        xl = x[i-1];
        yl = yy[i-1];
      }
      --i;
      //std::cout << j << "   " << uj[0] << "  " << uj[1] << "   " << uj[2] << std::endl;
      //std::cout << i << "   " << sj[0] << "  " << sj[1] << "   " << sj[2] << std::endl;
      //std::cout <<"         " << xl << "   " << yl << "     " << sum << std::endl;
      //if (i >= 2){
      //  //std::cout << "go to 530" << std::endl; 
      //  //continue;
      //}
      //else {
      //  break;
      //}
      if ( i < 2 ){
        break;
      }
    }


   // std::cout << std::endl; 
   // std::cout << "580" << std::endl; 
    ++j;
    uj[j-1] = x[0];
    sj[j-1] = yy[0];
    int nmu = j;
    ubar[ie] = 0.0;
    sum = sum + 0.5*(yy[0]+yl)*(x[0]-xl);
    xsi[ie-1] = sum*0.5;
    for (int i = 2; i <= nmu; ++i){
      ubar[ie] = ubar[ie] + 0.5*(uj[i-1]-uj[i-2])*(sj[i-1]+sj[i-2])*(uj[i-1]+uj[i-2]);
      //std::cout << i << "    " <<  ubar[ie] << std::endl;
    }
    ubar[ie] = 0.5*ubar[ie]/sum;
   // std::cout << ubar[ie] << std::endl;


    std::vector<EnergyDistribution> energyDistributionsVec {
        EnergyDistribution( 1.0, { 4 }, { 2 }, { 1e-5, 1.1e+7, 1.147e+7, 3e+7 },
                                                { 0., 2., 4., 6. } ) };

    for (int il = 0; il < nmu; ++il ){
      std::vector<double> scr(1000,0.0);
      u = uj[il];
      //std::cout << "u      " << u << std::endl;
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

      int istart = 1;
      int k = 8;
       
      // 595
      //std::cout << " -------  595 " << std::endl; 
      int iend = nep;
      if (iend-istart >= 153){ iend = istart + 153 - 1; }
      j = k - 1;
      int ib = istart - 1;


      std::vector<double>   EpVec(iend-ib);
      std::vector<double> ProbVec(iend-ib);
      // 596
      do { 
        //std::cout << " -----------------  596 " << std::endl; 
        j  += 2;
        ib += 1;
        scr[j-1] = yu[  2*ib];
        scr[j]   = yu[1+2*ib]*2/sum;
        EpVec[ib-1]   = yu[  2*ib];
        ProbVec[ib-1] = yu[1+2*ib]*2/sum;
        //std::cout << " -----------------  596     " << scr[j-1] << "     " << 
        //                                               scr[j] << std::endl; 
      } while (ib < iend);

      scr.resize(j+1);

      int n2 = j+3;


      std::vector<long> interpolants = {(long) EpVec.size()},
                        boundaries = {2};

      //if (ie==0 and il == 0){std::cout << (EpVec|ranges::view::all) << std::endl;}
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

    //break;

  }

  std::vector< long > boundaries = { (long) esi.size() };
  std::vector< long > interpolants = { 1 };

  LaboratoryAngleEnergy laboratoryAngleEnergy( 
    std::move( boundaries ),
    std::move( interpolants ),
    std::move( angularDistVec) );

  return laboratoryAngleEnergy;

} 


