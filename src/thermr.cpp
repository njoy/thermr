#include "ENDFtk.hpp"
#include <range/v3/all.hpp>
//#include "calcem/calcem_util/e_mu_ep.h"
#include "calcem/calcem_util/e_ep_mu.h"
#include "calcem/calcem.h"
#include "coh/coh.h"

using namespace njoy::ENDFtk;
using ScatteringLaw = section::Type< 7, 4 >::ScatteringLaw;
using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
using AnalyticalFunctions = section::Type< 7, 4 >::AnalyticalFunctions;
using Tabulated = section::Type< 7, 4 >::Tabulated;
using ScatteringFunction = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
using EffectiveTemperature = section::Type< 7, 4 >::EffectiveTemperature;
using MF7 = njoy::ENDFtk::file::Type<7>;
using ContinuumEnergyAngle  = section::Type<6>::ContinuumEnergyAngle;
using ThermalScatteringData = section::Type<6>::ContinuumEnergyAngle::ThermalScatteringData;
using Variant = section::Type< 6 >::ContinuumEnergyAngle::Variant;
using ReactionProduct = section::Type< 6 >::ReactionProduct;

/*
std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5,
   1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 0.000253, 0.000297, 0.000350, 0.00042, 
   0.000506, 0.000615, 0.00075, 0.00087, 0.001012, 0.00123, 0.0015, 
   0.0018, 0.00203, 0.002277, 0.0026, 0.003, 0.0035, 0.004048, 0.0045, 
   0.005, 0.0056, 0.006325, 0.0072, 0.0081, 0.009108, 0.01, 0.01063, 
   0.0115, 0.012397, 0.0133, 0.01417, 0.015, 0.016192, 0.0182, 0.0199, 
   0.020493, 0.0215, 0.0228, 0.0253, 0.028, 0.030613, 0.0338, 0.0365, 
   0.0395, 0.042757, 0.0465, 0.050, 0.056925, 0.0625, 0.069, 0.075, 
   0.081972, 0.09, 0.096, 0.1035, 0.111573, 0.120, 0.128, 0.1355, 0.145728, 
   0.160, 0.172, 0.184437, 0.20, 0.2277, 0.2510392, 0.2705304, 0.2907501, 
   0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 0.5032575, 0.56, 
   0.625, 0.70, 0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 1.70, 1.855,
   2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 4.90, 5.35,
   5.85, 6.40, 7.00, 7.65, 8.40, 9.15, 9.85, 10.00 };
   */
std::vector<double> egrid {1.e-5,1.78e-5,10.0};




template <typename Range, typename Float>
std::optional<section::Type<6>>  
  thermr( int matde, int matdp, int nbin, int iinc, int icoh, int iform,
  int natom, int mtref, Range temps, Float tol, Float emax,
  MF7 leapr_MF7 ){

  double kb = 8.6173303e-5;
  std::cout.precision(15);

  njoy::ENDFtk::section::Type<7,4> leapr_MT4 = leapr_MF7.MT(4_c);
  auto za  = leapr_MT4.ZA();
  auto awr = leapr_MT4.AWR();
  auto lat = leapr_MT4.LAT();
  auto lasym = leapr_MT4.LASYM();

  auto constants = leapr_MT4.constants();
  auto table = std::get<Tabulated>(leapr_MT4.scatteringLaw());

  std::vector<double> alphas = table.betas()[0].alphas(),
                      betas(table.numberBetas() );


  for (int ibeta = 0; ibeta < table.numberBetas(); ++ibeta){ 
      auto value = table.betas()[ibeta];
      betas[ibeta] = value.beta();
  }

  auto freeCrossSections = constants.totalFreeCrossSections();
  std::vector<double> boundCrossSections(freeCrossSections.size(),0.0);
  boundCrossSections[0] = freeCrossSections[0]*std::pow((awr+1)/awr,2)/natom;
  if (constants.analyticalFunctionTypes().size() > 0){
    if (constants.analyticalFunctionTypes()[0] == 0){
      auto awr2 = constants.atomicWeightRatios()[0];
      boundCrossSections[1] = freeCrossSections[1]*std::pow((awr2+1)/awr2,2)/natom; // also divided by scr(18) whatever that is
    }
  }
    

  std::vector<double> eftemp(temps.size(),0.0), eftmp2(temps.size(),0.0);

  if (matde == 0){
    for (size_t itemp = 0; itemp < temps.size(); ++itemp){
      eftemp[itemp] = temps[itemp];
    }
  }

  for (size_t itemp = 0; itemp < temps.size(); ++itemp){

    std::vector<double> sab (alphas.size()*betas.size());
 
    for (size_t ibeta = 0; ibeta < betas.size(); ++ibeta){
      for (size_t ialpha = 0; ialpha < alphas.size(); ++ialpha){
        sab[ialpha*betas.size()+ibeta] = 
            log(table.betas()[ibeta].thermalScatteringValues()[itemp][ialpha]);
      }
    }

    auto temp = temps[itemp];
    auto tev  = temp*kb;
    auto teff = eftemp[itemp];

    
    // compute incoherent inelastic cross sections
    if (iinc != 0){

       if (iform == 0){
         // E E' mu
         auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();
         teff = effectiveTemp.effectiveTemperatures()[0]*kb;

         auto out = e_ep_mu( egrid, tev, tol, lat,  iinc, lasym, alphas, betas, 
                             sab, awr, boundCrossSections, teff, nbin, temp );

         auto incidentEnergies = std::get<0>(out);
         auto totalSCR     = std::get<1>(out);
         auto totalOutput  = std::get<2>(out);
         int n2 = nbin+2;

         auto firstSCR = totalSCR[0];
         ThermalScatteringData chunk( incidentEnergies[0], n2, std::move(firstSCR) );
         std::vector<Variant> chunks {chunk};
         //chunks.resize(incidentEnergies.size());
         for ( size_t j = 1; j < incidentEnergies.size(); ++j){
           auto scratch = totalSCR[j];
           ThermalScatteringData chunk( incidentEnergies[j], n2, std::move(scratch) );
           chunks.push_back(chunk);
         }

      //int lang = 3;
      long lep = 1;
      std::vector<long>   boundaries = {3},
                        interpolants = {2};
      ContinuumEnergyAngle continuumChunk( lep, std::move( boundaries ),
                                  std::move( interpolants ),
                                  std::move( chunks ) );
      

      std::vector<ReactionProduct> products {ReactionProduct(
        // multiplicity
        { 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, { 1., 1. }},
        // distribution
         continuumChunk )};

         int jp = 0, lct = 1;
         return section::Type<6> (mtref, za, awr, jp, lct, std::move(products));
       }
       else {
         // E mu E' 
       }
    }



    if (icoh > 0 and icoh <= 10){
      std::cout << "COH" << std::endl;
      coh( temp, lat, emax, natom, egrid, tol);
    }
    else if (icoh > 10){
      std::cout << "INCOH" << std::endl;
    }

  }

  
  return {};
  std::cout << za << awr << lat << lasym << std::endl;
  std::cout << matde+matdp+nbin+iinc+icoh+iform+natom+mtref<< std::endl;
  std::cout << tol + emax + temps[0] << std::endl;
  //std::cout << tev << std::endl;
  std::cout << alphas.size()<< std::endl;
  std::cout << betas.size()<< std::endl;
  std::cout << kb << std::endl;

}


