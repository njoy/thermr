#include "ENDFtk.hpp"
//#include "ENDFtk/section/6.hpp"
#include <range/v3/all.hpp>
//#include "calcem/e_mu_ep/e_mu_ep.h"
#include "inelastic/e_ep_mu.h"
//#include "inelastic/calcem.h"
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
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using Unknown = section::Type<6>::Unknown;

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


template <typename Range, typename Float>
//std::optional<section::Type<6>>  
auto  
  thermr( int matde, int matdp, int nbin, int iinc, int icoh, int iform,
  int natom, const int mtref, Range temps, Float tol, Float emax,
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
      boundCrossSections[1] = freeCrossSections[1]*std::pow((awr2+1)/awr2,2)/natom; 
      // also divided by scr(18) whatever that is
    }
  }
    

  std::vector<double> eftemp(temps.size(),0.0), eftmp2(temps.size(),0.0);

  if (matde == 0){
    for (size_t itemp = 0; itemp < temps.size(); ++itemp){
      eftemp[itemp] = temps[itemp];
    }
  }



  std::vector<file::Type<6>> MF6_vec {};

  for (size_t itemp = 0; itemp < temps.size(); ++itemp){

    std::vector<ReactionProduct> reactionProducts;

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
 
        //std::cout << "A" << std::endl;
        auto incidentEnergies = std::get<0>(out);
        auto totalSCR     = std::get<1>(out);
        auto totalOutput  = std::get<2>(out);
        int n2 = nbin+2;
 
        // Resize the energies and the scattering out to abide by emax
        for (size_t i = 0; i < incidentEnergies.size(); ++i){
          if (incidentEnergies[i] >= emax){
            incidentEnergies.resize(i+2);
            totalSCR.resize(i+2);
            break;
          }
        }


        auto firstSCR = totalSCR[0];
        ThermalScatteringData chunk( incidentEnergies[0], n2, std::move(firstSCR) );
        std::vector<Variant> chunks {chunk};
        for ( size_t j = 1; j < incidentEnergies.size(); ++j){
          auto scratch = totalSCR[j];
          ThermalScatteringData chunk( incidentEnergies[j], n2, std::move(scratch) );
          chunks.push_back(chunk);
        }

        long lep = 1;
        std::vector<long>  boundaries = {(long) incidentEnergies.size()},
                           interpolants = {2};
        ContinuumEnergyAngle continuumChunk( lep, std::move( boundaries ),
                                    std::move( interpolants ),
                                    std::move( chunks ) );
  
        reactionProducts.push_back(ReactionProduct(
            // multiplicity                                      // distribution
          { 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, { 1., 1. }}, continuumChunk ));

        /*
        std::vector<ReactionProduct> products {ReactionProduct(
            // multiplicity                                      // distribution
          { 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, { 1., 1. }}, continuumChunk )};
  
        int jp = 0, lct = 1;
        //return section::Type<6> (mtref, za, awr, jp, lct, std::move(products));
        MF6_vec.push_back(section::Type<6> (mtref, za, awr, jp, lct, std::move(products)));
        */
      }
      else {
        // E mu E' 
 //       auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();
 //       teff = effectiveTemp.effectiveTemperatures()[0]*kb;

        //auto out = e_mu_ep( egrid, tev, tol, lat,  iinc, lasym, alphas, betas, 
        //                    sab, awr, boundCrossSections, teff );
// 
      }
    }



    if (icoh > 0 and icoh <= 10){
      std::cout << "COH" << std::endl;
      if (leapr_MF7.hasMT(2)){
        njoy::ENDFtk::section::Type<7,2> leapr_MT2 = leapr_MF7.MT(2_c);
        auto MT2_law = std::get<CoherentElastic>(leapr_MT2.scatteringLaw());
        int nbragg = MT2_law.numberBraggEdges();
 
        //int jp = 0, lct = 1;
        reactionProducts.push_back(ReactionProduct(
            // multiplicity                                      // distribution
          { 1., 1, -nbragg, 0, {2}, {2}, { 1.e-5, emax }, { 1., 1. }}, Unknown() ));

        /*
        std::vector<ReactionProduct> products {ReactionProduct(
            // multiplicity                                      // distribution
          { 1., 1, -nbragg, 0, {2}, {2}, { 1.e-5, emax }, { 1., 1. }}, Unknown() )};
 
        MF6_vec.push_back(section::Type<6>(mtref, za, awr, jp, lct, 
                                           std::move(products))
                         );
                         */
        //section::Type<6> test(mtref, za, awr, jp, lct, std::move(products));
        //std::string buffer;
        //auto output = std::back_inserter(buffer);
        //test.print(output,425,6);
        //std::cout << buffer << std::endl;
        //std::cout << std::endl;
        //std::cout << std::endl;

      }
      //coh( temp, lat, emax, natom, egrid, tol);
    }
    else if (icoh > 10){
      std::cout << "INCOH" << std::endl;
    }
    int jp = 0, lct = 1;

    std::vector<ReactionProduct> inelasticReaction {reactionProducts[0]};
    std::vector<ReactionProduct>   elasticReaction {reactionProducts[1]};
    //auto elasticReaction = reactionProducts[1];
    //section::Type<6,mtref> inelastic_MF6(mtref, za, awr, jp, lct, std::move(inelasticReaction));
    //section::Type<6,mtref+1> elastic_MF6(mtref, za, awr, jp, lct, std::move(  elasticReaction));
    //MF6_vec.push_back(section::Type<6>(std::move(inelastic_MF6),std::move(elastic_MF6));
    MF6_vec.push_back(section::Type<6>(mtref, za, awr, jp, lct, 
                                       std::move(inelasticReaction)));
    MF6_vec.push_back(section::Type<6>(mtref+1, za, awr, jp, lct, 
                                       std::move(elasticReaction)));



  }

  
  return MF6_vec;
  std::cout << za << awr << lat << lasym << std::endl;
  std::cout << matde+matdp+nbin+iinc+icoh+iform+natom+mtref<< std::endl;
  std::cout << tol + emax + temps[0] << std::endl;
  //std::cout << tev << std::endl;
  std::cout << alphas.size()<< std::endl;
  std::cout << betas.size()<< std::endl;
  std::cout << kb << std::endl;

}


