#include "ENDFtk.hpp"
//#include "ENDFtk/section/6.hpp"
#include <range/v3/all.hpp>
#include "inelastic/e_ep_mu.h"
#include "inelastic/e_mu_ep.h"
#include "coherentElastic/coherentElastic.h"
#include "incoherentElastic/incoherentElastic.h"
#include "generalTools/constants.h"
using namespace njoy;
#include "lipservice.hpp"

using namespace njoy::ENDFtk;
using Tabulated = section::Type< 7, 4 >::Tabulated;
using ContinuumEnergyAngle  = section::Type<6>::ContinuumEnergyAngle;
using LaboratoryAngleEnergy = section::Type<6>::LaboratoryAngleEnergy;
using ThermalScatteringData = section::Type<6>::ContinuumEnergyAngle::ThermalScatteringData;
using Variant = section::Type< 6 >::ContinuumEnergyAngle::Variant;
using ReactionProduct = section::Type< 6 >::ReactionProduct;
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;
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






auto finalTHERMR( nlohmann::json jsonInput, 
  njoy::ENDFtk::syntaxTree::Tape<std::string> leaprTape,
  njoy::ENDFtk::syntaxTree::Tape<std::string> pendfTape ){

  njoy::ENDFtk::file::Type<7> leapr_MF7 = leaprTape.materialNumber(
                      int(jsonInput["matde"])).front().fileNumber(7).parse<7>();
  njoy::ENDFtk::file::Type<1> pendfFile = pendfTape.materialNumber(
                      int(jsonInput["matdp"])).front().fileNumber(1).parse<1>();

  njoy::ENDFtk::section::Type<7,4> leapr_MT4 = leapr_MF7.MT(4_c);
  int nbin   = jsonInput["nbin"];
  int iinc   = jsonInput["iin"];
  int icoh   = jsonInput["icoh"];
  int iform  = jsonInput["iform"];
  int natom  = jsonInput["natom"];
  int mtref  = jsonInput["mtref"];
  std::vector<double> temps = jsonInput["tempr"];
  double tol    = jsonInput["tol"];
  double emax   = jsonInput["emax"];
  
  auto za        = pendfFile.section( 451_c ).ZA();
  auto awr       = leapr_MT4.AWR();
  auto lat       = leapr_MT4.LAT();
  auto lasym     = leapr_MT4.LASYM();
  auto constants = leapr_MT4.constants();
  auto table     = std::get<Tabulated>(leapr_MT4.scatteringLaw());

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
    
  std::vector<file::Type<6>> MF6_vec {};

  for (size_t itemp = 0; itemp < temps.size(); ++itemp){

    std::vector<section::Type<6>> section6Vec {};

    std::vector<double> sab (alphas.size()*betas.size());
    for (size_t ibeta = 0; ibeta < betas.size(); ++ibeta){
      for (size_t ialpha = 0; ialpha < alphas.size(); ++ialpha){
        sab[ialpha*betas.size()+ibeta] = 
            log(table.betas()[ibeta].thermalScatteringValues()[itemp][ialpha]);
      }
    }

    auto temp = temps[itemp];
    
    // compute incoherent inelastic cross sections
    if (iinc != 0){

      if (iform == 0){
        // E E' mu
        auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();
        auto teff = effectiveTemp.effectiveTemperatures()[itemp]*kb;
        auto tev  = temp*kb;
        std::vector<double> initialEnergies;
        for (size_t i = 0; i < egrid.size(); ++i){
          initialEnergies.push_back(egrid[i]);
          if (egrid[i] > emax){ break; }
        }

        auto out = e_ep_mu( initialEnergies, tev, tol, lat,  iinc, lasym, 
                            alphas, betas, sab, awr, boundCrossSections, teff, 
                            nbin, temp );
 
        auto incidentEnergies = std::get<0>(out);
        auto totalSCR     = std::get<1>(out);
        auto totalOutput  = std::get<2>(out);
        int n2 = nbin+2;

        auto firstSCR = totalSCR[0];
        ThermalScatteringData chunk( incidentEnergies[0], n2, std::move(firstSCR) );
        std::vector<Variant> chunks {chunk};
        for ( size_t j = 1; j < incidentEnergies.size(); ++j){
          auto scratch = totalSCR[j];
          ThermalScatteringData chunk(incidentEnergies[j], n2, std::move(scratch));
          chunks.push_back(chunk);
        }

        long lep = 1;
        ContinuumEnergyAngle continuumChunk( lep, 
              {(long) incidentEnergies.size()}, {2}, std::move( chunks ) );
  
        int jp = 0, lct = 1;
        auto pendf_awr       = pendfFile.section( 451_c ).AWR();
        std::vector<ReactionProduct> products = 
          {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                           { 1., 1. }}, continuumChunk )};
        section6Vec.push_back(section::Type<6>(mtref, za, pendf_awr, jp, lct, 
                                           std::move(products)));

      }
      else {
        // E mu E' 
        auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();
        auto teff = effectiveTemp.effectiveTemperatures()[itemp]*kb;

        LaboratoryAngleEnergy labAngleEnergy = e_mu_ep( alphas, betas, sab, iinc, 
            egrid, temp, emax, tol, lat, lasym, awr, boundCrossSections, teff );

        int jp = 0, lct = 1;
        std::vector<ReactionProduct> products = 
          //               zap  awp lip law
          {ReactionProduct({ 1., 1, -1, 7, {2}, {2}, { 1.e-5, emax }, 
                           { 1., 1. }}, std::move(labAngleEnergy))};

        auto pendf_awr       = pendfFile.section( 451_c ).AWR();
        section6Vec.push_back(section::Type<6>(mtref, za, pendf_awr, jp, lct, 
                                           std::move(products)));

      }
    }

    // Coherent Elastic
    if (icoh > 0 and icoh <= 10){
      if (leapr_MF7.hasMT(2)){
        njoy::ENDFtk::section::Type<7,2> leapr_MT2 = leapr_MF7.MT(2_c);
        auto MT2_law = std::get<CoherentElastic>(leapr_MT2.scatteringLaw());
        int nbragg = MT2_law.numberBraggEdges();
        int jp = 0, lct = 1;
        std::vector<ReactionProduct> products {ReactionProduct(
            // multiplicity                                      // distribution
          { 1., 1, -nbragg, 0, {2}, {2}, { 1.e-5, emax }, { 1., 1. }}, Unknown() )};
 
        auto pendf_awr       = pendfFile.section( 451_c ).AWR();
        section6Vec.push_back(section::Type<6>(mtref+1, za, pendf_awr, jp, lct, 
                                           std::move(products)) );

      }
    }
    // Incoherent Elastic
    else if (icoh > 10){

      njoy::ENDFtk::section::Type<7,2> leapr_MT2 = leapr_MF7.MT(2_c);
      auto incoh_law = std::get<IncoherentElastic>(leapr_MT2.scatteringLaw());
      auto chunk = section6Vec[0];

      double debyeWallerFactor = 0.0;

      std::vector<double> temperatures = incoh_law.temperatures(),
                          debyeWaller  = incoh_law.debyeWallerValues();
        
      if ( temp < 0.9*temperatures[0] or 
           temp > 1.1*temperatures[temperatures.size()-1] ){
        std::cout << "oh no put an error here" << std::endl;
      }
      else {
        debyeWallerFactor = interpolate(temperatures, debyeWaller, temp);
      }

      auto law = std::get<ContinuumEnergyAngle>(chunk.products()[0].distribution());

      auto chunks = incoherentElastic( law, nbin, debyeWallerFactor );
    /*
      std::vector<double> esi;
      for (const auto& subsection : law.subsections()){
        esi.push_back(std::get<ThermalScatteringData>(subsection).energy());
      }

      std::vector<std::vector<double>> equiProbCosinesVec;
      std::vector<double> equiprobCosines(nbin+2,1.0);

      for ( const double& E : esi ){
        getIncohElasticDataSingleEnergy( E, debyeWallerFactor, equiprobCosines );
        equiProbCosinesVec.push_back(equiprobCosines);
      }
      //return equiProbCosinesVec;
      int n2 = nbin+2;
      auto firstSCR = equiProbCosinesVec[0];
      ThermalScatteringData chunky( esi[0], n2, std::move(firstSCR) );
      std::vector<Variant> chunks {chunky};
      for ( size_t j = 1; j < esi.size(); ++j){
        auto scratch = equiProbCosinesVec[j];
        ThermalScatteringData chunky( esi[j], n2, std::move(scratch) );
        chunks.push_back(chunky);
      }
      */

      long lep = 1;
      ContinuumEnergyAngle continuumChunk( lep, 
        {(long) chunks.size()}, {2}, std::move(chunks) );

      auto pendf_awr       = pendfFile.section( 451_c ).AWR();
      std::vector<ReactionProduct> iel_products = 
        {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                         { 1., 1. }}, continuumChunk )};
      /*
                         */
      int jp = 0, lct = 1;
      //auto iel_products = incoherentElastic( incoh_law, temp, emax, nbin, section6Vec[0] );
      //auto pendf_awr       = pendfFile.section( 451_c ).AWR();
      section6Vec.push_back(section::Type<6>(mtref+1, za, pendf_awr, jp, lct, 
                                         std::move(iel_products)));

    }
  
    MF6_vec.emplace_back(std::move(section6Vec));


  }

  return MF6_vec;
}


