#include "ENDFtk.hpp"
//#include "ENDFtk/section/6.hpp"
using namespace njoy;
#include <range/v3/all.hpp>
#include "inelastic/e_ep_mu.h"
#include "inelastic/e_mu_ep.h"
#include "coherentElastic/coherentElastic.h"
#include "incoherentElastic/incoherentElastic.h"
#include "generalTools/constants.h"
#include "lipservice.hpp"

using namespace njoy::ENDFtk;
using Tabulated = section::Type< 7, 4 >::TabulatedFunctions;
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


void prepareMF3_inelastic( std::vector<double>& MF3_energies, 
  std::vector<double>& MF3_XS, section::Type<3> MF3_1, double emax,
  std::vector<double> initialEnergies, std::vector<double> xsi ) {

  std::vector<double> desiredEnergies; 
  std::vector<double> finalXS;

  for (auto energy : MF3_1.energies()){
    if (energy >= emax or std::fabs((energy-emax)/emax) < 1e-10){ 
        if (energy > emax*(1+1e-10)){ 
            energy = emax; 
        }
        desiredEnergies.emplace_back(energy); 
        if (std::fabs(energy-emax) < 1e-10*emax){ 
          desiredEnergies.emplace_back(1.00001*energy); 
        }
        break; 
    }
    desiredEnergies.emplace_back(energy);
  }
  for ( const auto& energy : desiredEnergies ){
    finalXS.emplace_back(terp(initialEnergies,xsi,energy,5));
  }
  desiredEnergies.emplace_back(2e7);
  finalXS.emplace_back(0.0);

  if (desiredEnergies[desiredEnergies.size()-2] > emax){ 
      finalXS[desiredEnergies.size()-2] = 0.0; }

  MF3_energies = desiredEnergies;
  MF3_XS       = finalXS;
}




section::Type<3> prepareMF3_cohElastic( std::vector<double> MF3CohElasticEnergies,
  std::vector<double> MF3CohElasticCrossSections, int mtref, int za, 
  double pendf_awr, double temp, std::vector<double>& MF3_energies, 
  std::vector<double>& MF3_XS ){

  std::vector<double> MF3_energies_New = MF3CohElasticEnergies;
  std::vector<double> MF3_XS_New;

  std::vector<long> interpolants = {2},
                    boundaries = {long(MF3CohElasticEnergies.size())};

        // Write out bragg peaks
        section::Type< 3 > MF3_cohEl( mtref+1, za, pendf_awr, temp, 0.0, 0,
                                std::move( boundaries ),
                                std::move( interpolants ),
                                std::move( MF3CohElasticEnergies), 
                                std::move( MF3CohElasticCrossSections) );


        // Write ensure that other cross sections are on same grid as bragg
        int numSec3Energies = MF3_energies.size();

        for (const auto& E : MF3_energies_New){
          int begin = 0, end = 0;

          for ( int j = 0; j < numSec3Energies; ++j){
            if (MF3_energies[j] >= E){ 
              begin = j-2;
              end   = j+3;
              if (j-3 < 0){
                begin = 0;
                end   = 5;
              }

              if (j+3 >= numSec3Energies ){
                end = numSec3Energies - 1;
                begin = end - 5;
              }
              break; 
            }
          }
          int interpOrder = (begin == 0 or end == numSec3Energies - 1) ? 3 : 4;
          int nl          = (begin == 0 or end == numSec3Energies - 1) ? 4 : 5;


          std::vector<double> 
            temporary1 (MF3_energies.begin()+begin,MF3_energies.begin()+end),
            temporary2 (MF3_XS.begin()      +begin,MF3_XS.begin()      +end);

          MF3_XS_New.push_back(terp(temporary1,temporary2,E,interpOrder,nl));
    }
    MF3_XS_New[MF3_XS_New.size()-2] = 0.0;
    MF3_XS_New[MF3_XS_New.size()-1] = 0.0;

    MF3_energies = MF3_energies_New;
    MF3_XS = MF3_XS_New;

  return MF3_cohEl;
}











auto finalTHERMR( const nlohmann::json jsonInput ){
                  //std::ostream& ,//output,
                  //std::ostream& ,//error,
                  //const nlohmann::json&){

  std::cout << std::fixed << std::showpoint;
  std::cout << std::setprecision(10);


  tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                      std::to_string(int(jsonInput["nendf"]))));
  tree::Tape<std::string> pendfTape ( utility::slurpFileToMemory("tape" + 
                                      std::to_string(int(jsonInput["nin"]))));

  file::Type<7> MF7 = leaprTape.material(
                      int(jsonInput["matde"])).front().file(7).parse<7>();
  file::Type<1> MF1 = pendfTape.material(
                      int(jsonInput["matdp"])).front().file(1).parse<1>();
  file::Type<3> MF3 = pendfTape.material(
                      int(jsonInput["matdp"])).front().file(3).parse<3>();
  section::Type<3> MF3_1 = MF3.MT(1_c);

  section::Type<7,4> leapr_MT4 = MF7.MT(4_c);
  int nbin   = jsonInput["nbin"];
  int iinc   = jsonInput["iin"];
  int icoh   = jsonInput["icoh"];
  int natom  = jsonInput["natom"];
  int mtref  = jsonInput["mtref"];
  double tol    = jsonInput["tol"];
  double emax   = jsonInput["emax"];
  
  auto za        = MF1.section( 451_c ).ZA();
  auto awr       = leapr_MT4.AWR();
  auto lat       = leapr_MT4.LAT();
  auto lasym     = leapr_MT4.LASYM();
  auto constants = leapr_MT4.constants();
  auto table     = std::get<Tabulated>(leapr_MT4.scatteringLaw());
  auto pendf_awr = MF1.section( 451_c ).AWR();

  std::vector<double> alphas = table.scatteringFunctions()[0].alphas(),
                      betas  = table.betas();

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
    
  std::vector<Material>        materials {};
  std::vector<DirectoryRecord> index     {};

  for (size_t itemp = 0; itemp < jsonInput["tempr"].size(); ++itemp){

    std::vector<section::Type<6>> section6Vec {};
    std::vector<section::Type<3>> section3Vec {};

    std::vector<double> MF3_energies;
    std::vector<double> MF3_XS;



    std::vector<double> sab (alphas.size()*betas.size());
    for (size_t ibeta = 0; ibeta < betas.size(); ++ibeta){
      for (size_t ialpha = 0; ialpha < alphas.size(); ++ialpha){
        sab[ialpha*betas.size()+ibeta] = 
            log(table.scatteringFunctions()[ibeta].thermalScatteringValues()[itemp][ialpha]);
      }
    }


    double temp = jsonInput["tempr"][itemp];
    
    // compute incoherent inelastic cross sections
    if (iinc != 0){

      std::vector<double> initialEnergies;
      for (size_t i = 0; i < egrid.size(); ++i){
        initialEnergies.push_back(egrid[i]);
        if (egrid[i] > emax){ break; }
      }

      if (jsonInput["iform"] == 0){
        // E E' mu
        auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();
        auto teff = effectiveTemp.effectiveTemperatures()[itemp]*kb;
        //auto tev  = temp*kb;

        auto E_Ep_Mu_output = e_ep_mu( initialEnergies, temp*kb, tol, lat,  iinc, lasym, 
                            alphas, betas, sab, awr, boundCrossSections, teff, 
                            nbin, temp );
 
        auto incidentEnergies = std::get<0>(E_Ep_Mu_output);
        auto totalSCR         = std::get<1>(E_Ep_Mu_output);
        auto totalOutput      = std::get<2>(E_Ep_Mu_output);

        std::vector<Variant> chunks;
        for ( size_t j = 0; j < incidentEnergies.size(); ++j){
          std::vector<double> scratch = totalSCR[j];
          ThermalScatteringData chunk(incidentEnergies[j], nbin+2, std::move(scratch));
          chunks.push_back(chunk);
        }

        ContinuumEnergyAngle continuumChunk(1, {(long) incidentEnergies.size()}, 
                                           {2}, std::move( chunks ) );
  
        std::vector<ReactionProduct> products = 
          {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                           { 1., 1. }}, continuumChunk )};
        section::Type<6> inelastic(mtref, za, pendf_awr, 0, 1, std::move(products));
        index.emplace_back( 6, mtref, inelastic.NC(), 0 );
        section6Vec.push_back( std::move(inelastic) );



        // Do the MF3 stuff
        std::vector<double> xsi;
        for (const auto& entry : totalOutput){
          xsi.emplace_back(entry[0]);
        }

        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_1, emax, initialEnergies, xsi );


      }
      else {
        // E mu E' 
        auto effectiveTempInfo = leapr_MT4.principalEffectiveTemperature();
        auto effectiveTemp = effectiveTempInfo.effectiveTemperatures()[itemp]*kb;

        auto output = e_mu_ep( alphas, betas, sab, iinc, egrid, temp, emax, 
                      tol, lat, lasym, awr, boundCrossSections, effectiveTemp );

        LaboratoryAngleEnergy labAngleEnergy = std::get<0>(output);
        std::vector<double> esi = std::get<1>(output),
                            xsi = std::get<2>(output);

        std::vector<ReactionProduct> products = 
          {ReactionProduct({ 1., 1, -1, 7, {2}, {2}, { 1.e-5, emax }, 
                           { 1., 1. }}, std::move(labAngleEnergy))};

        section::Type<6> inelastic(mtref, za, pendf_awr, 0, 1, std::move(products));
        index.emplace_back( 6, mtref, inelastic.NC(), 0 );
        section6Vec.push_back( std::move(inelastic) );

        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_1, emax, initialEnergies, xsi );


      }
    }

    // Coherent Elastic
    if (icoh > 0 and icoh <= 10){
      if (MF7.hasMT(2)){
        njoy::ENDFtk::section::Type<7,2> leapr_MT2 = MF7.MT(2_c);
        auto MT2_law = std::get<CoherentElastic>(leapr_MT2.scatteringLaw());
        int nbragg = MT2_law.numberBraggEdges();

        std::vector<ReactionProduct> products {ReactionProduct(
          { 1., 1, -nbragg, 0, {2}, {2}, {1.e-5, emax}, {1., 1.}}, Unknown() )};
 
        auto pendf_awr       = MF1.section( 451_c ).AWR();
        section::Type<6> cohElastic(mtref+1, za, pendf_awr, 0, 1, std::move(products));
        index.emplace_back( 6, mtref+1, cohElastic.NC(), 0 );
        section6Vec.push_back( std::move(cohElastic) );

        std::vector<double> cohElasticEnergies = MT2_law.energies(),
                            cohElastic1 = MT2_law.thermalScatteringValues()[0];
        auto out = coh( temp, lat, emax, natom, MF3_energies, tol, 
                        cohElasticEnergies, cohElastic1 );

        std::vector<double> MF3CohElasticEnergies      = std::get<0>(out);
        std::vector<double> MF3CohElasticCrossSections = std::get<1>(out);

        section::Type< 3 > MF3_cohEl = prepareMF3_cohElastic( MF3CohElasticEnergies,
          MF3CohElasticCrossSections, mtref, za, pendf_awr, temp,
          MF3_energies, MF3_XS );

        section3Vec.emplace_back(std::move(MF3_cohEl));


      }
    }
    // Incoherent Elastic
    else if (icoh > 10){

      njoy::ENDFtk::section::Type<7,2> leapr_MT2 = MF7.MT(2_c);
      auto incoh_law = std::get<IncoherentElastic>(leapr_MT2.scatteringLaw());
      auto chunk = section6Vec[0];


      std::vector<double> temperatures = incoh_law.temperatures(),
                          debyeWaller  = incoh_law.debyeWallerValues();
        
      if ( temp < 0.9*temperatures[0] or 
           temp > 1.1*temperatures[temperatures.size()-1] ){
        std::cout << "oh no put an error here" << std::endl;
        // cannot get dwf 
      }

      double debyeWallerFactor = interpolate(temperatures, debyeWaller, temp);

      auto law = std::get<ContinuumEnergyAngle>( 
                                   chunk.reactionProducts()[0].distribution());

      auto chunks = incoherentElastic( law, nbin, debyeWallerFactor );

      long lep = 1;
      ContinuumEnergyAngle continuumChunk( lep, 
        {(long) chunks.size()}, {2}, std::move(chunks) );

      auto pendf_awr       = MF1.section( 451_c ).AWR();
      std::vector<ReactionProduct> iel_products = 
        {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                         { 1., 1. }}, continuumChunk )};

      section::Type<6> incElastic(mtref+1, za, pendf_awr, 0, 1, 
                                  std::move(iel_products));
      index.emplace_back( 6, mtref+1, incElastic.NC(), 0 );
      section6Vec.push_back( std::move(incElastic) );

    }

    // Get MF3 ready to write
    std::vector< long > interpolants = { 2 };
    std::vector< long > boundaries = { long(MF3_energies.size()) };

    section::Type< 3 > MF3( mtref, za, pendf_awr, temp, 0.0, 0,
                            std::move( boundaries ),
                            std::move( interpolants ),
                            std::move( MF3_energies), 
                            std::move( MF3_XS) );

    section3Vec.emplace_back(std::move(MF3));

    file::Type<1> MF1_copy = MF1;
    Material material( int(jsonInput["matdp"]), std::move(MF1_copy), 
                       file::Type<3>(std::move(section3Vec)),
                       file::Type<6>(std::move(section6Vec)) );
    materials.emplace_back(std::move(material));
  
  }

  Tape tape( TapeIdentification(std::string(pendfTape.TPID().text()),1), 
             std::move(materials) );
  std::string buffer;
  auto materialOutput = std::back_inserter( buffer );
  tape.print( materialOutput );
  std::string name = "tape"+std::to_string(int(jsonInput["nout"]));
  std::ofstream out(name);
  out << buffer;
  out.close();

}


