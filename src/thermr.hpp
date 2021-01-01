#include "ENDFtk.hpp"
using namespace njoy;
#include <range/v3/all.hpp>
#include "inelastic/e_ep_mu.h"
#include "inelastic/e_mu_ep.h"
#include "coherentElastic/coherentElastic.h"
#include "incoherentElastic/incoherentElastic.h"
#include "generalTools/constants.h"
#include "write/writingFuncs.h"
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




auto readData( double& awr, int& lat, int& lasym, std::vector<double>& freeCrossSections,
  std::vector<double>& analyticalFunctionTypes, std::vector<double>& atomicWeightRatios,
  std::vector<double>& boundCrossSections, std::vector<double>& alphas, 
  std::vector<double>& betas, std::vector<double>& effectiveTemps, 
  const double& pendf_awr, const nlohmann::json& jsonInput ){
  std::vector<std::vector<double>> sabVec;

  int nendf  = jsonInput["nendf"];
  int matde  = jsonInput["matde"];
  int natom  = jsonInput["natom"];
  int iinc   = jsonInput["iin"];
  if ( iinc == 2 ){ 
    // Use S(a,b) - Read in from LEAPR output
    tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                        std::to_string(nendf)));
    file::Type<7> MF7 = leaprTape.material(
                        matde).front().file(7).parse<7>();

    section::Type<7,4> leapr_MT4 = MF7.MT(4_c);
    awr       = leapr_MT4.AWR();
    lat       = leapr_MT4.LAT();
    lasym     = leapr_MT4.LASYM();
    auto constants = leapr_MT4.constants();
    auto table     = std::get<Tabulated>(leapr_MT4.scatteringLaw());
    freeCrossSections = constants.totalFreeCrossSections();
    analyticalFunctionTypes = constants.analyticalFunctionTypes();
    atomicWeightRatios = constants.atomicWeightRatios();

    alphas = table.scatteringFunctions()[0].alphas();
    betas  = table.betas();

    boundCrossSections.resize(freeCrossSections.size());

    boundCrossSections[0] = freeCrossSections[0]*std::pow((awr+1)/awr,2)/natom;
    if (analyticalFunctionTypes.size() > 0){
      if (analyticalFunctionTypes[0] == 0){
        auto awr2 = atomicWeightRatios[0];
        boundCrossSections[1] = freeCrossSections[1]*std::pow((awr2+1)/awr2,2)/natom; 
        // also divided by scr(18) whatever that is
      }
    }


    auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();

    for (size_t itemp = 0; itemp < jsonInput["tempr"].size(); ++itemp){
      std::vector<double> sab (alphas.size()*betas.size());
      for (size_t ibeta = 0; ibeta < betas.size(); ++ibeta){
        for (size_t ialpha = 0; ialpha < alphas.size(); ++ialpha){
          sab[ialpha*betas.size()+ibeta] = 
              log(table.scatteringFunctions()[ibeta].thermalScatteringValues()[itemp][ialpha]);
        }
      }
      sabVec.push_back(sab);
      auto teff = effectiveTemp.effectiveTemperatures()[itemp]*kb;
      effectiveTemps.push_back(teff);

    }



  }
  else if (iinc == 1) {
    // Free Gas - Do not read in from LEAPR output
    awr = pendf_awr;
    lat = 0;
    lasym = 0;
    boundCrossSections = {std::pow(((pendf_awr +1)/pendf_awr ),2)};
    analyticalFunctionTypes = {0};
    freeCrossSections = {0};

    alphas = {0};
    betas  = { 0, 0.1, 2, 4, 6, 8, 10, 15, 25, 30, 35, 40, 45, 50, 55, 60, 
               65, 70, 75, 80, 100, 120, 140, 160, 180, 200, 250, 300, 350, 
               400, 500, 600, 700, 800, 900, 1000, 1250, 1500, 1750, 2000, 
               2250, 2500, 2750, 3000, 3500 };

    for (size_t itemp = 0; itemp < jsonInput["tempr"].size(); ++itemp){
      std::vector<double> sab(betas.size(),0.0);
      sabVec.push_back(sab);
      effectiveTemps.push_back(1.0);
    }
  }
  return sabVec;
}











namespace njoy {
namespace THERMR {

class THERMR {

public:
void operator()( const nlohmann::json& jsonInput,
                 std::ostream& output,
                 std::ostream& error,
                 const nlohmann::json& ){

  output << "Input arguments:\n" << jsonInput.dump(2) << std::endl;

  tree::Tape<std::string> pendfTape ( utility::slurpFileToMemory("tape" + 
                                      std::to_string(int(jsonInput["nin"]))));

  file::Type<1> MF1 = pendfTape.material(
                      int(jsonInput["matdp"])).front().file(1).parse<1>();
  file::Type<3> MF3 = pendfTape.material(
                      int(jsonInput["matdp"])).front().file(3).parse<3>();
  section::Type<3> MF3_2 = MF3.MT(2_c);

  int matde  = jsonInput["matde"];
  int nendf  = jsonInput["nendf"];
  int nbin   = jsonInput["nbin"];
  int iinc   = jsonInput["iin"];
  int icoh   = jsonInput["icoh"];
  int natom  = jsonInput["natom"];
  int mtref  = jsonInput["mtref"];
  double tol    = jsonInput["tol"];
  double emax   = jsonInput["emax"];
  
  auto za = MF1.section( 451_c ).ZA();
  auto pendf_awr = MF1.section( 451_c ).AWR();

  double awr;
  int lat, lasym;
  std::vector<double> freeCrossSections, analyticalFunctionTypes, atomicWeightRatios;
  std::vector<double> boundCrossSections, alphas, betas, effectiveTemps;



   auto sabVec = readData( awr, lat, lasym, freeCrossSections, 
                           analyticalFunctionTypes, atomicWeightRatios,
                           boundCrossSections, alphas, betas, effectiveTemps,
                           pendf_awr, jsonInput );


  std::vector<Material>        materials {};
  std::vector<DirectoryRecord> index     {};

  for (size_t itemp = 0; itemp < jsonInput["tempr"].size(); ++itemp){

    std::vector<section::Type<6>> section6Vec {};
    std::vector<section::Type<3>> section3Vec {};

    std::vector<double> MF3_energies, MF3_XS;

    auto sab = sabVec[itemp];


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

        auto effectiveTemp = effectiveTemps[itemp];
        auto E_Ep_Mu_output = e_ep_mu( initialEnergies, temp*kb, tol, lat,  iinc, lasym, 
                            alphas, betas, sab, awr, boundCrossSections, effectiveTemp, 
                            nbin, temp );
 
        auto incidentEnergies = std::get<0>(E_Ep_Mu_output);
        auto totalSCR         = std::get<1>(E_Ep_Mu_output);
        auto totalOutput      = std::get<2>(E_Ep_Mu_output);



        section::Type<6> inelastic = prepareMF6_inelastic( incidentEnergies, 
                                    totalSCR, nbin, emax, mtref, za, pendf_awr);


        index.emplace_back( 6, mtref, inelastic.NC(), 0 );
        section6Vec.push_back( std::move(inelastic) );

        // Prepare and write inelastic MF3
        std::vector<double> xsi;
        for (const auto& entry : totalOutput){
          xsi.emplace_back(entry[0]);
        }

        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_2, emax, initialEnergies, xsi, iinc );


      }
      else {
        // E mu E' 
        auto output = e_mu_ep( alphas, betas, sab, iinc, egrid, temp, emax, 
                      tol, lat, lasym, awr, boundCrossSections, effectiveTemps[itemp] );

        LaboratoryAngleEnergy labAngleEnergy = std::get<0>(output);
        std::vector<double>   xsi            = std::get<2>(output);

        std::vector<ReactionProduct> products = 
          {ReactionProduct({ 1., 1, -1, 7, {2}, {2}, { 1.e-5, emax }, 
                           { 1., 1. }}, std::move(labAngleEnergy))};

        section::Type<6> inelastic(mtref, za, pendf_awr, 0, 1, std::move(products));
        index.emplace_back( 6, mtref, inelastic.NC(), 0 );
        section6Vec.push_back( std::move(inelastic) );

        // Prepare and write inelastic MF3
        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_2, emax, initialEnergies, xsi, iinc );


      }
    }

    // Coherent Elastic
    if (icoh > 0 and icoh <= 10 and matde != 0){
      tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                          std::to_string(nendf)));
      file::Type<7> MF7 = leaprTape.material(matde).front().file(7).parse<7>();

      if (MF7.hasMT(2)){
        njoy::ENDFtk::section::Type<7,2> leapr_MT2 = MF7.MT(2_c);
        auto MT2_law = std::get<CoherentElastic>(leapr_MT2.scatteringLaw());
        int nbragg = MT2_law.numberBraggEdges();

        std::vector<ReactionProduct> products {ReactionProduct(
          { 1., 1, -nbragg, 0, {2}, {2}, {1.e-5, emax}, {1., 1.}}, Unknown() )};
 
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
          MF3CohElasticCrossSections, mtref, za, pendf_awr, temp, MF3_energies, 
          MF3_XS );
 
        section3Vec.emplace_back(std::move(MF3_cohEl));
  
      }
    }
    // Incoherent Elastic
    else if (icoh > 10){

      if (matde != 0){
        tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                            std::to_string(nendf)));

        file::Type<7> MF7 = leaprTape.material(matde).front().file(7).parse<7>();

        njoy::ENDFtk::section::Type<7,2> leapr_MT2 = MF7.MT(2_c);
        auto incoh_law = std::get<IncoherentElastic>(leapr_MT2.scatteringLaw());
        auto chunk = section6Vec[0];
  
  
        std::vector<double> temperatures = incoh_law.temperatures(),
                            debyeWaller  = incoh_law.debyeWallerValues();
          
        if ( temp < 0.9*temperatures[0] or 
             temp > 1.1*temperatures[temperatures.size()-1] ){
          error << "Cannot obtain Debye-Waller factor " <<
                   "- Temperature out of interpolation range" << std::endl;
        }
  
        double debyeWallerFactor = interpolate(temperatures, debyeWaller, temp);
  
        auto law = std::get<ContinuumEnergyAngle>( 
                                     chunk.reactionProducts()[0].distribution());
  
        auto chunks = incoherentElastic( law, nbin, debyeWallerFactor );
  
        ContinuumEnergyAngle continuumChunk( long(1), 
          {(long) chunks.size()}, {2}, std::move(chunks) );
  
        std::vector<ReactionProduct> iel_products = 
          {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                           { 1., 1. }}, continuumChunk )};
  
        section::Type<6> incElastic(mtref+1, za, pendf_awr, 0, 1, 
                                    std::move(iel_products));
        index.emplace_back( 6, mtref+1, incElastic.NC(), 0 );
        section6Vec.push_back( std::move(incElastic) );

      }
    }

    // Get MF3 ready to write
    section::Type< 3 > MF3( mtref, za, pendf_awr, temp, 0.0, 0,
                            std::vector<long>(1,long(MF3_energies.size())),
                            std::vector<long>(1,2),
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

};
}
}
