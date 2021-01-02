#include "ENDFtk.hpp"
using namespace njoy;
#include <range/v3/all.hpp>
#include "inelastic/e_ep_mu.h"
#include "inelastic/e_mu_ep.h"
#include "coherentElastic/coherentElastic.h"
#include "incoherentElastic/incoherentElastic.h"
#include "generalTools/constants.h"
#include "readWrite/readingFuncs.h"
#include "readWrite/writingFuncs.h"
#include "lipservice.hpp"

using namespace njoy::ENDFtk;
using ContinuumEnergyAngle  = section::Type<6>::ContinuumEnergyAngle;
using LaboratoryAngleEnergy = section::Type<6>::LaboratoryAngleEnergy;
using ThermalScatteringData = section::Type<6>::ContinuumEnergyAngle::ThermalScatteringData;
using Variant = section::Type< 6 >::ContinuumEnergyAngle::Variant;
using ReactionProduct = section::Type< 6 >::ReactionProduct;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;

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

  int matde   = jsonInput["matde"], nendf   = jsonInput["nendf"],
      nbin    = jsonInput["nbin" ], iinc    = jsonInput["iin"  ],
      icoh    = jsonInput["icoh" ], natom   = jsonInput["natom"],
      mtref   = jsonInput["mtref"], lat, lasym;
  double tol  = jsonInput["tol"  ],
         emax = jsonInput["emax" ],
         za   = MF1.section( 451_c ).ZA(),
         pendf_awr = MF1.section( 451_c ).AWR(),
         awr;


  // Read in the data from LEAPR output [or set to default values if free gas]
  std::vector<double> freeXS, analyticalFunctionTypes, awrVec, boundXS, alphas, 
                      betas, effectiveTemps;
  auto sabVec = readLEAPRData( awr, lat, lasym, freeXS, analyticalFunctionTypes, 
       awrVec, boundXS, alphas, betas, effectiveTemps, pendf_awr, jsonInput );


  std::vector<Material>        materials {};
  std::vector<DirectoryRecord> index     {};

  for (size_t itemp = 0; itemp < jsonInput["tempr"].size(); ++itemp){
    double temp = jsonInput["tempr"][itemp];

    std::vector<section::Type<6>> section6Vec {};
    std::vector<section::Type<3>> section3Vec {};
    std::vector<double> MF3_energies, MF3_XS, sab = sabVec[itemp];
    
    // compute incoherent inelastic cross sections
    if (iinc != 0){
      std::vector<double> initialEnergies;
      for (size_t i = 0; i < egrid.size(); ++i){
        initialEnergies.push_back(egrid[i]);
        if (egrid[i] > emax){ break; }
      }

      if (jsonInput["iform"] == 0){
        // E E' mu
        auto E_Ep_Mu_output = e_ep_mu( initialEnergies, temp*kb, tol, lat,  iinc, 
             lasym, alphas, betas, sab, awr, boundXS, effectiveTemps[itemp], 
             nbin, temp );
 
        prepareMF6_E_Ep_Mu( E_Ep_Mu_output, nbin, emax, mtref, za, pendf_awr, 
                            index, section6Vec);

        // Prepare and write inelastic MF3
        std::vector<double> xsi;
        auto totalOutput = std::get<2>(E_Ep_Mu_output);
        for (const auto& entry : totalOutput){ xsi.emplace_back(entry[0]); }
        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_2, emax, initialEnergies, xsi, iinc );

      }
      else {
        // E mu E' 
        auto output = e_mu_ep( alphas, betas, sab, iinc, egrid, temp, emax, 
                      tol, lat, lasym, awr, boundXS, effectiveTemps[itemp] );

        LaboratoryAngleEnergy labAngleEnergy = std::get<0>(output);
        std::vector<double>   xsi            = std::get<2>(output);

        // Prepare and write inelastic MF6
        prepareMF6_E_Mu_Ep( emax, labAngleEnergy, mtref, za, pendf_awr, index, 
                            section6Vec );

        // Prepare and write inelastic MF3
        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_2, emax, 
                              initialEnergies, xsi, iinc );
      }
    }

    // Coherent Elastic
    if (icoh > 0 and icoh <= 10 and matde != 0){
      tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                          std::to_string(nendf)));
      file::Type<7> MF7 = leaprTape.material(matde).front().file(7).parse<7>();
      if (MF7.hasMT(2)){   

        // MF6 Coherent Elastic
        auto cohElasticInfo = prepareMF6_cohElastic( MF7, emax, mtref, za, 
                                                pendf_awr, index, section6Vec );
        // Pull out Bragg Edges and create MF3
        std::vector<double> 
          cohElasticEnergies = cohElasticInfo.energies(),
          cohElasticWeights  = cohElasticInfo.thermalScatteringValues()[0];

        auto out = coh( temp, lat, emax, natom, MF3_energies, tol, 
                        cohElasticEnergies, cohElasticWeights );

        std::vector<double> MF3CohElasticEnergies = std::get<0>(out),
                            MF3CohElasticXS       = std::get<1>(out);

        section3Vec.emplace_back(prepareMF3_cohElastic( out, mtref, za, 
                                       pendf_awr, temp, MF3_energies, MF3_XS ));
      } 
    } 

    else if (icoh > 10 and matde != 0){
      // Incoherent Elastic
      prepareMF6_incohElastic( nendf, matde, temp, emax, nbin, mtref, za, 
                               pendf_awr, section6Vec, error, index );
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
