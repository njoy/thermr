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


template <typename Sec3, typename Sec6>
void write2( Sec3 sec3_vec, Sec6 sec6_vec, std::vector<double> temperatures,
        nlohmann::json json){
  tree::Tape<std::string> tape (utility::slurpFileToMemory("tape" + 
                                std::to_string(int(json["nin"]))));

  std::string newtape;
  auto output = std::back_inserter( newtape );
  tape.TPID().print( output, 1, 0, 0 );

  double currentTemperature;
  bool addedMF3 = false;
  bool addedMF6 = false;
  //std::cout << "added mf3? mf6?      " << addedMF3 << "   " << addedMF6 << std::endl;

  for ( const auto& material : tape.materials() ) {

    //std::cout << "starting material" << std::endl;
    // std::cout << "-----------------" << std::endl;
    int MAT = material.MAT();

    std::string newmaterial;
    std::vector< DirectoryRecord > newindex;
    std::string description_451;
    double za,awr;

    for ( const auto& file : material.files() ) {
      std::cout << "starting new file in material" << std::endl;
      switch ( file.MF() ) {
        case 1 : { 
          std::cout << "Doing 1" << std::endl;
          section::Type<1,451> section_451 = file.parse(1_c).section(451_c);
          currentTemperature = section_451.TEMP();
          description_451 = section_451.description();
          za  = section_451.ZA();
          awr = section_451.AWR();
          continue;
        }
        case 3 : { 
          std::cout << "Doing 3" << std::endl;
          std::vector<section::Type<3>> sections = file.parse(3_c).sections();
          for (size_t i = 0; i < temperatures.size(); ++i){ 
            if (std::fabs(currentTemperature-temperatures[i]) < 1e-2){
              for (auto& sec3_chunk : sec3_vec[i]){sections.emplace_back(sec3_chunk);}
              addedMF3 = true; 
              break;
            }
          }
          file::Type< 3 > newfile( std::move( sections ) );
          for ( const auto& section : newfile.sections() ) {
            newindex.emplace_back( 3, section.MT(), section.NC(), 0 );
          }
          output = std::back_inserter( newmaterial );
          newfile.print( output, MAT );
          continue;
        }
        case 6: { std::cout << "Doing 6" << std::endl;
          std::vector< section::Type< 6 > > sections = file.parse( 6_c ).sections();
          for (size_t i = 0; i < temperatures.size(); ++i){ 
            if (std::fabs(currentTemperature-temperatures[i]) < 1e-2){
              for (auto& sec6_chunk : sec6_vec[i]){sections.emplace_back(sec6_chunk);}
              addedMF6 = true; 
              break;
            }
          }
          file::Type< 6 > newfile( std::move( sections ) );
          for ( const auto& section : newfile.sections() ) {
            newindex.emplace_back( 6, section.MT(), section.NC(), 0 );
          }
          output = std::back_inserter( newmaterial );
          newfile.print( output, MAT );
          continue;
        }
        default : { std::cout << "Doing other" << std::endl;
          std::string string( file.buffer() );
          newmaterial += string;
          for ( const auto& section : file.sections() ) {
            newindex.emplace_back(file.MF(),section.MT(),ranges::count(section.buffer(),'\n'),0);
          }
          continue;
        }
      }
    }
    std::cout << "added mf3? mf6?      " << addedMF3 << "   " << addedMF6 << std::endl;

    if (not addedMF3){ std::cout << "Doing 3 (later)" << std::endl;
      std::vector<section::Type<3>> sections {};
      for (size_t i = 0; i < temperatures.size(); ++i){ 
        if (std::fabs(currentTemperature-temperatures[i]) < 1e-2){
          for (auto& sec3_chunk : sec3_vec[i]){
            sections.emplace_back(sec3_chunk);
          }
          break;
        }
      }

      file::Type<3> newfile(std::move(sections));
      for (const auto& section : newfile.sections()){
        newindex.emplace_back(3,section.MT(),section.NC(),0);
      }
      output = std::back_inserter(newmaterial);
      newfile.print(output,MAT);
    }
    
    if (not addedMF6){ std::cout << "Doing 6 (later)" << std::endl;
      std::vector<section::Type<6>> sections {};
      for (size_t i = 0; i < temperatures.size(); ++i){ 
        if (std::fabs(currentTemperature-temperatures[i]) < 1e-2){
          for (auto& sec6_chunk : sec6_vec[i]){
            sections.emplace_back(sec6_chunk);
          }
          break;
        }
      }

      file::Type<6> newfile(std::move(sections));
      for (const auto& section : newfile.sections()){
        newindex.emplace_back(6,section.MT(),section.NC(),0);
      }
      output = std::back_inserter(newmaterial);
      newfile.print(output,MAT);
    }

    section::Type<1,451> new451 (za, awr, 2, 0, 0, 0, 0, 0, 0, 0, 6, 1, 2e7, 
                                 0, 10, 8, currentTemperature, 1, 
                                 description_451, std::move(newindex) );
    file::Type<1> newMF1 (std::move(new451));

    std::string buffer_451;
    auto output2 = std::back_inserter( buffer_451 );
    newMF1.print( output2, MAT );

    newtape += buffer_451 + newmaterial;

    output = std::back_inserter( newtape );
    MEND().print( output );
  }

  output = std::back_inserter( newtape );
  TEND().print( output );

  std::string name = "tape"+std::to_string(int(json["nout"]));
  std::ofstream out(name);
  out << newtape;
 
}




namespace njoy {
namespace THERMR {

class THERMR {

public:
void operator()( const nlohmann::json& json, std::ostream& output,
                 std::ostream& error, const nlohmann::json& ){

  tree::Tape<std::string> pendf (utility::slurpFileToMemory("tape" + 
                                      std::to_string(int(json["nin"]))));

  file::Type<1> MF1 = pendf.material(int(json["matdp"])).front().file(1).parse<1>();
  file::Type<3> MF3 = pendf.material(int(json["matdp"])).front().file(3).parse<3>();
  section::Type<3> MF3_2 = MF3.MT(2_c);

  int matde  = json["matde"], nendf  = json["nendf"],
      nbin   = json["nbin" ], iinc   = json["iin"  ],
      icoh   = json["icoh" ], natom  = json["natom"],
      mtref  = json["mtref"], lat, lasym;

  double tol       = json["tol"  ],
         emax      = json["emax" ],
         za        = MF1.section( 451_c ).ZA(),
         pendf_awr = MF1.section( 451_c ).AWR(),
         awr;

  std::vector<double> initialEnergies;
  int i = 0;
  do { initialEnergies.push_back(egrid[i]); }
  while (egrid[i++] <= emax);

  // Read in the data from LEAPR output [or set to default values if free gas]
  std::vector<double> freeXS, analyticalFunctionTypes, awrVec, boundXS, alphas, 
                      betas, effectiveTemps;
  auto sabVec = readLEAPRData( awr, lat, lasym, freeXS, analyticalFunctionTypes, 
       awrVec, boundXS, alphas, betas, effectiveTemps, pendf_awr, json );


  //std::vector<Material>        materials {};

  std::vector<std::vector<section::Type<3>>> Sec3Vec {};
  std::vector<std::vector<section::Type<6>>> Sec6Vec {};
  std::vector<double> temperatures;;

  for (size_t itemp = 0; itemp < json["tempr"].size(); ++itemp){
    double temp = json["tempr"][itemp];
    temperatures.emplace_back(temp);
    std::vector<DirectoryRecord> index     {};

    std::vector<section::Type<6>> section6Vec {};
    std::vector<section::Type<3>> section3Vec {};
    std::vector<double> MF3_energies, MF3_XS, sab = sabVec[itemp],
                        inelasticEnergies;
    
    // compute incoherent inelastic cross sections
    if (iinc != 0){

      if (json["iform"] == 0){
        // E E' mu
        auto E_Ep_Mu_output = e_ep_mu( initialEnergies, temp*kb, tol, lat,  iinc, 
             lasym, alphas, betas, sab, awr, boundXS, effectiveTemps[itemp], 
             nbin, temp );
 
        prepareMF6_E_Ep_Mu( E_Ep_Mu_output, nbin, emax, mtref, za, pendf_awr, 
                            index, section6Vec);

        // Prepare and write inelastic MF3
        std::vector<double> xsi;
        inelasticEnergies = std::get<0>(E_Ep_Mu_output);
        auto totalOutput = std::get<2>(E_Ep_Mu_output);
        for (const auto& entry : totalOutput){ xsi.emplace_back(entry[0]); }

        prepareMF3_inelastic( MF3_energies, MF3_XS,  MF3_2, emax, 
                       initialEnergies, xsi, iinc, mtref, za, pendf_awr, temp );


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
                       initialEnergies, xsi, iinc, mtref, za, pendf_awr, temp );


      }
    }

    // Coherent Elastic
    if (icoh > 0 and icoh <= 10 and matde != 0){

      // MF6 Coherent Elastic
      auto cohElastic = prepareMF6_cohElastic( nendf, matde, emax, mtref, za, 
                                              pendf_awr, index, section6Vec );

      // Pull out Bragg Edges and create MF3
      auto out = coh( temp, lat, emax, natom, MF3_energies, tol, 
              cohElastic.energies(), cohElastic.thermalScatteringValues()[0] );

      section3Vec.emplace_back(prepareMF3_cohElastic( out, mtref, za, 
                               pendf_awr, temp, MF3_energies, MF3_XS, index ));
    } 

    else if (icoh > 10 and matde != 0){
      // Incoherent Elastic
      double debyeWallerFactor = prepareMF6_incohElastic( nendf, matde, temp, 
              emax, nbin, mtref, za, pendf_awr, section6Vec, error, index );

      section3Vec.emplace_back(prepareMF3_incohElastic( inelasticEnergies,
        MF3_energies, boundXS[0], debyeWallerFactor, mtref, natom, za, 
        pendf_awr, temp, index));

    } 


    section::Type<3> MF3( mtref, za, pendf_awr, temp, 0.0, 0,
      std::vector<long>(1,long(MF3_energies.size())), std::vector<long>(1,2), 
      std::move( MF3_energies ), std::move( MF3_XS ) );

    index.emplace_back( 3, mtref, MF3.NC(), 0 );

    section3Vec.emplace_back(std::move(MF3));


    Sec3Vec.emplace_back(std::move(section3Vec));
    Sec6Vec.emplace_back(std::move(section6Vec));


    //file::Type<1> MF1_copy = MF1;
    //section::Type<1,451> section_451 = MF1.MT(451_c);
    //for (const auto& val : section_451.index()){
    //  index.emplace_back(val.MF(), val.MT(), val.NC(), 0);
    //}

    //section::Type<1,451> new451 (za, awr, 2, 0, 0, 0, 0, 0, 0, 0, 6, 1, 2e7, 
    //                             0, 10, 8, temp, 1, section_451.description(), 
    //                             std::move(index) );
    //file::Type<1> newMF1 (std::move(new451));

    //Material material( int(json["matdp"]), std::move(newMF1), 
    //                   file::Type<3>(std::move(section3Vec)),
    //                   file::Type<6>(std::move(section6Vec)) );
    //materials.emplace_back(std::move(material));
  
  }

  //writeTape( materials, pendf, int(json["nout"]) );
  //write(section3Vec,section6Vec,json);
  write2(Sec3Vec,Sec6Vec,temperatures,json);

  return;

  output << "Input arguments:\n" << json.dump(2) << std::endl;

}

};
}
}







