#include "ENDFtk.hpp"
#include <range/v3/all.hpp>
#include "calcem/calcem_util/e_mu_ep.h"
#include "calcem/calcem_util/e_ep_mu.h"

using namespace njoy::ENDFtk;

using ScatteringLaw = section::Type< 7, 4 >::ScatteringLaw;
using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
using AnalyticalFunctions = section::Type< 7, 4 >::AnalyticalFunctions;
using Tabulated = section::Type< 7, 4 >::Tabulated;
using ScatteringFunction = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
using EffectiveTemperature = section::Type< 7, 4 >::EffectiveTemperature;

using namespace njoy::ENDFtk;
int thermr( int mat, std::string fileName, int iform ){


  std::string contents = njoy::utility::slurpFileToMemory(fileName);
  njoy::ENDFtk::syntaxTree::Tape< std::string > tape = 
    njoy::ENDFtk::syntaxTree::Tape< std::string >( contents );

  decltype(auto) materials = tape.MAT( mat ); 
  // MAT( ... ) returns all materials with that nuber 
  // - a requirement due to multitemperature tapes
  decltype(auto) material = materials.front(); // this is the first material, the one you need
  decltype(auto) file7 = material.MF(7).parse< 7 >();
     
  decltype(auto) mt4 = file7.section( 4_c ); // MT( 4_c ) would achieve the same thing
  // mt4 is now a fully parser MF7 MT4 from which you can now request the data in it:

  auto LAT = mt4.temperatureOption(); // LAT()
  auto LASYM = mt4.symmetryOption(); // LASYM
  auto constants = mt4.constants(); // the ScatteringLawConstants instance inside mt4

  auto totalFreeXS = constants.totalFreeCrossSections();
  auto atomicWeightRatios = constants.atomicWeightRatios();
  std::cout << LAT << std::endl;
  std::cout << LASYM << std::endl;
  std::cout << (totalFreeXS|ranges::view::all) << std::endl;
  std::cout << (atomicWeightRatios|ranges::view::all) << std::endl;

  auto law = mt4.scatteringLaw(); // a variant of either a tabulated law or an analytical law
  auto table = std::get< Tabulated >( law );
  auto value = table.betas()[0];
  auto nbeta = table.numberBetas();
  std::vector<double> betas (nbeta);
  for (size_t i = 0; i < (unsigned) nbeta; ++i){
    betas[i] = table.betas()[i].beta();
  }
  std::vector<double> alphas = table.betas()[0].alphas();

  std::cout << (alphas|ranges::view::all) << std::endl;
  std::cout << (betas|ranges::view::all) << std::endl;



  //std::cout << value.beta() << std::endl;
  //std::cout << value.LT() << std::endl;
  //std::cout << value.temperatureDependenceFlag() << std::endl;
  //std::cout << value.NT() << std::endl;
  //std::cout << value.numberTemperatures() << std::endl;

  //std::cout << betas << std::endl;

  //std::cout << typeid(betas).name() << std::endl;
  //std::cout << (betas|ranges::view::all) << std::endl;


  //auto table = ;


  if (iform == 0){
    // E E' mu
    //auto out = //std::make_tuple(eVec,total_SCR,total_OutputData);
    //  e_ep_mu( eVec, tev, tol, lat, iinc, lasym, alphas, betas, sab, az, sigma_b, sigma_b2, teff, nbin, temp );
    //auto eVec = std::get<0>(out);
    //auto total_SCR = std::get<1>(out);
    //auto total_OutputData = std::get<2>(out);



  }
  else {
    // E mu E'
  }

  /*
chunk.LLN() );
chunk.sabStorageType() );
chunk.NI() );
chunk.numberConstants() );
chunk.NS() );
chunk.numberNonPrincipalScatterers() );
chunk.epsilon() ) );
chunk.upperEnergyLimit() ) );
chunk.totalFreeCrossSections().size() );
chunk.totalFreeCrossSections()[0] ) );
chunk.totalFreeCrossSections()[1] ) );
chunk.totalFreeCrossSections()[2] ) );
chunk.atomicWeightRatios().size() );
chunk.atomicWeightRatios()[0] ) );
chunk.atomicWeightRatios()[1] ) );
chunk.atomicWeightRatios()[2] ) );
chunk.numberAtoms().size() );
chunk.numberAtoms()[0] ) );
chunk.numberAtoms()[1] ) );
chunk.numberAtoms()[2] ) );
chunk.analyticalFunctionTypes().size() );
chunk.analyticalFunctionTypes()[0] );
chunk.analyticalFunctionTypes()[1] );
*/


  
  //std::cout << constants << std::endl;
      //auto principalTemp = mt4.principalEffectiveTemperature(); // Teff for principal scatterer
      //auto secondaryTemps = mt4.secondaryEffectiveTemperatures(); // Teff for every secondary scaterer (if any)

      //std::cout << typeid(law).name() << std::endl;

      //std::cout << law << std::endl;


  return 0;
}

