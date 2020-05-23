#include "ENDFtk.hpp"
#include <range/v3/all.hpp>
//#include "calcem/calcem_util/e_mu_ep.h"
#include "calcem/calcem_util/e_ep_mu.h"
#include "calcem/calcem.h"

using namespace njoy::ENDFtk;
using ScatteringLaw = section::Type< 7, 4 >::ScatteringLaw;
using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
using AnalyticalFunctions = section::Type< 7, 4 >::AnalyticalFunctions;
using Tabulated = section::Type< 7, 4 >::Tabulated;
using ScatteringFunction = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
using EffectiveTemperature = section::Type< 7, 4 >::EffectiveTemperature;

using namespace njoy::ENDFtk;
int thermr( int mat, std::string fileName, int iform, int iinc, int nbin,
  std::vector<double> temperatures, double tol, double emax, int natom){


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


   double kb = 8.6173303e-5;
   //double tev = T*kb;


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
  auto ZA = mt4.ZA();
  auto LAT = mt4.temperatureOption(); // LAT()
  auto LASYM = mt4.symmetryOption(); // LASYM
  auto constants = mt4.constants(); // the ScatteringLawConstants instance inside mt4


  auto temp = mt4.principalEffectiveTemperature();
  std::vector<double> effectiveTemperatures = temp.effectiveTemperatures();
  std::cout << "EFFECTIVE TEMPERATURES    " <<  (effectiveTemperatures|ranges::view::all) << std::endl;
  auto moderatorTemperatures = temp.moderatorTemperatures();
  std::cout << "MODERATOR TEMPERATURES    " <<  (moderatorTemperatures|ranges::view::all) << std::endl;
  //auto secondaryTemp = mt4.secondaryEffectiveTemperatures()[0].value();
  //auto secondaryTemp = mt4.secondaryEffectiveTemperatures()[0].value();
  //std::cout << secondaryTemp.numberTemperatures() << std::endl;
  //auto secondaryEffectiveTemperatures = secondaryTemp.effectiveTemperatures();
  //std::cout << "SEC. EFFECTIVE TEMPERATURES    " <<  (secondaryEffectiveTemperatures|ranges::view::all) << std::endl;
  //std::vector<double> secondaryEffectiveTemperatures {0.0};




  auto totalFreeXS = constants.totalFreeCrossSections();
  auto atomicWeightRatios = constants.atomicWeightRatios();
  auto law = mt4.scatteringLaw(); // a variant of either a tabulated law or an analytical law
  auto table = std::get< Tabulated >( law );
  auto value = table.betas()[0];
  auto nbeta = table.numberBetas();
  std::vector<double> betas (nbeta);
  for (size_t i = 0; i < (unsigned) nbeta; ++i){
    betas[i] = table.betas()[i].beta();
  }
  std::vector<double> alphas = table.betas()[0].alphas();

  //auto sab = table.betas()[0].thermalScatteringValues();

  //std::cout << value.beta() << std::endl;
  //std::cout << value.LT() << std::endl;
  //std::cout << value.temperatureDependenceFlag() << std::endl;
  //std::cout << value.NT() << std::endl;
  //std::cout << value.numberTemperatures() << std::endl;

  //std::cout << betas << std::endl;

  //std::cout << typeid(betas).name() << std::endl;
  //std::cout << (betas|ranges::view::all) << std::endl;


  /*
  std::vector<double> sab;
  std::vector<std::vector<double>> chunkySAB;
  for (size_t i = 0; i < (unsigned) nbeta; ++i){
    betas[i] = table.betas()[i].beta();
    chunkySAB.push_back(table.betas()[i].thermalScatteringValues()[0]);
  }

  for (size_t j =0; j < chunkySAB[0].size(); ++j){
    for (size_t i = 0; i < chunkySAB.size(); ++i){
      sab.push_back(std::log(chunkySAB[i][j]));
    }
  }

  int a, b;
  a = 10; 
  b = 5;
  std::cout << sab[a*nbeta+b] << std::endl;


  return 0;

  */


  

  if (iform == 0){
    // E E' mu
    calcem(iinc, mat, nbin, natom, tol, temperatures[0], mt4, effectiveTemperatures[0]);
    //auto out = //std::make_tuple(eVec,total_SCR,total_OutputData);
    //e_ep_mu( egrid, tev, tol, LAT, iinc, LASYM, alphas, betas, sab, ZA, , teff, nbin, temp );
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
  std::cout << (totalFreeXS|ranges::view::all) << std::endl;
  std::cout << kb +tol+emax << "   " << temperatures[0] << std::endl;
  std::cout << ZA << std::endl;
  std::cout << LAT << std::endl;
  std::cout << LASYM << std::endl;
  std::cout << (atomicWeightRatios|ranges::view::all) << std::endl;
  std::cout << (alphas|ranges::view::all) << std::endl;
  std::cout << ( betas|ranges::view::all) << std::endl;
  //std::cout << sab.size() << std::endl;




}


