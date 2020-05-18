#include "ENDFtk.hpp"

using namespace njoy::ENDFtk;
using ScatteringLaw = section::Type< 7, 4 >::ScatteringLaw;
using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
using AnalyticalFunctions = section::Type< 7, 4 >::AnalyticalFunctions;
using Tabulated = section::Type< 7, 4 >::Tabulated;
using ScatteringFunction = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
using EffectiveTemperature = section::Type< 7, 4 >::EffectiveTemperature;




template <typename Tape>
auto calcem( int iinc, int mat, int nbin, int natom, double tol, double temp, 
  Tape mt4, double effectiveTemperature){//, int itemp  ){

  double kb = 8.6173303e-5;
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


  if (iinc == 1){
      std::cout << "some free gas stuff" << std::endl;
  }

  /*decltype(auto) materials = tape.MAT( mat ); 
  decltype(auto) material = materials.front(); 
  decltype(auto) file7 = material.MF(7).parse< 7 >();
  decltype(auto) mt4 = file7.section( 4_c );
  */

  auto za = mt4.ZA();
  auto awr = mt4.AWR();
  //std::cout << awr << std::endl;
  auto lat = mt4.temperatureOption(); // LAT()
  auto lasym= mt4.symmetryOption(); // LASYM
  auto constants = mt4.constants(); // the ScatteringLawConstants instance inside mt4
  //std::cout << mat << std::endl;

  std::vector<double> totalFreeXs = constants.totalFreeCrossSections();
  //auto atomicWeightRatios = constants.atomicWeightRatios();

  double smz = totalFreeXs[0]/natom;
  //std::cout << "--------" << std::endl;
  //std::cout << smz << std::endl;
  //std::cout << za << std::endl;
  //std::cout << "--------" << std::endl;
  double sb = smz*std::pow(((awr+1)/awr),2);
  std::cout << "SB  " << sb << std::endl;

  double tev = temp * kb;

  auto law = mt4.scatteringLaw(); // a variant of either a tabulated law or an analytical law
  auto table = std::get< Tabulated >( law );
  auto value = table.betas()[0];
  auto nbeta = table.numberBetas();
  std::vector<double> betas (nbeta);
  std::vector<double> sab;//(nbeta);
  for (size_t i = 0; i < (unsigned) nbeta; ++i){
    betas[i] = table.betas()[i].beta();
    //std::vector<double> sabVal = table.betas()[i].thermalScatteringValues()[0];
    for (auto& sabVal : table.betas()[i].thermalScatteringValues()[0]){
        sab.push_back(sabVal);
    }
    //sab.push_back(table.betas()[i].thermalScatteringValues()[0]);
    //sab[i] = table.betas()[i].thermalScatteringValues()[0];

  }
  std::vector<double> alphas = table.betas()[0].alphas();


  std::cout << (egrid|ranges::view::all) << std::endl;
  std::cout << tev << std::endl;
  std::cout << tol << std::endl;
  std::cout << lat << std::endl;
  std::cout << iinc << std::endl;
  std::cout << lasym << std::endl;
  std::cout << alphas[0] << "    " << alphas[10] << "     " << alphas.size() << std::endl;
  std::cout << betas[0] << "    " << betas[10] << "     " << betas.size() << std::endl;
  std::cout << sab[0] << "    " << sab.size() << std::endl;
  std::cout << za  << std::endl;
  std::cout << (totalFreeXs|ranges::view::all) << std::endl;
  std::cout << effectiveTemperature  << std::endl;
  std::cout <<  nbin << std::endl;
  std::cout <<  temp<< std::endl;

  //e_ep_mu(egrid, tev, tol, lat, iinc, lasym, alphas, betas, sab, za, totalFreeXs, effectiveTemperature, nbin, temp);
  return;

  std::cout << effectiveTemperature << std::endl;
  std::cout << tol << std::endl;
  std::cout << za << std::endl;
  std::cout << lat << std::endl;
  std::cout << lasym << std::endl;
  std::cout << tev << std::endl;
  std::cout << nbin << std::endl;
  std::cout << betas.size() << std::endl;
  std::cout << alphas.size() << std::endl;
  std::cout << sab.size() << std::endl;
  std::cout << mat << std::endl;


  //std::cout << (effectiveTemperatures.size()) << std::endl;
  //std::cout << (secondaryEffectiveTemperatures.size()) << std::endl;
}


