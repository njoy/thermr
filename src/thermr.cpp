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

using namespace njoy::ENDFtk;

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
auto thermr( int matde, int matdp, int nbin, int iinc, int icoh, int iform,
  int natom, int mtref, Range temps, Float tol, Float emax,
  MF7 leapr_MF7 ){

  double kb = 8.6173303e-5;
  std::cout.precision(15);

  //njoy::ENDFtk::section::Type<7,2> leapr_MT2 = leapr_MF7.MT(2_c);
  njoy::ENDFtk::section::Type<7,4> leapr_MT4 = leapr_MF7.MT(4_c);
  //std::cout << leapr_MT2.ZA() << std::endl;
  //std::cout << leapr_MT4.ZA() << std::endl;
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
  
  //std::cout << constants.LLN() << std::endl;
  //std::cout << constants.sabStorageType() << std::endl;
  //std::cout << constants.NI() << std::endl;
  //std::cout << constants.numberConstants() << std::endl;
  //std::cout << constants.NS() << std::endl;
  //std::cout << constants.numberNonPrincipalScatterers() << std::endl;

  //std::string buffer;
  //auto output = std::back_inserter(buffer);
  //leapr_MT4.print(output,27);
  //std::cout << buffer << std::endl;
  

  //std::cout << (alphas|ranges::view::all) << std::endl;
  //std::cout << (betas|ranges::view::all) << std::endl;



  //int iverf = 6;
  // if endf6 format, then iverf = 4 or 5???
  

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
    if (0 > 1){
    std::cout << (sab|ranges::view::all) << std::endl;
         std::cout << awr << std::endl;
    }
    //std::cout << table.betas()[0].beta() << std::endl;
    //std::cout << table.betas()[0].thermalScatteringValues()[itemp][0] << std::endl;
    //std::cout << table.betas()[0].thermalScatteringValues()[itemp][1] << std::endl;
    //std::cout << table.betas()[1].thermalScatteringValues()[itemp][0] << std::endl;
    //std::cout << table.betas()[1].thermalScatteringValues()[itemp][1] << std::endl;


    auto temp = temps[itemp];
    auto tev  = temp*kb;
    auto teff = eftemp[itemp];

    
    // compute incoherent cross sections
    if (iinc != 0){
       //double teff  = 0;
       //double teff2 = 0;

       if (iform == 0){
         // E E' mu
         //std::cout << "E E' mu" << std::endl;
         //std::cout << tev << std::endl;
         auto effectiveTemp = leapr_MT4.principalEffectiveTemperature();
         teff = effectiveTemp.effectiveTemperatures()[0]*kb;

         //std::cout << "THIS TEFF    " << teff << std::endl;
         auto out = e_ep_mu( egrid, tev, tol, lat,  iinc, lasym, alphas, betas, 
                             sab, awr, boundCrossSections, teff, nbin, temp );

         auto incidentEnergies = std::get<0>(out);
         //std::cout << "INCIDENT ENERGY    " << (incidentEnergies |ranges::view::all) << std::endl;
         auto totalSCR     = std::get<1>(out);
         auto totalOutput  = std::get<2>(out);
         //std::cout << outputEnergy[0] << "   " << outputEnergy[1] << std::endl;
         //std::cout << totalOutput[0][0] << "   " << totalOutput[0][1] << std::endl;
         //std::cout << totalSCR[0][0] << "   " << totalSCR[0][1] << std::endl;
         int n2 = nbin+2;
         //std::cout << "N2    " << n2 << std::endl;
         //std::cout << (totalSCR[0]|ranges::view::all) << std::endl;
         //std::cout << (chunkSCR[0]|ranges::view::all) << std::endl;
         //std::cout << (chunkSCR[1]|ranges::view::all) << std::endl;
         //std::cout << (chunkSCR[2]|ranges::view::all) << std::endl;
         //std::cout << (chunkSCR[3]|ranges::view::all) << std::endl;

         auto firstSCR = totalSCR[0];
         ThermalScatteringData chunk( incidentEnergies[0], n2, std::move(firstSCR) );
         std::cout << chunk.LANG() << std::endl;
         std::cout << chunk.LTT() << std::endl;
         std::cout << chunk.energy() << std::endl;
         std::cout << chunk.NW() << std::endl;
         std::cout << chunk.N2() << std::endl;
         std::cout << std::endl;
         std::cout << std::endl;
         //for ( size_t i = 0; i < incidentEnergies.size(); ++i ){
           //auto chunkSCR = totalSCR[i] | ranges::view::chunk(n2);
           //ThermalScatteringData chunk( incidentEnergies[i], n2, std::move(chunkSCR) );
         //}      
         /*

         */


         /*
      double energy = 1e-5; // INCOMING ENERGY
      int n2 = 6;                  // Eout, pp/totalpp?, mu1, mu2, mu3 Eout pp/totalpp? mu1...
      std::vector< double > data = { 0., 0., 0., 0., 0., 0.,
                                     9.999999e-6, 9.477167e+1, -5.379121e-1,
                                     0.21062848, 0.70490082, 9.552579e-1,
                                     1.265100e-1, 0., 0., 0., 0., 0. };
                                     */

         //for (int i = 0; i < 20; ++i){
         //  std::cout << " ----  " << totalSCR[0][i] << std::endl;
         //}
       }
       else {
         // E mu E' 
       }
       //if (sz2 > 0) {teff2=eftmp2[itemp];}
       //if (iinc == 2) temp=t;
       //calcem(temp,itemp,iold,inew,np,nex)
    }



    if (icoh > 0 and icoh <= 10){
      std::cout << "COH" << std::endl;
      coh( temp, lat, emax, natom, egrid, tol);
    }
    else if (icoh > 10){
      std::cout << "INCOH" << std::endl;
    }

  }

  
  return;
  std::cout << za << awr << lat << lasym << std::endl;
  std::cout << matde+matdp+nbin+iinc+icoh+iform+natom+mtref<< std::endl;
  std::cout << tol + emax + temps[0] << std::endl;
  //std::cout << tev << std::endl;
  std::cout << alphas.size()<< std::endl;
  std::cout << betas.size()<< std::endl;
  std::cout << kb << std::endl;

}




/*

int THERMR( int mat, std::string fileName, int iform, int iinc, int nbin,
  std::vector<double> temperatures, double tol, double emax, int natom){


   //std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5,
   //   1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 0.000253, 0.000297, 0.000350, 0.00042, 
   //   0.000506, 0.000615, 0.00075, 0.00087, 0.001012, 0.00123, 0.0015, 
   //   0.0018, 0.00203, 0.002277, 0.0026, 0.003, 0.0035, 0.004048, 0.0045, 
   //   0.005, 0.0056, 0.006325, 0.0072, 0.0081, 0.009108, 0.01, 0.01063, 
   //   0.0115, 0.012397, 0.0133, 0.01417, 0.015, 0.016192, 0.0182, 0.0199, 
   //   0.020493, 0.0215, 0.0228, 0.0253, 0.028, 0.030613, 0.0338, 0.0365, 
   //   0.0395, 0.042757, 0.0465, 0.050, 0.056925, 0.0625, 0.069, 0.075, 
   //   0.081972, 0.09, 0.096, 0.1035, 0.111573, 0.120, 0.128, 0.1355, 0.145728, 
   //   0.160, 0.172, 0.184437, 0.20, 0.2277, 0.2510392, 0.2705304, 0.2907501, 
   //   0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 0.5032575, 0.56, 
   //   0.625, 0.70, 0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 1.70, 1.855,
   //   2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 4.90, 5.35,
   //   5.85, 6.40, 7.00, 7.65, 8.40, 9.15, 9.85, 10.00 };


   double kb = 8.6173303e-5;
   //double tev = T*kb;


  //std::string contents = njoy::utility::slurpFileToMemory(fileName);
  //njoy::ENDFtk::syntaxTree::Tape< std::string > tape = 
  //  njoy::ENDFtk::syntaxTree::Tape< std::string >( contents );

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


  //std::vector<double> sab;
  //std::vector<std::vector<double>> chunkySAB;
  //for (size_t i = 0; i < (unsigned) nbeta; ++i){
  //  betas[i] = table.betas()[i].beta();
  //  chunkySAB.push_back(table.betas()[i].thermalScatteringValues()[0]);
 // }

  //for (size_t j =0; j < chunkySAB[0].size(); ++j){
  //  for (size_t i = 0; i < chunkySAB.size(); ++i){
  //    sab.push_back(std::log(chunkySAB[i][j]));
  //  }
  //}

  //int a, b;
  //a = 10; 
  //b = 5;
  //std::cout << sab[a*nbeta+b] << std::endl;


  //return 0;
  

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

//chunk.LLN() );
//chunk.sabStorageType() );
//chunk.NI() );
//chunk.numberConstants() );
//chunk.NS() );
//chunk.numberNonPrincipalScatterers() );
//chunk.epsilon() ) );
//chunk.upperEnergyLimit() ) );
//chunk.totalFreeCrossSections().size() );
//chunk.totalFreeCrossSections()[0] ) );
//chunk.totalFreeCrossSections()[1] ) );
//chunk.totalFreeCrossSections()[2] ) );
//chunk.atomicWeightRatios().size() );
//chunk.atomicWeightRatios()[0] ) );
//chunk.atomicWeightRatios()[1] ) );
//chunk.atomicWeightRatios()[2] ) );
//chunk.numberAtoms().size() );
//chunk.numberAtoms()[0] ) );
//chunk.numberAtoms()[1] ) );
//chunk.numberAtoms()[2] ) );
//chunk.analyticalFunctionTypes().size() );
//chunk.analyticalFunctionTypes()[0] );
//chunk.analyticalFunctionTypes()[1] );


  
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
*/


