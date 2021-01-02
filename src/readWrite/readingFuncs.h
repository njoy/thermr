#include "ENDFtk.hpp"
#include "lipservice.hpp"

using namespace njoy::ENDFtk;
using Tabulated = section::Type< 7, 4 >::TabulatedFunctions;


auto readLEAPRData( double& awr, int& lat, int& lasym, 
  std::vector<double>& freeXS, std::vector<double>& analyticalFunctionTypes, 
  std::vector<double>& awrVec, std::vector<double>& boundXS, 
  std::vector<double>& alphas, std::vector<double>& betas, 
  std::vector<double>& effectiveTemps, const double& pendf_awr, 
  const nlohmann::json& jsonInput ){

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
    freeXS = constants.totalFreeCrossSections();
    analyticalFunctionTypes = constants.analyticalFunctionTypes();
    awrVec = constants.atomicWeightRatios();

    alphas = table.scatteringFunctions()[0].alphas();
    betas  = table.betas();

    boundXS.resize(freeXS.size());

    boundXS[0] = freeXS[0]*std::pow((awr+1)/awr,2)/natom;
    if (analyticalFunctionTypes.size() > 0){
      if (analyticalFunctionTypes[0] == 0){
        auto awr2 = awrVec[0];
        boundXS[1] = freeXS[1]*std::pow((awr2+1)/awr2,2)/natom; 
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
    boundXS = {std::pow(((pendf_awr +1)/pendf_awr ),2)};
    analyticalFunctionTypes = {0};
    freeXS = {0};

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









