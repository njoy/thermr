#include "ENDFtk.hpp"
#include <typeinfo>       // operator typeid

using namespace njoy::ENDFtk;

using ScatteringLaw = section::Type< 7, 4 >::ScatteringLaw;
using ScatteringLawConstants = section::Type< 7, 4 >::ScatteringLawConstants;
using AnalyticalFunctions = section::Type< 7, 4 >::AnalyticalFunctions;
using Tabulated = section::Type< 7, 4 >::Tabulated;
using ScatteringFunction = section::Type< 7, 4 >::Tabulated::ScatteringFunction;
using EffectiveTemperature = section::Type< 7, 4 >::EffectiveTemperature;

using namespace njoy::ENDFtk;
int thermr( int mat, std::string fileName ){


    std::string contents = njoy::utility::slurpFileToMemory(fileName);//"/Users/ameliajo/thermr/src/test/tape24");
     njoy::ENDFtk::syntaxTree::Tape< std::string > tape = 
         njoy::ENDFtk::syntaxTree::Tape< std::string >( contents );

     //int mat = 101;
     decltype(auto) materials = tape.MAT( mat ); // MAT( ... ) returns all materials with that nuber - a requirement due to multitemperature tapes
     decltype(auto) material = materials.front(); // this is the first material, the one you need
     decltype(auto) file7 = material.MF(7).parse< 7 >();
     

     decltype(auto) mt4 = file7.section( 4_c ); // MT( 4_c ) would achieve the same thing

      // mt4 is now a fully parser MF7 MT4 from which you can now request the data in it:

      auto LAT = mt4.temperatureOption(); // LAT()
      std::cout << LAT << std::endl;
      //auto LASYM = mt4.symmetryOption(); // LASYM
      //auto constants = mt4.constants(); // the ScatteringLawConstants instance inside mt4
      //auto principalTemp = mt4.principalEffectiveTemperature(); // Teff for principal scatterer
      //auto secondaryTemps = mt4.secondaryEffectiveTemperatures(); // Teff for every secondary scaterer (if any)
      auto law = mt4.scatteringLaw(); // a variant of either a tabulated law or an analytical law

      std::cout << typeid(law).name() << std::endl;

      //std::cout << law << std::endl;


  return 0;
}

