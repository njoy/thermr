#include "catch.hpp"
#include "iel/incoherentElastic.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"

TEST_CASE( "get equiprobable scattering angles" ){
  GIVEN( "an initial energy (which is also the final energy here)" ){
    double E = 1e-2, debyeWaller = 10.0, sigma_b = 2.5, xs;
    int numAtoms = 1, numAngles = 8;
    std::cout.precision(15);
    std::vector<double> cosines ( numAngles, 0.0 );
    std::vector<double> correctCosines(numAngles);


    E = 1e-5;
    correctCosines = {-0.8749771,-0.6249396,-0.3749146,-0.1249021, 
                       0.1250979, 0.3750854, 0.6250604, 0.8750229 };
    getIncohElasticDataSingleEnergy(E,debyeWaller,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );



    E = 1e-4;
    correctCosines = {-0.8747706,-0.6243953,-0.3741454,-0.1240207, 
                       0.1259790, 0.3758537, 0.6256037, 0.8752289 };
    getIncohElasticDataSingleEnergy(E,debyeWaller,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );

    E = 1e-3;
    correctCosines = {-0.8726826,-0.6189088,-0.3664165,-0.1151928, 
                       0.1347749, 0.3834992, 0.6309923, 0.8772663 };
    getIncohElasticDataSingleEnergy(E,debyeWaller,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );

    E = 1e-2;
    correctCosines = {-0.8493610,-0.5596117,-0.2857415,-2.6099821E-2, 
                       0.2207195, 0.4559238, 0.6805574, 0.8955298 };
    getIncohElasticDataSingleEnergy(E,debyeWaller,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );

    E = 1e-1;
    correctCosines = {-0.3270719, 0.1933194, 0.4351759, 0.5967253, 
                       0.7184233, 0.8161444, 0.8978172, 0.9679844 };
    getIncohElasticDataSingleEnergy(E,debyeWaller,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );



  } // GIVEN
} // TEST CASE



TEST_CASE( "incoherent elastic average scattering xs" ){
    double debyeWaller = 10.0, sigma_b = 2.5, xs;
    int numAtoms = 1;

    std::vector<double> 
      energyVec { 1e-5, 1e-4, 1e-3, 1e-2, 1e-1 },
      correctXS { 2.4995001, 2.4950067, 2.4506601, 2.0604997, 0.6135527 };

    for ( size_t i = 0; i < energyVec.size(); ++i ){
      REQUIRE( correctXS[i] == Approx(getIncohElasticXS(sigma_b,debyeWaller,
                                      energyVec[i],numAtoms)).epsilon(1e-6) );
    }
    

} // TEST CASE
