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

    getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,cosines);
    std::vector<double> correctCosines(numAngles);


    E = 1e-5;
    correctCosines = {-0.8749771,-0.6249396,-0.3749146,-0.1249021, 
                       0.1250979, 0.3750854, 0.6250604, 0.8750229 };
    xs = getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );
    REQUIRE( 2.4995001 == Approx(xs).epsilon(1e-6) );


    E = 1e-4;
    correctCosines = {-0.8747706,-0.6243953,-0.3741454,-0.1240207, 
                       0.1259790, 0.3758537, 0.6256037, 0.8752289 };
    xs = getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );
    REQUIRE( 2.4950067 == Approx(xs).epsilon(1e-6) );

    E = 1e-3;
    correctCosines = {-0.8726826,-0.6189088,-0.3664165,-0.1151928, 
                       0.1347749, 0.3834992, 0.6309923, 0.8772663 };
    xs = getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );
    REQUIRE( 2.4506601 == Approx(xs).epsilon(1e-6) );

    E = 1e-2;
    correctCosines = {-0.8493610,-0.5596117,-0.2857415,-2.6099821E-2, 
                       0.2207195, 0.4559238, 0.6805574, 0.8955298 };
    xs = getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );
    REQUIRE( 2.0604997 == Approx(xs).epsilon(1e-6) );

    E = 1e-1;
    correctCosines = {-0.3270719, 0.1933194, 0.4351759, 0.5967253, 
                       0.7184233, 0.8161444, 0.8978172, 0.9679844 };
    xs = getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,cosines);
    REQUIRE( ranges::equal( cosines, correctCosines, equal ) );
    REQUIRE(  0.6135527 == Approx(xs).epsilon(1e-6) );




  } // GIVEN
} // TEST CASE
