#include "catch.hpp"
#include "iel/incoherentElastic.h"


TEST_CASE( "get equiprobable scattering angles" ){
  GIVEN( "an initial energy (which is also the final energy here)" ){
    double E = 1e-5, debyeWaller = 0.1, sigma_b = 2.5;
    int numAtoms = 1, numAngles = 8;

    getIncohElasticEquiprobableAngles(E,debyeWaller,sigma_b,numAtoms,numAngles);


  } // GIVEN
} // TEST CASE
