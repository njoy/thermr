#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/e_ep_mu.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"




TEST_CASE( "E-E'-mu" ){
  REQUIRE( true );

  double T = 296.0, teff = 1397.671, teff2 = 0.0;
  e_ep_mu( T, teff, teff2 );


}

