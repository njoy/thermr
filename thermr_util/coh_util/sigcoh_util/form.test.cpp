#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "form.h"
#include <tuple>


TEST_CASE( "form" ){
  int lat, l1, l2, l3;
  GIVEN( "the material of interest is graphite" ){
    lat = 1;
    THEN( "computed form factors are correct" ){
      std::vector<std::tuple<int,int,int>> inputs { { 0, 0, 0 }, { 1, 0, 0 }, 
        { 0, 1, 0 }, { 0, 0, 1 }, { 1, 1, 0 }, { 1, 0, 1 }, { 0, 1, 1 }, 
        { 1, 1, 1 }, {-1, 0, 0 }, { 0,-1, 0 }, { 0, 0,-1 }, {-1, 1, 0 },
        { 1, 0,-1 }, { 0,-1, 1 }, {-1, 1,-1 }, { 2, 3,-4 }, { 3, 1, 0 }, 
        { 0,-2,-3 }, { 1,-2, 1 } };
      std::vector<double> correctVals { 4, 0.25, 0.25, 0, 4, 0.75, 0.75, 0, 
        0.25, 0.25, 0, 0.25, 0.75, 0.75, 0.75, 0.25, 0.25, 0.75, 0};
      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( correctVals[i] == Approx( form( lat, std::get<0>(inputs[i]), 
          std::get<1>(inputs[i]), std::get<2>(inputs[i]) ) ).epsilon(1e-6) );
      }
      
    } // THEN
  } // GIVEN
} // TEST CASE
