#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "160.h"



TEST_CASE( "tausq" ){
  GIVEN( "inputs" ){
    std::vector<std::tuple<int,int,int,double,double>> inputs
      { {0,0,0,2,4}, {1,0,0,2,4}, {0,1,0,2,4}, {0,0,1,2,4}, {1,1,0,2,4}, 
	{1,0,1,2,4}, {0,1,1,2,4}, {1,1,1,2,4}, {1,2,3,4,5}, {5,3,6,4,5},
        {8,7,9,.1,.2} };
    std::vector<double> output { 0, 78.956835, 78.956835, 157.91367041, 
	236.870505, 236.870505, 236.870505, 394.78417604, 2881.924485, 
	14843.885019, 1306.735642 };
    
    for ( auto i = 0; i < output.size(); ++i ){
      REQUIRE( tausq( std::get<0>(inputs[i]), std::get<1>(inputs[i]),
        std::get<2>(inputs[i]), std::get<3>(inputs[i]), 
        std::get<4>(inputs[i]) ) == Approx(output[i]).epsilon(1e-6) );
    }
  } // GIVEN
} // TEST CASE


TEST_CASE( "computeCrossSections" ){
  GIVEN( "" ){
      
    THEN( "" ){
      REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
