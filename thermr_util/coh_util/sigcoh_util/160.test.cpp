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
  int l1, l2, l3, k, lat, nw;
  double c1, c2, w1, w2, w3, t2, wint;
  std::vector<double> wrk(10);

  GIVEN( "" ){
    WHEN( "wint value is zero, eliminating exponential contribution" ){
    l1 = 1, l2 = 0, l3 = 0, k = 0, lat = 1;
    c1 = 2, c2 = 4;
    w1 = 1, w2 = 2, w3 = 3;
    t2 = 5;
    wint = 0;
    nw = 10;
    do160( lat, w1, w2, w3, l1, l2, l3, k, c1, c2, t2, wrk, wint, nw );
    THEN( "" ){
      REQUIRE( 78.95683520871 == Approx( wrk[0] ).epsilon(1e-6) );
      REQUIRE( 0.168809309279 == Approx( wrk[1] ).epsilon(1e-6) );
      for ( size_t i = 2; i < wrk.size(); ++i ){ 
        REQUIRE( 0 == Approx( wrk[i] ).epsilon(1e-6) );
      }
    } // THEN
    } // WHEN
    WHEN( "wint value is significantly small to suppress exponential term" ){
      wint = 0.01;
      do160( lat, w1, w2, w3, l1, l2, l3, k, c1, c2, t2, wrk, wint, nw );
    THEN( "" ){
      REQUIRE( 78.95683520871 == Approx( wrk[0] ).epsilon(1e-6) );
      REQUIRE( 3.257395853E-3 == Approx( wrk[1] ).epsilon(1e-6) );
      for ( size_t i = 2; i < wrk.size(); ++i ){ 
        REQUIRE( 0 == Approx( wrk[i] ).epsilon(1e-6) );
      }
    } // THEN
    } // WHEN

  } // GIVEN
  
  GIVEN( "" ){
    WHEN( "wint value is zero, eliminating exponential contribution" ){
    l1 = 1, l2 = 2, l3 = 3, k = 0, lat = 1;
    c1 = 2, c2 = 4;
    w1 = 1, w2 = 2, w3 = 3;
    t2 = 5;
    wint = 0.01;
    nw = 10;
    do160( lat, w1, w2, w3, l1, l2, l3, k, c1, c2, t2, wrk, wint, nw );
    THEN( "" ){
      REQUIRE(  1658.0935393830089 == Approx( wrk[0] ).epsilon(1e-6) );
      REQUIRE(  1.5205534520116722E-66 == Approx( wrk[1] ).epsilon(1e-6) );
      for ( size_t i = 2; i < wrk.size(); ++i ){ 
        REQUIRE( 0 == Approx( wrk[i] ).epsilon(1e-6) );
      }
    } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "lat value greater than 1" ){
    WHEN( "wint value is zero, eliminating exponential contribution" ){
    l1 = 1, l2 = 2, l3 = 3, k = 0, lat = 2;
    c1 = 2, c2 = 4;
    w1 = 1, w2 = 2, w3 = 3;
    t2 = 5;
    wint = 0.01;
    nw = 10;
    do160( lat, w1, w2, w3, l1, l2, l3, k, c1, c2, t2, wrk, wint, nw );
    THEN( "" ){
      REQUIRE(  1658.0935393830089 == Approx( wrk[0] ).epsilon(1e-6) );
      for ( size_t i = 1; i < wrk.size(); ++i ){ 
        REQUIRE( 0 == Approx( wrk[i] ).epsilon(1e-6) );
      }
    } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
