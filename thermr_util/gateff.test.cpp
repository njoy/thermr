#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "gateff.h"


TEST_CASE( "gateff" ){
  std::vector<double> matVals { 1002 , 1004 , 1064 , 1099 , 1065 , 1114 , 1096 , 1097 , 1095 , 1167 };
  GIVEN( "one temperature" ){
    std::vector<double> temp = {296}, eftemp(1);

    WHEN( "we loop through all ENDF thermal scattering materials" ){
      std::vector<double> effectiveTemp { 1396.8, 940.91, 405.64, 596.4, 713.39, 1222,   317.27, 806.79, 1165.9, 296 };
      for ( size_t i = 0; i < matVals.size(); ++i ){
        gateff( temp, eftemp, matVals[i] );
        REQUIRE( effectiveTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
        eftemp[0] = 0;
      }
    } // WHEN

    WHEN( "we loop again through all ENDF thermal scattering materials" ){
      temp = {650};
      for ( size_t i = 0; i < matVals.size(); ++i ){
        gateff( temp, eftemp, matVals[i] );
        REQUIRE( temp[0] == Approx( eftemp[0] ).epsilon(1e-6) );
        eftemp[0] = 0;
      }
    } // WHEN
  } // GIVEN

  GIVEN( "multiple temperature" ){
    std::vector<double> temp { 295, 300, 400, 500 }, eftemp(4);

    WHEN( "we loop through all ENDF thermal scattering materials" ){
      std::vector<std::tuple<double,double,double,double>> effectiveTemp { 
        {1411.6, 1396.8, 1427.4, 1464.1 } };
      for ( size_t i = 0; i < matVals.size(); ++i ){
        gateff( temp, eftemp, matVals[i] );
        std::cout << eftemp[0] << std::endl;
        //REQUIRE( std::get<0>(effectiveTemp[i]) == Approx( eftemp[0] ).epsilon(1e-6) );
        //REQUIRE( effectiveTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
        //REQUIRE( effectiveTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
        //REQUIRE( effectiveTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
        break;
        eftemp[0] = 0;
      }
    } // WHEN

  } // GIVEN



} // TEST CASE
