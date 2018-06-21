#include "catch.hpp"
#include "temp_eff/temp_eff_util/gateff.h"


TEST_CASE( "gateff" ){
  std::vector<double> matVals { 1002, 1004, 1064, 1099, 1065, 1114, 1096, 1097,
    1095, 1167 };

  GIVEN( "one material temperature" ){
    std::vector<double> temp(1), eftemp(1), effTemp(10);
    WHEN( "we loop through all ENDF thermal scattering materials" ){
      AND_WHEN( "material temperture is low" ){
        temp[0] = 296;
        THEN( "effective temperature is correctly pulled from table" ){
          effTemp = { 1396.8, 940.91, 405.64, 596.4, 713.39, 1222, 317.27, 
          806.79, 1165.9, 296 };
          for ( size_t i = 0; i < matVals.size(); ++i ){
            gateff( temp, eftemp, matVals[i] );
            REQUIRE( effTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
            eftemp = { 0 };
          }
        } // THEN
      } // AND WHEN
      AND_WHEN( "material temperature is high" ){
        temp = {650};
        THEN( "material temperature is used as the effective temperature" ){
          for ( size_t i = 0; i < matVals.size(); ++i ){
            gateff( temp, eftemp, matVals[i] );
            REQUIRE( temp[0] == Approx( eftemp[0] ).epsilon(1e-6) );
            eftemp = { 0 };
          }
        } // THEN
      } // AND WHEN
    } // WHEN
  } // GIVEN
  GIVEN( "multiple temperature" ){
    std::vector<double> temp { 200, 300, 400, 500 }, eftemp(4);
    WHEN( "we loop through all ENDF thermal scattering materials" ){
      std::vector<std::tuple<double,double,double,double>> effTemp { 
        { 200, 1396.8, 1427.4, 1464.1 }, { 200, 940.91, 982.93, 1030.9 },
        { 200, 405.64, 484.22, 568.53 }, { 200, 596.4,  643.9,  704.6 },
        { 200, 713.39, 754.68, 806.67 }, { 200, 1222,   400,    500 },
        { 200, 317.27, 416.29, 513.22 }, { 200, 806.79, 829.98, 868.44 },
        { 200, 1165.9, 1191.4, 1226   }, { 200, 300,    400,    500 }
      };
      for ( size_t i = 0; i < matVals.size(); ++i ){
        gateff( temp, eftemp, matVals[i] );
        REQUIRE( std::get<0>(effTemp[i]) == Approx( eftemp[0] ).epsilon(1e-6) );
        REQUIRE( std::get<1>(effTemp[i]) == Approx( eftemp[1] ).epsilon(1e-6) );
        REQUIRE( std::get<2>(effTemp[i]) == Approx( eftemp[2] ).epsilon(1e-6) );
        REQUIRE( std::get<3>(effTemp[i]) == Approx( eftemp[3] ).epsilon(1e-6) );
        eftemp = { 0, 0, 0, 0 };
      }
    } // WHEN
  } // GIVEN
} // TEST CASE


