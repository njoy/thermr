#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "gateff.h"


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





TEST_CASE( "gatef2" ){
  int mat;
  GIVEN( "only one material temperature is considered" ){
    std::vector<double> temp(1), eftemp(1), tempVals(8), effTemp(8);

    WHEN( "material considered is Benzine" ){
      mat = 1095;
      tempVals = { 200, 250, 296, 300, 315, 350, 600, 743 };
      effTemp = { 200, 250, 685.54, 685.54, 315, 712.02, 866.63, 743 };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp[0] = tempVals[i];
          gatef2( temp, eftemp, mat );
          REQUIRE( effTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
          eftemp = { 0 };
        } // for
      } // THEN
    } // WHEN

    WHEN( "material considered is BeO" ){
      mat = 1099;
      tempVals = { 200, 250, 296, 300, 315, 350, 600, 743 };
      effTemp = { 200, 250, 427.8, 427.8, 315, 350, 671.3, 743 };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp[0] = tempVals[i];
          gatef2( temp, eftemp, mat );
          REQUIRE( effTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
          eftemp = { 0 };
        } // for
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "multiple material temperatures are considered" ){
    std::vector<double> temp(4), eftemp(4);
    std::vector<std::tuple<double,double,double,double>> tempVals(4), effTemp(4);

    WHEN( "material considered is Benzine" ){
      mat = 1095;
      tempVals = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 296.00, 300.00, 301.00, 405.00 }, { 300.00, 315.00, 320.00, 345.00 },
        { 400.00, 455.00, 520.00, 645.00 } };
      effTemp  = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 685.54, 685.54, 685.54, 738.97 }, { 685.54, 315.00, 320.00, 712.02 },
        { 738.97, 768.10, 520.00, 645.00 } };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp = { std::get<0>(tempVals[i]), std::get<1>(tempVals[i]), 
                   std::get<2>(tempVals[i]), std::get<3>(tempVals[i]) };
          gatef2( temp, eftemp, mat );
          REQUIRE( std::get<0>(effTemp[i]) == Approx( eftemp[0] ).epsilon(1e-6) );
          REQUIRE( std::get<1>(effTemp[i]) == Approx( eftemp[1] ).epsilon(1e-6) );
          REQUIRE( std::get<2>(effTemp[i]) == Approx( eftemp[2] ).epsilon(1e-6) );
          REQUIRE( std::get<3>(effTemp[i]) == Approx( eftemp[3] ).epsilon(1e-6) );
          eftemp = { 0, 0, 0, 0 };
        } // end for
      } // THEN
    } // WHEN

    WHEN( "material considered is BeO" ){
      mat = 1099;
      tempVals = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 296.00, 300.00, 301.00, 405.00 }, { 300.00, 315.00, 320.00, 345.00 },
        { 400.00, 455.00, 520.00, 645.00 } };
      effTemp  = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 427.80, 427.80, 427.80, 502.80 }, { 427.80, 315.00, 320.00, 345.00 },
        { 502.80, 455.00, 520.00, 645.00 } };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp = { std::get<0>(tempVals[i]), std::get<1>(tempVals[i]), 
                   std::get<2>(tempVals[i]), std::get<3>(tempVals[i]) };
          gatef2( temp, eftemp, mat );
          REQUIRE( std::get<0>(effTemp[i]) == Approx( eftemp[0] ).epsilon(1e-6) );
          REQUIRE( std::get<1>(effTemp[i]) == Approx( eftemp[1] ).epsilon(1e-6) );
          REQUIRE( std::get<2>(effTemp[i]) == Approx( eftemp[2] ).epsilon(1e-6) );
          REQUIRE( std::get<3>(effTemp[i]) == Approx( eftemp[3] ).epsilon(1e-6) );
          eftemp = { 0, 0, 0, 0 };
        } // end for
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
